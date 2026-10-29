
from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.special import expit


def log10_reads(reads):
    """Transform nonnegative read counts, retaining zeros and missing values."""
    reads = np.asarray(reads, dtype=float)
    if np.any(reads < 0):
        raise ValueError("Read counts must be nonnegative.")
    return np.log10(reads + 1)


def inverse_log10_reads(values):
    return np.power(10.0, values) - 1


@dataclass(frozen=True)
class LinearRegression:
    x: np.ndarray
    y: np.ndarray
    slope: float
    intercept: float
    rvalue: float
    pvalue: float
    stderr: float
    intercept_stderr: float
    residual_variance: float
    confidence_level: float

    @property
    def n(self):
        return len(self.x)

    @property
    def df_resid(self):
        return self.n - 2

    @property
    def critical_value(self):
        return stats.t.ppf((1 + self.confidence_level) / 2, self.df_resid)

    def predict(self, x=None):
        """Return x, fitted mean, lower and upper pointwise mean CI."""
        x = np.linspace(self.x.min(), self.x.max(), 200) if x is None else np.asarray(x, dtype=float)
        mean = self.intercept + self.slope * x
        sxx = np.sum((self.x - self.x.mean()) ** 2)
        se = np.sqrt(self.residual_variance * (1 / self.n + (x - self.x.mean()) ** 2 / sxx))
        margin = self.critical_value * se
        return x, mean, mean - margin, mean + margin

    def summary(self):
        result = dict(n=self.n, method="OLS with intercept", confidence_level=self.confidence_level,
                      interval="pointwise mean; Student-t", df_resid=self.df_resid,
                      intercept=self.intercept, slope=self.slope, r_value=self.rvalue,
                      r_squared=self.rvalue ** 2, p_value=self.pvalue,
                      std_error_slope=self.stderr, intercept_std_error=self.intercept_stderr)
        for name, se in (("slope", self.stderr), ("intercept", self.intercept_stderr)):
            result[f"{name}_ci_lower"] = getattr(self, name) - self.critical_value * se
            result[f"{name}_ci_upper"] = getattr(self, name) + self.critical_value * se
        return result


def fit_linear_regression(x, y, confidence_level=0.95):
    """Fit finite paired observations; require n >= 3 and a varying predictor."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be one-dimensional arrays of equal length.")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must lie between 0 and 1.")
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 3:
        raise ValueError("OLS confidence intervals require at least three finite pairs.")
    if np.ptp(x) == 0:
        raise ValueError("OLS requires variation in the predictor.")
    fit = stats.linregress(x, y)
    residual_variance = np.sum((y - fit.intercept - fit.slope * x) ** 2) / (len(x) - 2)
    return LinearRegression(x, y, fit.slope, fit.intercept, fit.rvalue, fit.pvalue,
                            fit.stderr, fit.intercept_stderr, residual_variance, confidence_level)


@dataclass(frozen=True)
class LogisticRegression:
    x: np.ndarray
    y: np.ndarray
    slope: float
    intercept: float
    pvalue: float
    stderr: float
    intercept_stderr: float
    covariance: np.ndarray
    confidence_level: float

    @property
    def n(self):
        return len(self.x)

    @property
    def critical_value(self):
        return stats.norm.ppf((1 + self.confidence_level) / 2)

    @property
    def odds_ratio(self):
        """Odds ratio per one-unit increase in x."""
        return np.exp(self.slope)

    def predict(self, x=None):
        """Return x, probability and pointwise CI, computed on the logit scale."""
        x = np.linspace(self.x.min(), self.x.max(), 200) if x is None else np.asarray(x, dtype=float)
        design = np.stack([np.ones_like(x), x], axis=-1)
        logit = self.intercept + self.slope * x
        variance = np.einsum('...i,ij,...j->...', design, self.covariance, design)
        margin = self.critical_value * np.sqrt(np.maximum(variance, 0))
        return x, expit(logit), expit(logit - margin), expit(logit + margin)

    def summary(self):
        result = dict(n=self.n, method='Binary logit with intercept',
                      confidence_level=self.confidence_level,
                      interval='pointwise probability; normal on logit scale',
                      intercept=self.intercept, slope=self.slope,
                      p_value=self.pvalue, std_error_slope=self.stderr,
                      intercept_std_error=self.intercept_stderr,
                      odds_ratio=self.odds_ratio)
        for name, se in (('slope', self.stderr), ('intercept', self.intercept_stderr)):
            result[f'{name}_ci_lower'] = getattr(self, name) - self.critical_value * se
            result[f'{name}_ci_upper'] = getattr(self, name) + self.critical_value * se
        result['odds_ratio_ci_lower'] = np.exp(result['slope_ci_lower'])
        result['odds_ratio_ci_upper'] = np.exp(result['slope_ci_upper'])
        return result


def fit_logistic_regression(x, y, confidence_level=0.95):
    import warnings
    from statsmodels.discrete.discrete_model import Logit

    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError('x and y must be one-dimensional arrays of equal length.')
    if not 0 < confidence_level < 1:
        raise ValueError('confidence_level must lie between 0 and 1.')
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 3:
        raise ValueError('Logistic regression requires at least three finite pairs.')
    if not np.all(np.isin(y, [0, 1])) or len(np.unique(y)) != 2:
        raise ValueError('Logistic regression requires a binary outcome containing both 0 and 1.')
    if np.ptp(x) == 0:
        raise ValueError('Logistic regression requires variation in the predictor.')
    zeros, ones = x[y == 0], x[y == 1]
    if zeros.max() <= ones.min() or ones.max() <= zeros.min():
        raise ValueError('Complete or quasi-complete separation: no finite logit fit exists.')
    design = np.column_stack([np.ones_like(x), x])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            fit = Logit(y, design).fit(disp=False, maxiter=100)
            covariance = np.asarray(fit.cov_params())
    except (np.linalg.LinAlgError, RuntimeWarning) as error:
        raise ValueError('Logistic regression could not estimate finite coefficients and uncertainty.') from error
    if (not fit.mle_retvals['converged'] or not np.isfinite(fit.params).all()
            or not np.isfinite(covariance).all()):
        raise ValueError('Logistic regression did not converge to a finite fit.')
    return LogisticRegression(x, y, fit.params[1], fit.params[0], fit.pvalues[1],
                              fit.bse[1], fit.bse[0], covariance, confidence_level)
