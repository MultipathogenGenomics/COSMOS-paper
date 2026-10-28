from .linregress_abundance_viruses_coast import make_df_for_TG
from .mosquito_virus_corr_analysis import linear_regression_with_random_effect, clean_df, log_transform_desired_columns
import polars as pl
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from umap import UMAP
from skbio.diversity import beta_diversity
from pathlib import Path
from ...common.mosquito_ids import shorted_genus
from ...common.regression import fit_linear_regression, fit_logistic_regression
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib.gridspec as gridspec
from importlib import import_module
publication_figres = import_module("analysis.publication-figures.process")
import os

def brays_curtis_transform(data: np.ndarray, ids: list) -> np.ndarray:
    # compute Bray-Curtis dissimilarity matrix
    bc_dm = beta_diversity("braycurtis", data, ids=ids)
    # convert to a numpy array
    bc_array = bc_dm.to_data_frame().to_numpy()
    return bc_array

def hellinger_transform(data: np.ndarray) -> np.ndarray:
    # Hellinger transformation: square root of relative abundance
    data = np.asarray(data, dtype=float)
    if not np.isfinite(data).all() or (data < 0).any():
        raise ValueError('Hellinger transformation requires finite, nonnegative abundances')
    row_sums = data.sum(axis=1, keepdims=True)
    if (row_sums <= 0).any():
        raise ValueError('Hellinger transformation requires positive total abundance in each sample')
    relative_abundance = data / row_sums
    hellinger_data = np.sqrt(relative_abundance)
    return hellinger_data

def run_PCA(features, n_components=2):
    pca = PCA(n_components=n_components)
    pca.fit(features)
    return pca

def run_PCA_on_mosquito_abundance(df, n_components=2):
    # get rid of any rows with NaN values
    df = df.drop_nulls()
    # assume all columns except "sample" are for PCA
    samples = df.select("sample").to_numpy().flatten()
    species_columns = [col for col in df.columns if col != "sample"]
    features = df.select([col for col in df.columns if col != "sample"]).to_numpy()
    transformed_features = hellinger_transform(features)
    pca = run_PCA(transformed_features, n_components=n_components)
    df_PCA = pl.DataFrame(pca.transform(transformed_features), schema=[f"PC{i+1}" for i in range(n_components)])
    # Plot and model -PC1, while retaining original PC1 species loadings.
    df_PCA = df_PCA.with_columns(-pl.col("PC1"), pl.Series("sample", samples))
    # save PCA results and loadings to tables
    outpath = "/Users/afor0966/Documents/code/COSMOS-analysis/results/abundance-geospatial/PCA_loadings_with_virus/PCA/tables"
    os.makedirs(outpath, exist_ok=True)
    df_PCA.select(["sample"] + [f"PC{i+1}" for i in range(n_components)]).write_csv(os.path.join(outpath, "PCA_scores.csv"))
    pl.DataFrame(pca.components_.T, schema=[f"PC{i+1}" for i in range(n_components)]).with_columns(pl.Series("features", species_columns)).select(["features"] + [f"PC{i+1}" for i in range(n_components)]).write_csv(os.path.join(outpath, "PCA_loadings.csv"))

    return df_PCA, pca.explained_variance_ratio_, pca

def add_elevation_and_rainfall_to_data(df, filepath="data/COSMOS_metadata/cosmos_elevation_and_rainfall.csv"):
    # load elevation and rainfall data
    elevation_rainfall_df = pl.read_csv(filepath)
    # merge with the main dataframe on "sample"
    merged_df = df.join(elevation_rainfall_df, left_on="sample", right_on="sample_id", how="left")
    return merged_df

def plot_scatter(df, x_col, y_col, hue_col=None, title=None, save_path=None):
    plt.figure(figsize=(10, 8))
    if hue_col:
        sns.scatterplot(data=df.to_pandas(), x=x_col, y=y_col, hue=hue_col)
    else:
        sns.scatterplot(data=df.to_pandas(), x=x_col, y=y_col)
    if title:
        plt.title(title)
    if save_path:
        plt.savefig(save_path)
    plt.close()

def plot_multipanel_PCA(df, n_components=5, save_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/multipanel_PCA.png",
                        explained_variance_ratio=None):
    if not 2 <= n_components <= 5:
        raise ValueError('Choose 2–5 components for the four-panel figure')
    if explained_variance_ratio is not None and len(explained_variance_ratio) < n_components:
        raise ValueError('Explained variance is required for every plotted component')
    # A4 portrait width; readable at native print size, with vector PDF text.
    with plt.rc_context({'font.size': 9, 'axes.labelsize': 10, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'pdf.fonttype': 42}):
        fig, axes = plt.subplots(2, 2, figsize=(8.27, 8.0), layout='constrained')
        plot_df = df.to_pandas()
        # Input scores use -PC1 for the main figures; restore original PC1 here.
        plot_df['PC1'] = -plot_df['PC1']
        for i, ax in enumerate(axes.flat):
            if i >= n_components - 1:
                ax.set_visible(False)
                continue
            sns.scatterplot(data=plot_df, x=f'PC{i+1}', y=f'PC{i+2}', ax=ax,
                            s=30, color='grey', alpha=0.8,
                            edgecolor='white', linewidth=0.35, zorder=3)
            for component, setter in [(i, ax.set_xlabel), (i + 1, ax.set_ylabel)]:
                label = f'PC{component + 1}'
                if explained_variance_ratio is not None:
                    label += f' ({explained_variance_ratio[component]:.1%} explained variance)'
                setter(label, labelpad=7)
            # ax.set_title(f'{chr(65 + i)}  PC{i+1} vs PC{i+2}', loc='left',
            #              fontsize=11, fontweight='bold', pad=10)
            ax.axhline(0, color='0.75', linewidth=0.6, linestyle='--', zorder=1)
            ax.axvline(0, color='0.75', linewidth=0.6, linestyle='--', zorder=1)
            ax.set_box_aspect(1)
            ax.margins(0.08)
            ax.tick_params(length=3, width=0.6)
            sns.despine(ax=ax)
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300)
        fig.savefig(out.with_suffix('.pdf'))
        plt.close(fig)

def plot_PCA_loadings(pca, feature_names, save_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/PCA_loadings.png"):
    # Preserve scaled loadings: eigenvectors multiplied by sqrt(eigenvalue).
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    if len(feature_names) != loadings.shape[0]:
        raise ValueError('Feature names must match the species columns used to fit PCA')
    labels = []
    for name in feature_names:
        if name == 'ab_S-AUS_Cx_annu':
            labels.append('Cx. annulirostris (S-AUS)')
        elif name == 'ab_AUS_Cx_annu':
            labels.append('Cx. annulirostris (AUS)')
        else:
            labels.append(shorted_genus(name.removeprefix('ab_').replace('_', ' ')))
    n_components = loadings.shape[1]
    n_rows = (n_components + 1) // 3
    with plt.rc_context({'font.size': 8, 'pdf.fonttype': 42}):
        fig, axes = plt.subplots(n_rows, 3, figsize=(8.27, 3.7 * n_rows),
                                 squeeze=False, layout='constrained')
        limit = np.max(np.abs(loadings)) * 1.1
        for i, ax in enumerate(axes.flat):
            if i >= n_components:
                ax.set_visible(False)
                continue
            # Strongest species first; retain signs so direction is visible.
            order = np.argsort(-np.abs(loadings[:, i]), kind='stable')
            values = loadings[order, i]
            ax.barh(np.arange(len(labels)), values,
                    color=np.where(values >= 0, '#0072B2', '#F28E2B'))
            ax.set_yticks(np.arange(len(labels)), labels=[labels[j] for j in order],
                          fontsize=7.5, fontstyle='italic')
            ax.invert_yaxis()
            ax.axvline(0, color='0.4', linewidth=0.6)
            ax.set_xlim(-limit, limit)
            ax.set_title(f'PC{i+1} ({pca.explained_variance_ratio_[i]:.1%} explained variance)',
                         fontsize=10, loc='left', pad=8)
            ax.set_xlabel('Scaled loading')
            ax.tick_params(axis='y', length=0)
            sns.despine(ax=ax)
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300)
        fig.savefig(out.with_suffix('.pdf'))
        pl.DataFrame(loadings, schema=[f'PC{i+1}' for i in range(n_components)]).with_columns(
            pl.Series('species', feature_names)
        ).select(['species'] + [f'PC{i+1}' for i in range(n_components)]).write_csv(
            out.with_suffix('.tsv'), separator='\t',
        )
        plt.close(fig)

def associate_PCs_with_variables(df, pca_df, continuous_variables=None, categorical_variables=None, grouping_variable=None, return_coefficients=False):
    """Fit PC associations with continuous predictors z-scored within each model."""
    import warnings
    continuous_variables = continuous_variables or []
    categorical_variables = categorical_variables or []
    if set(continuous_variables) & set(categorical_variables):
        raise ValueError('A variable cannot be both continuous and categorical')
    merged_df = df.join(pca_df, on='sample', how='inner', validate='1:1').to_pandas()
    results, coefficients = [], []
    variables = [(var, True) for var in categorical_variables] + [
        (var, False) for var in continuous_variables
    ]
    for var, categorical in variables:
        for pc in [col for col in pca_df.columns if col.startswith('PC')]:
            columns = [pc, var] + ([grouping_variable] if grouping_variable else [])
            data = merged_df.dropna(subset=columns).copy()
            numeric = [pc] if categorical else [pc, var]
            data = data.loc[np.isfinite(data[numeric].to_numpy(dtype=float)).all(axis=1)]
            formula = f'{pc} ~ C({var})' if categorical else f'{pc} ~ {var}'
            row = {
                'variable': var, 'variable_type': 'categorical' if categorical else 'continuous',
                'PC': pc, 'formula': formula + (f' + (1 | {grouping_variable})' if grouping_variable else ''),
                'method': 'Linear mixed model (site random intercept)' if grouping_variable else 'OLS',
                'n': len(data), 'n_sites': int(data[grouping_variable].nunique()) if grouping_variable else None,
                'p_value': None, 'test': 'Joint Wald chi-square' if grouping_variable else 'Overall F-test',
                'r_squared': None, 'r_squared_marginal': None, 'r_squared_conditional': None,
                'converged': False, 'fit_notes': '',
                'effect_unit': 'vs reference category' if categorical else 'per 1 SD increase',
                'predictor_mean': None, 'predictor_sd': None,
            }
            if len(data) < 3 or data[var].nunique() < 2:
                row['fit_notes'] = 'Insufficient complete observations or predictor variation'
                results.append(row)
                continue
            if not categorical:
                mean, sd = float(data[var].mean()), float(data[var].std(ddof=1))
                if not np.isfinite(sd) or sd <= 0:
                    row['fit_notes'] = 'Continuous predictor has no finite, positive standard deviation'
                    results.append(row)
                    continue
                data[var] = (data[var] - mean) / sd
                row.update(predictor_mean=mean, predictor_sd=sd)
            model = None
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                try:
                    if grouping_variable:
                        if data[grouping_variable].nunique() < 2:
                            raise ValueError('Random intercept requires at least two sites')
                        model = smf.mixedlm(formula, data=data, groups=data[grouping_variable]).fit()
                        # Joint test of all predictor coefficients, excluding intercept and site variance.
                        fixed_terms = list(model.fe_params.index)
                        tested = [term for term in fixed_terms if term != 'Intercept']
                        restriction = np.zeros((len(tested), len(model.params)))
                        for i, term in enumerate(tested):
                            restriction[i, model.params.index.get_loc(term)] = 1
                        row['p_value'] = float(model.wald_test(restriction, scalar=True).pvalue)
                        fixed_variance = float(np.var(model.model.exog @ model.fe_params.to_numpy(), ddof=1))
                        site_variance = float(model.cov_re.iloc[0, 0])
                        total = fixed_variance + site_variance + float(model.scale)
                        row['r_squared_marginal'] = fixed_variance / total
                        row['r_squared_conditional'] = (fixed_variance + site_variance) / total
                        row['converged'] = bool(model.converged)
                    else:
                        model = smf.ols(formula, data=data).fit()
                        row['p_value'] = float(model.f_pvalue)
                        row['r_squared'] = float(model.rsquared)
                        row['converged'] = True
                except (ValueError, np.linalg.LinAlgError) as error:
                    row['fit_notes'] = str(error)
            notes = '; '.join(dict.fromkeys(str(item.message) for item in caught))
            row['fit_notes'] = '; '.join(note for note in [row['fit_notes'], notes] if note)
            if model is not None and row['p_value'] is not None:
                fixed = model.fe_params if grouping_variable else model.params
                terms = [term for term in fixed.index if term != 'Intercept']
                intervals = model.conf_int()
                reference = (next(str(level) for level in data[var].unique()
                                  if f'C({var})[T.{level}]' not in terms)
                             if categorical else None)
                for term in terms:
                    low, high = intervals.loc[term]
                    coefficients.append({
                        'PC': pc, 'variable': var, 'term': term, 'reference_level': reference,
                        'beta': float(fixed[term]), 'ci_lower': float(low), 'ci_upper': float(high),
                        'std_error': float(model.bse[term]), 'p_value': float(model.pvalues[term]),
                        'overall_p_value': row['p_value'],
                        'r_squared': row['r_squared_marginal'] if grouping_variable else row['r_squared'],
                        'n': len(data),
                        'n_sites': row['n_sites'], 'converged': row['converged'], 'fit_notes': row['fit_notes'],
                        'effect_unit': row['effect_unit'],
                        'predictor_mean': row['predictor_mean'], 'predictor_sd': row['predictor_sd'],
                    })
            results.append(row)
    associations = pl.DataFrame(results, schema_overrides={
        'p_value': pl.Float64, 'r_squared': pl.Float64,
        'r_squared_marginal': pl.Float64, 'r_squared_conditional': pl.Float64,
        'n_sites': pl.Int64,
        'predictor_mean': pl.Float64, 'predictor_sd': pl.Float64,
    }).sort('p_value', nulls_last=True)
    # One overall test per variable and PC; pool only PC1 and PC2 for BH FDR.
    eligible = associations.filter(
        pl.col('PC').is_in(['PC1', 'PC2']) & pl.col('p_value').is_finite()
    )
    adjusted = multipletests(eligible['p_value'].to_numpy(), method='fdr_bh')[1] if eligible.height else []
    correction = eligible.select('PC', 'variable').with_columns(
        pl.Series('p_value_fdr', adjusted, dtype=pl.Float64)
    )
    associations = associations.join(correction, on=['PC', 'variable'], how='left')
    if return_coefficients:
        coefficient_table = pl.DataFrame(coefficients, schema_overrides={
            'reference_level': pl.String, 'predictor_mean': pl.Float64, 'predictor_sd': pl.Float64,
        })
        coefficient_table = coefficient_table.join(correction, on=['PC', 'variable'], how='left')
        return associations, coefficient_table
    return associations

def plot_virus_positives_on_PCA_with_colours(df, pca_scores, fig, ax, colour_dict=None):
    # plot PCA scores with points coloured by virus presence/absence
    if colour_dict is None:
        return
    for virus, colour in colour_dict.items():
        virus_col = f"{virus}_detected"
        virus_samples = df.filter(pl.col(virus_col) == 1).select('sample').to_numpy().flatten()
        virus_pca_scores = pca_scores.filter(pl.col('sample').is_in(virus_samples)).to_pandas()
        ax.scatter(virus_pca_scores['PC1'], virus_pca_scores['PC2'], label=virus, color=colour, s=30, alpha=0.8, edgecolor='white', linewidth=0.35, zorder=4)

def figure_for_main_text_v1(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names):
    # Two wide loading rows with narrow association panels, then the scatter.
    fig = plt.figure(figsize=(8.27, 7.27), constrained_layout=True)
    fig.get_layout_engine().set(w_pad=0.01, wspace=0.01)
    gs = fig.add_gridspec(3, 2, width_ratios=[4, 1], height_ratios=[1, 1, 3], wspace=0.03)
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    if len(feature_names) != loadings.shape[0]:
        raise ValueError('Labels must match the species columns used for PCA')
    lineage_labels = {'ab_S-AUS_Cx_annu': 'Cx. annulirostris (S-AUS)',
                      'ab_AUS_Cx_annu': 'Cx. annulirostris (AUS)'}
    species_labels = [lineage_labels[name] if name in lineage_labels else
                      shorted_genus(name.removeprefix('ab_').replace('_', ' '))
                      for name in feature_names]

    def contrast_label(row):
        if row['reference_level'] is not None:
            level = row['term'].split('[T.', 1)[1].removesuffix(']')
            variable = 'Location' if row['variable'] == 'location_type' else 'Season'
            return f"{variable}: {level}\nvs {row['reference_level']}".replace('_', ' ')
        return {'distance_to_coast_km_log10': 'Log10 coast distance',
                'elevation_m': 'Elevation (m)', 'rainfall_mm_3w': 'Rainfall (3 weeks)'}.get(row['variable'], row['variable'])

    loading_axes, association_axes = [], []
    order = np.argsort(loadings[:, 0], kind='stable')[::-1]
    for i in range(2):
        loading_ax = fig.add_subplot(gs[i, 0],
                                     sharex=loading_axes[0] if loading_axes else None,
                                     sharey=loading_axes[0] if loading_axes else None)
        association_ax = fig.add_subplot(gs[i, 1], sharex=association_axes[0] if association_axes else None)
        loading_axes.append(loading_ax)
        association_axes.append(association_ax)
        loading_ax.bar(np.arange(len(order)), loadings[order, i], color='dimgray', width=0.8)
        loading_ax.set_xticks(np.arange(len(order)), labels=[species_labels[j] for j in order],
                              rotation=90, fontsize=6.5, fontstyle='italic')
        loading_ax.set_ylabel(f'PC{i+1} scaled loading', fontsize=9)
        loading_ax.axhline(0, color='0.7', linewidth=0.6)
        loading_ax.margins(x=0.01)
        betas = assoc_coefficients.filter(
            (pl.col('PC') == f'PC{i+1}') & pl.col('variable').is_in(variables_to_plot)
        ).sort('beta')
        positions = np.arange(betas.height)
        beta = betas['beta'].to_numpy()
        association_ax.errorbar(
            beta, positions, xerr=[beta - betas['ci_lower'].to_numpy(),
                                   betas['ci_upper'].to_numpy() - beta],
            fmt='o', markersize=3, color='black', ecolor='lightgrey', elinewidth=2,
        )
        association_ax.set_yticks(positions,
                                  labels=[contrast_label(row) for row in betas.iter_rows(named=True)],
                                  fontsize=6.5)
        association_ax.yaxis.tick_right()
        association_ax.axvline(0, color='0.7', linewidth=0.6, linestyle='--')
        association_ax.set_xlabel('Beta (95% CI)\nContinuous: per 1 SD', fontsize=7)
        association_ax.margins(y=0.08)

    loading_axes[0].tick_params(axis='x', labelbottom=False)
    scatter_ax = fig.add_subplot(gs[2, 0])
    sns.scatterplot(data=pca_scores.to_pandas(), x='PC1', y='PC2', ax=scatter_ax,
                    s=25, color='grey', alpha=0.8, edgecolor='white', linewidth=0.35, zorder=3)
    virus_colours = {
        'MVEV': 'tab:red', 'JEV': 'tab:orange', 'KUN': 'tab:green', 'RRV': 'tab:blue',
        'BFV': 'tab:purple', 'SINDV': 'tab:brown', 'STRV': 'tab:pink',
        'KOK': 'tab:cyan', 'EHV': 'tab:olive', 'LSV': 'black',
    }
    plot_virus_positives_on_PCA_with_colours(df, pca_scores, fig=fig, ax=scatter_ax,
                                           colour_dict=virus_colours)
    scatter_ax.set_xlabel(f'-PC1 ({explained_variance[0]:.1%})', fontsize=9)
    scatter_ax.set_ylabel(f'PC2 ({explained_variance[1]:.1%})', fontsize=9)
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                          markersize=5, label=virus) for virus, color in virus_colours.items()]
    scatter_ax.legend(handles=handles, title='Virus detected', fontsize=6.5,
                      title_fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    for ax in loading_axes + association_axes + [scatter_ax]:
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(axis='both', length=0)
        ax.tick_params(axis='y' if ax in loading_axes else 'x', labelsize=7)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    fig.savefig(out.with_suffix('.pdf'))
    plt.close(fig)

def plot_loading_barplot(loadings, feature_names, pc_index=0, ax=None, orientation='vertical'):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    order = np.argsort(loadings[:, pc_index], kind='stable')[::-1]
    species_labels = fix_species_labels(feature_names)
    labels = [species_labels[j] for j in order]
    if orientation == 'vertical':
        ax.bar(np.arange(len(order)), loadings[order, pc_index], color='dimgray', width=0.8)
        ax.set_xticks(np.arange(len(order)), labels=labels, rotation=90, fontsize=6.5, fontstyle='italic')
        ax.set_ylabel(f'loading', fontsize=9)
        ax.axhline(0, color='0.7', linewidth=0.6)
    else:
        ax.barh(np.arange(len(order)), loadings[order, pc_index], color='dimgray', height=0.8)
        ax.set_yticks(np.arange(len(order)), labels=labels, fontsize=6.5, fontstyle='italic')
        ax.set_xlabel(f'loading', fontsize=9)
        ax.axvline(0, color='0.7', linewidth=0.6)
    # ax.margins(x=0.01)
    # remove spines, keep ticks
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ax

def fix_species_labels(feature_names):
    lineage_labels = {'ab_S-AUS_Cx_annu': 'Cx. annulirostris (S-AUS)',
                      'ab_AUS_Cx_annu': 'Cx. annulirostris (AUS)'}
    species_labels = [lineage_labels[name] if name in lineage_labels else
                      shorted_genus(name.removeprefix('ab_').replace('_', ' '))
                      for name in feature_names]
    return species_labels

def annotate_scatter_with_virus_centroids(df, pca_scores, ax, colour_dict=None, label_locations=None):
    if colour_dict is None:
        return
    if label_locations is None:
        # (horizontal, vertical) offsets in points; +x is right, +y is up
        # on the page, including when the scatter's x-axis is inverted.
        label_locations = {
            'MVEV': (-23, 4), 'JEV': (-6, 18), 'KUN': (8, 8),
            'RRV': (0, 18), 'BFV': (22, 4), 'SINDV': (0, 18),
            'STRV': (10, 18), 'KOK': (6, -18), 'EHV': (-25, -12), 'LSV': (-10, -18),
        }
    for virus, colour in colour_dict.items():
        virus_col = f"{virus}_detected"
        virus_samples = df.filter(pl.col(virus_col) == 1).select('sample').to_numpy().flatten()
        virus_pca_scores = pca_scores.filter(pl.col('sample').is_in(virus_samples)).to_pandas()
        if not virus_pca_scores.empty:
            centroid_x = virus_pca_scores['PC1'].mean()
            centroid_y = virus_pca_scores['PC2'].mean()
            n = len(virus_pca_scores)
            label = f'{virus}\n(n={n})'
            ax.scatter(centroid_x, centroid_y, s=180, color=colour, alpha=0.65,
                       edgecolor='black', linewidth=1.2, zorder=5)
            dx, dy = label_locations.get(virus, (8, 8))
            ax.annotate(label, xy=(centroid_x, centroid_y),
                        xytext=(dx, dy), textcoords='offset points',
                        ha='center', va='center',
                        fontsize=7, fontweight='bold', color=colour, zorder=6,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

def draw_PCA_on_main_figure(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names, ax):
    # plot PCA scatter
    sns.scatterplot(data=pca_scores.to_pandas(), x='PC1', y='PC2', ax=ax, s=10, color='black', alpha=0.8, edgecolor='white', linewidth=0.35, zorder=3)
    virus_colours = {
        'JEV': 'tab:orange', 'KUN': 'tab:green', 'RRV': 'tab:blue',
        'BFV': 'tab:purple', 'SINDV': 'tab:brown', 'STRV': 'tab:pink',
        'KOK': 'tab:cyan', 'EHV': 'tab:olive', 'LSV': 'black','MVEV': 'tab:red'
    }
    plot_virus_positives_on_PCA_with_colours(df, pca_scores, fig=ax.figure, ax=ax, colour_dict=virus_colours)
    ax.set_xlabel(f'-PC1 ({explained_variance[0]:.1%})', fontsize=9)
    ax.set_ylabel(f'PC2 ({explained_variance[1]:.1%})', fontsize=9, labelpad=15)
    # remove spines, ticks and tick labels, add faint grey lines at x=0 and y=0
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis='both', length=0, labelsize=7)
    original_xlim, original_ylim = ax.get_xlim(), ax.get_ylim()
    # Fixed data-coordinate segments keep the grey lines out of the padding.
    ax.plot(original_xlim, [0, 0], color='lightgrey', linewidth=0.5, zorder=1)
    ax.plot([0, 0], original_ylim, color='lightgrey', linewidth=0.5, zorder=1)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    annotate_scatter_with_virus_centroids(df, pca_scores, ax, colour_dict=virus_colours)
    padding = 0.0  # Extra space on each side, relative to the original axis span.
    dx = (original_xlim[1] - original_xlim[0]) * padding
    dy = (original_ylim[1] - original_ylim[0]) * padding
    ax.set_xlim(original_xlim[0] - dx, original_xlim[1] + dx)
    ax.set_ylim(original_ylim[0] - dy, original_ylim[1] + dy)

def draw_beta_forrest_plot(assoc_coefficients, pc_index=0, variables_to_plot=None, ax=None):
    from matplotlib.transforms import blended_transform_factory
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))
    pc = f'PC{pc_index + 1}'
    data = assoc_coefficients.filter(pl.col('PC') == pc)
    if variables_to_plot is not None:
        data = data.filter(pl.col('variable').is_in(variables_to_plot))
    groups = {
        'Environmental': ['distance_to_coast_km_log10', 'E_or_W_of_GDR', 'location_type', 'elevation_m', 'rainfall_mm_3w'],
        'Temporal': ['season', 'week'],
    }
    labels = {'distance_to_coast_km_log10': 'Coast distance (log₁₀)',
              'elevation_m': 'Elevation', 'rainfall_mm_3w': 'Rainfall', 'week': 'Week',
              'E_or_W_of_GDR': 'W vs E of GDR',}
    rows, headings = [], []
    y = 0.0
    for group, variables in groups.items():
        group_rows = []
        for variable in variables:
            for row in data.filter(pl.col('variable') == variable).sort('term').iter_rows(named=True):
                if row['reference_level'] is not None:
                    level = row['term'].split('[T.', 1)[1].removesuffix(']').replace('_', ' ')
                    reference = row['reference_level'].replace('_', ' ')
                    if variable == 'season':
                        level = f'20{level[:2]}-{level[2:]}'
                        reference = f'20{reference[:2]}-{reference[2:]}'
                    label = f'{level} vs {reference}'
                    if variable == 'location_type':
                        label = label.replace('Metropolitan', 'Met.').replace('Met Fringe', 'fringe')
                    elif variable == 'E_or_W_of_GDR':
                        label = f'GDR: {level} vs {reference}'
                else:
                    label = labels.get(variable, variable)
                group_rows.append((row, label))
        if not group_rows:
            continue
        headings.append((group, y))
        y += 0.85
        for row, label in group_rows:
            rows.append((row, label, y))
            y += 1
        y += 0.6  # Small gap before the next group.
    ax.set_axis_off()
    if not rows:
        ax.text(0.5, 0.5, f'No coefficients for {pc}', ha='center', transform=ax.transAxes)
        return ax
    # Keep labels, forest, R² and p inside the supplied panel, rather than
    # letting long tick labels push neighbouring panels apart.
    forest = ax.inset_axes([0.52, 0.09, 0.23, 0.80])
    forest.set_ylim(y - 0.2, -0.9)
    text_transform = blended_transform_factory(ax.transAxes, forest.transData)
    ax.text(0.02, -1.15, '-PC1' if pc == 'PC1' else pc, fontsize=9, fontweight='bold', transform=text_transform)
    ax.text(0.88, -1.15, r'$R^2$', fontsize=8, ha='right', transform=text_transform)
    ax.text(1.04, -1.15, 'Adj. p', fontsize=8, ha='right', transform=text_transform)
    for group, position in headings:
        ax.text(0.02, position - 0.2, group, fontsize=8, fontweight='bold',
                va='center', transform=text_transform)
    for row, label, position in rows:
        beta, low, high = row['beta'], row['ci_lower'], row['ci_upper']
        significant = (row['p_value_fdr'] is not None and row['p_value_fdr'] < 0.05
                       and (low > 0 or high < 0))
        forest.errorbar(beta, position, xerr=[[beta-low], [high-beta]],
                        fmt='o', markersize=6 if significant else 4, color='black', ecolor='0.65', elinewidth=1.4,
                        markerfacecolor='black' if significant else 'white',
                        markeredgecolor='black')
        ax.text(0.02, position, label, fontsize=7, va='center', transform=text_transform)
    # Model-level values appear once, centred across each variable's contrasts.
    for variable in dict.fromkeys(row['variable'] for row, _, _ in rows):
        variable_rows = [(row, position) for row, _, position in rows if row['variable'] == variable]
        row = variable_rows[0][0]
        first, last = variable_rows[0][1], variable_rows[-1][1]
        centre = (first + last) / 2
        r_squared, p_value = row['r_squared'], row['p_value_fdr']
        ax.text(0.88, centre, f'{r_squared:.2f}'.removeprefix('0') if r_squared is not None else 'NA',
                fontsize=7, ha='right', va='center', transform=text_transform)
        p_label = 'NA' if p_value is None else ('<.001' if p_value < 0.001 else f'{p_value:.3f}'.removeprefix('0'))
        ax.text(1.04, centre, p_label, fontsize=7, ha='right', va='center', transform=text_transform)
        if row['reference_level'] is not None and len(variable_rows) > 1:
            # Square bracket to the left of the shared model statistics.
            for x in [0.80]:
                ax.plot([x - 0.012, x, x, x - 0.012], [first, first, last, last],
                        color='0.35', linewidth=0.7, transform=text_transform, clip_on=False)
    forest.axvline(0, color='0.7', linewidth=0.6, linestyle='--')
    forest.set_yticks([])
    forest.tick_params(axis='x', labelsize=7, length=2)
    from matplotlib.ticker import MaxNLocator
    forest.xaxis.set_major_locator(MaxNLocator(2))
    forest.set_xlabel('Beta (95% CI)', fontsize=7, labelpad=2)
    for spine in forest.spines.values():
        spine.set_visible(False)
    return ax


def figure_for_main_text_v2(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names):
    # gridspec 3 cols and 4 rows
    # Fixed panel positions prevent automatic layout from widening the gaps.
    fig = plt.figure(figsize=(11.69, 7.5))
    gs = fig.add_gridspec(4, 3, width_ratios=[0.3, 1.7, 1.4],
                          height_ratios=[1, 1, 1, 1], left=0.27, right=0.98,
                          bottom=0.30, top=0.94, wspace=0.13, hspace=0.35)
    # first col, first three rows
    PC2_loadings_ax = fig.add_subplot(gs[0:3, 0])
    PC1_loadings_ax = fig.add_subplot(gs[3, 1])
    bottom_left_corner_ax = fig.add_subplot(gs[3, 0])
    PCA_ax = fig.add_subplot(gs[0:3, 1])
    PC1_assoc_ax = fig.add_subplot(gs[0:2, 2])
    PC2_assoc_ax = fig.add_subplot(gs[2:4, 2])
    # empty botom left corner
    bottom_left_corner_ax.axis('off')
    bottom_left_corner_ax.text(0.5, 0.5, 'Loading', ha='center', va='center',
                               fontsize=9, transform=bottom_left_corner_ax.transAxes)
    # Bottom strip: species along x, with bars extending vertically.
    plot_loading_barplot(pca.components_.T * np.sqrt(pca.explained_variance_), feature_names, pc_index=0, ax=PC1_loadings_ax, orientation='vertical')
    # Fix the labels
    plt.setp(PC1_loadings_ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    # Left strip: species along y, with bars extending horizontally.
    plot_loading_barplot(pca.components_.T * np.sqrt(pca.explained_variance_), feature_names, pc_index=1, ax=PC2_loadings_ax, orientation='horizontal')
    PC2_loadings_ax.invert_yaxis()
    plt.setp(PC2_loadings_ax.get_yticklabels(), rotation=0, fontsize=7)
    from matplotlib.ticker import MultipleLocator
    loading_values = pca.components_.T * np.sqrt(pca.explained_variance_)
    for i, loading_ax in enumerate([PC1_loadings_ax, PC2_loadings_ax]):
        low, high = loading_values[:, i].min(), loading_values[:, i].max()
        padding = (high - low) * 0.08
        if i == 0:
            loading_ax.set_ylim(low - padding, high + padding)
            loading_ax.yaxis.set_major_locator(MultipleLocator(0.1))
            loading_ax.tick_params(axis='y', labelsize=6)
            loading_ax.set_ylabel('')
        else:
            loading_ax.set_xlim(low - padding, high + padding)
            loading_ax.xaxis.set_major_locator(MultipleLocator(0.1))
            loading_ax.tick_params(axis='x', labelsize=6)
            loading_ax.set_xlabel('')
    # plot PCA scatter
    draw_PCA_on_main_figure(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names, ax=PCA_ax)


    draw_beta_forrest_plot(assoc_coefficients, pc_index=0,
                           variables_to_plot=variables_to_plot, ax=PC1_assoc_ax)
    draw_beta_forrest_plot(assoc_coefficients, pc_index=1,
                           variables_to_plot=variables_to_plot, ax=PC2_assoc_ax)

    # Scale every label consistently, including text in the inset forest axes.
    from matplotlib.text import Text
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize() * 1.25)

    # savefig
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)

def pairplot_variables(df, variables, out_path):
    df = df.select(variables).to_pandas()
    g = sns.pairplot(df, diag_kind='hist', plot_kws={'alpha': 0.5, 's': 10})
    g.fig.suptitle('Pairplot of selected variables', y=1.02)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    g.savefig(out, dpi=300, bbox_inches='tight')
    g.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(g.fig)

def plot_box_plot_of_location_type_vs_D1_with_regline(df, grouping_variable, out_path):
    # perform regression of location type on D1: simpsons_diversity_index ~ C(location_type) + (1| location-site)
    # filter out rows with missing values in the relevant columns
    df = df.drop_nulls(['simpson_diversity_index', 'location_type', 'distance_to_coast_km_log10', grouping_variable]).to_pandas()
    formula = 'simpson_diversity_index ~ C(location_type) + distance_to_coast_km_log10'
    model = smf.mixedlm(formula, data=df, groups=df[grouping_variable]).fit()
    # plot the results
    fig, axs = plt.subplots(figsize=(6, 2.5), ncols=2, nrows=1, gridspec_kw={'width_ratios': [1, 2], 'wspace': 0.3}, sharey=True)
    sns.boxplot(x='location_type', y='simpson_diversity_index', data=df, ax=axs[0], palette='Set2',
                width=0.65, linewidth=0.8, showfliers=False)
    # add stripplot of the points on top of the boxplot
    sns.stripplot(x='location_type', y='simpson_diversity_index', data=df, ax=axs[0], color='black', size=2, jitter=True, alpha=0.4)
    axs[0].set_xlabel('Location type', fontsize=8)
    axs[0].set_ylabel('Simpson diversity index (D1)', fontsize=8)
    # rotate x-axis labels for better readability
    plt.setp(axs[0].get_xticklabels(), rotation=35, ha='right', rotation_mode='anchor')
    # annotate p, R2, and n on the plot
    # Joint test of all location contrasts, adjusted for distance to coast.
    terms = [term for term in model.fe_params.index if term.startswith('C(location_type)')]
    restriction = np.zeros((len(terms), len(model.params)))
    for i, term in enumerate(terms):
        restriction[i, model.params.index.get_loc(term)] = 1
    p_value = float(model.wald_test(restriction, scalar=True).pvalue)
    fixed_variance = float(np.var(model.model.exog @ model.fe_params.to_numpy(), ddof=1))
    r_squared = fixed_variance / (fixed_variance + float(model.cov_re.iloc[0, 0]) + float(model.scale))
    n = int(model.nobs)
    axs[0].text(0, 1.06,
                f'Location type: p = {p_value:.3g}\nMarginal R² = {r_squared:.3f}; n = {n}\nCoast-adjusted; site random intercept',
                transform=axs[0].transAxes, fontsize=6.5, linespacing=1.4, ha='left', va='bottom')
    # also scatter D1 of y and distance to coast on x, hue by location type
    axs[1].set_xlabel('Distance to coast (log10 km)')
    axs[1].set_ylabel('')
    sns.scatterplot(x='distance_to_coast_km_log10', y='simpson_diversity_index', data=df, hue='location_type', ax=axs[1], palette='Set2', s=14, alpha=0.75, linewidth=0.3)
    axs[1].invert_xaxis()
    # Coordinates are already log10-transformed; label whole decades in km.
    coast_limits = axs[1].get_xlim()
    coast_ticks = np.arange(np.ceil(min(coast_limits)), np.floor(max(coast_limits)) + 1)
    axs[1].set_xticks(coast_ticks, labels=[f'{10.0 ** tick:g}' for tick in coast_ticks])
    axs[1].set_xlim(coast_limits)
    axs[1].set_xlabel('Distance to coast (km)', fontsize=8)
    axs[1].legend(title='Location type', loc='lower center', bbox_to_anchor=(0.5, 1.04),
                  ncol=2, frameon=False, fontsize=6.5, title_fontsize=7,
                  markerscale=0.8, handletextpad=0.4, columnspacing=1)
    fig.subplots_adjust(top=0.72, bottom=0.25)
    # despine both axes
    for ax in axs:
        ax.tick_params(axis='both', labelsize=7, length=3, width=0.6)
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='0.93', linewidth=0.5)
        sns.despine(ax=ax, trim=True, top=True, right=True)
    # save the figure
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)

def main():
    # load all data
    merged_df = make_df_for_TG()

    # add elevation and rainfall data for later
    merged_df = add_elevation_and_rainfall_to_data(merged_df)

    # add log10 transformed distance to coast column and other steps
    merged_df = log_transform_desired_columns(merged_df, ['distance_to_coast_km'])
    merged_df = clean_df(merged_df)

    # add E or W of the Great Dividing Range column
    merged_df = publication_figres.add_whether_east_or_west_of_GDR(merged_df)

    # run PCA on mosquito abundance
    mosquito_abundance_cols = [ col for col in merged_df.columns if col.startswith("ab_") and col != "ab_Culex_annulirostris" ]  # remove this column if it exists
    # mosquito_abundance_cols = [col for col in merged_df.columns if col.startswith("ab_") and col not in ["ab_S-AUS_Cx_annu", "ab_AUS_Cx_annu"]]  # remove these two columns if they exist
    assert len(mosquito_abundance_cols) == 24
    mosquito_abundance_df = merged_df.select(["sample"] + mosquito_abundance_cols)
    pca_scores, explained_variance, pca = run_PCA_on_mosquito_abundance(mosquito_abundance_df, n_components=5)
    # plot PC1-2, 2-3, 3-4, 4-5
    plot_multipanel_PCA(pca_scores, explained_variance_ratio=explained_variance,
                        save_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/multipanel_PCA.png")

    # plot PCA loadings
    plot_PCA_loadings(pca, feature_names=mosquito_abundance_cols, save_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/PCA_loadings.png")

    # simple association of PCs with variables
    # merged_df = merged_df.with_columns(
    #     (pl.col('distance_to_coast_km') + 1e-10).log10().alias('distance_to_coast_km_log10'),
    #     pl.col('season_week').str.split('_').list.get(0).cast(pl.Int32).alias('season'),
    #     pl.col('season_week').str.split('_').list.get(1).cast(pl.Int16).alias('week'),
    # )

    continuous_variables = ['distance_to_coast_km_log10', 'species_count_in_trap', 'simpson_diversity_index', 'elevation_m', 'rainfall_mm_3w', 'week', 'council_olg_code']
    categorical_variables = ['season', 'location_type', 'any_virus_detected', 'E_or_W_of_GDR']
    grouping_variable = 'location-site'
    association_results, association_coefficients = associate_PCs_with_variables(merged_df, pca_scores, continuous_variables=continuous_variables, categorical_variables=categorical_variables, grouping_variable=grouping_variable, return_coefficients=True)
    association_results.write_csv("results/abundance-geospatial/PCA_loadings_with_virus/PCA/association_results.csv")
    association_coefficients.write_csv("results/abundance-geospatial/PCA_loadings_with_virus/PCA/association_coefficients.csv")

    # plot figure for main text
    variables_to_plot = ['season', 'location_type', 'distance_to_coast_km_log10', 'elevation_m', 'rainfall_mm_3w', 'week', 'E_or_W_of_GDR']
    figure_for_main_text_v1(merged_df, pca_scores, explained_variance, pca, association_coefficients, variables_to_plot,
                         out_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/figure_for_main_text_v1.png",
                         feature_names=mosquito_abundance_cols)

    figure_for_main_text_v2(merged_df, pca_scores, explained_variance, pca, association_coefficients, variables_to_plot,
                         out_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/figure_for_main_text_v2.png",
                         feature_names=mosquito_abundance_cols)

    # make a big pairplot of all variables of interest
    pairplot_variables(merged_df, variables=continuous_variables + categorical_variables, out_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/pairplot_variables.png")

    # supp lementary figure: boxplot of location type vs D1 with regression line
    plot_box_plot_of_location_type_vs_D1_with_regline(merged_df, grouping_variable=grouping_variable, out_path="results/abundance-geospatial/PCA_loadings_with_virus/PCA/plots/location_type_vs_D1_boxplot.png")


if __name__ == "__main__":
    main()
