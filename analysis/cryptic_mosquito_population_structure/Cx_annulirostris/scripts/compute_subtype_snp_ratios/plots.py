import seaborn as sns
from matplotlib import pyplot as plt
import polars as pl

def scatter_snp_ratios_and_kma(df: pl.DataFrame, x_col: str, y_col: str, plot_file: str) -> None:
    plt.figure(figsize=(6, 6))
    sns.scatterplot(x=df[x_col].to_numpy(), y=df[y_col].to_numpy())
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.plot(df[x_col], df[x_col], color='dimgray', linestyle='--', alpha=0.5)
    sns.despine(top=True, right=True, trim=True)
    plt.tight_layout()
    plt.savefig(plot_file)