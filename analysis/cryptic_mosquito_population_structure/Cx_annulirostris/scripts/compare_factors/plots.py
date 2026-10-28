import seaborn as sns
from matplotlib import pyplot as plt
import polars as pl

def plot_scatter(df: pl.DataFrame, x: str, y: str, ax: plt.Axes, plot_x_eq_y: bool, hue_col: str = None) -> None:
    x_data = df[x].to_list()
    y_data = df[y].to_list()
    if hue_col is not None:
        hue_data = df[hue_col].to_list()
        sns.scatterplot(x=x_data, y=y_data, hue=hue_data, ax=ax)
    else:
        sns.scatterplot(x=x_data, y=y_data, ax=ax)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if plot_x_eq_y:
        ax.plot(df[x], df[x], color='dimgray', linestyle='--', alpha=0.5)

def scatter_lat_lon_coloured_by_snp_ratios(df: pl.DataFrame, lat_col: str, lon_col: str, snp_ratio_col: str, output_file: str) -> None:
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=df[lon_col].to_numpy(), y=df[lat_col].to_numpy(), hue=df[snp_ratio_col].to_numpy(), palette="viridis")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.legend(title=snp_ratio_col)
    sns.despine(top=True, right=True, trim=True)
    plt.tight_layout()
    plt.savefig(output_file)

def scatter_snp_ratios_and_spatial_factors(df: pl.DataFrame, spatial_factor_cols: list, snp_ratio_col: str, output_file: str) -> None:
    n_factors = len(spatial_factor_cols)
    n_cols = 3
    n_rows = (n_factors + n_cols - 1) // n_cols
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4))
    axs = axs.flatten()
    for i, spatial_col in enumerate(spatial_factor_cols):
        plot_scatter(df, x=spatial_col, y=snp_ratio_col, ax=axs[i], plot_x_eq_y=False)
    # remove unused axes
    for i in range(n_factors, len(axs)):
        fig.delaxes(axs[i])
    sns.despine(top=True, right=True, trim=True, fig=fig)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    # save
    fig.savefig(output_file)

def swarm_plot_snp_ratios_by_factor(df: pl.DataFrame, factor_col: str, snp_ratio_col: str, output_file: str) -> None:
    plt.figure(figsize=(6, 6))
    sns.swarmplot(x=df[factor_col].to_numpy(), y=df[snp_ratio_col].to_numpy())
    plt.xlabel(factor_col)
    plt.ylabel(snp_ratio_col)
    sns.despine(top=True, right=True, trim=True)
    plt.tight_layout()
    plt.savefig(output_file)

def virus_swarmplot_and_ct_scatter(df: pl.DataFrame, virus_col: str, snp_ratio_col: str, ct_col: str) -> plt.Figure:
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    sns.swarmplot(x=df[virus_col].to_numpy(), y=df[snp_ratio_col].to_numpy(), ax=axs[0])
    axs[0].set_xlabel(virus_col)
    axs[0].set_ylabel(snp_ratio_col)
    sns.scatterplot(x=df[ct_col].to_numpy(), y=df[snp_ratio_col].to_numpy(), ax=axs[1])
    axs[1].set_xlabel(ct_col)
    axs[1].set_ylabel(snp_ratio_col)
    sns.despine(top=True, right=True, trim=True, fig=fig)
    fig.tight_layout()
    return fig

def scatter_lat_lon_coloured_by_snp_ratio_split_by_season(df: pl.DataFrame, lat_col: str, lon_col: str, snp_ratio_col: str, season_col: str) -> plt.Figure:
    ncols = df[season_col].n_unique()
    fig, axs = plt.subplots(1, ncols, figsize=(6 * ncols, 6))
    for i, season in enumerate(df[season_col].unique()):
        season_df = df.filter(pl.col(season_col) == season)
        sns.scatterplot(x=season_df[lon_col].to_numpy(), y=season_df[lat_col].to_numpy(), hue=season_df[snp_ratio_col].to_numpy(), palette="viridis", ax=axs[i])
        axs[i].set_title(f"Season: {season}")
        axs[i].set_xlabel("Longitude")
        axs[i].set_ylabel("Latitude")
        axs[i].legend(title=snp_ratio_col)
    sns.despine(top=True, right=True, trim=True, fig=fig)
    fig.tight_layout()
    return fig