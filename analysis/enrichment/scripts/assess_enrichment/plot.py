# plot coverage across COX1

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import os
from analysis.common.io import parse_json
from analysis.common.mosquito_ids import shorted_genus
from analysis.virus_detection.scripts.plot import plot_line_of_best_fit_to_ax

def plot_lineplot_to_ax(data: pl.DataFrame, x: str, y: str, ax: plt.Axes, label=None, color=None, legend=True, alpha=0.8) -> None:
    sns.lineplot(data=data.to_pandas(), x=x, y=y, ax=ax, label=label, color=color, alpha=alpha)
    ax.set_xlabel(x)
    ax.set_ylabel(y)

def plot_scatter_with_line_of_best_fit_to_ax(data: pl.DataFrame, x: str, y: str, ax: plt.Axes, logx=False, logy=False) -> None:
    sns.scatterplot(data=data.to_pandas(), x=x, y=y, ax=ax, color='k')
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    plot_line_of_best_fit_to_ax(data, x=x, y=y, ax=ax, logx=logx, logy=logy)

def _format_log_axis(ax: plt.Axes, axis: str) -> None:
    log_axis = getattr(ax, f"{axis}axis")
    log_axis.set_major_locator(mticker.LogLocator(base=10))
    log_axis.set_major_formatter(mticker.FuncFormatter(lambda value, _: f"{int(value):,}" if value >= 1 else f"{value:g}"))
    log_axis.set_minor_locator(mticker.NullLocator())

def _style_ax(ax: plt.Axes, grid_axis: str = "both") -> None:
    ax.grid(True, which="major", axis=grid_axis, linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _add_cox1_boundary_lines(ax: plt.Axes, coverage_df: pl.DataFrame, position_col: str) -> None:
    boundaries = (
        coverage_df
        .filter(pl.col("is_cox1"))
        .group_by("ref_genome")
        .agg(
            pl.col(position_col).min().alias("start"),
            pl.col(position_col).max().alias("end"),
        )
        .select("start", "end")
        .unique()
    )
    if boundaries.height != 1:
        raise ValueError(f"Expected one shared COX1 interval, found: {boundaries.to_dicts()}")
    start, end = boundaries.row(0)
    ax.axvline(start, color="tab:red", linestyle="--", linewidth=2)
    ax.axvline(end, color="tab:red", linestyle="--", linewidth=2)


def _average_coverage_at_well_represented_positions(coverage_df: pl.DataFrame, position_col: str, min_sample_fraction: float) -> pl.DataFrame:
    total_samples = coverage_df.select(pl.col("sample").n_unique()).item()
    min_samples = max(1, int(np.ceil(total_samples * min_sample_fraction)))
    average_coverage = coverage_df.group_by(position_col).agg(
        pl.col("coverage").mean().alias("average_coverage"),
        pl.col("sample").n_unique().alias("samples_with_base"),
    )
    all_positions = pl.DataFrame({
        position_col: range(coverage_df[position_col].min(), coverage_df[position_col].max() + 1)
    })
    return (
        all_positions.join(average_coverage, on=position_col, how="left")
        .with_columns(
            pl.when(pl.col("samples_with_base") >= min_samples)
            .then(pl.col("average_coverage"))
            .otherwise(None)
            .alias("average_coverage")
        )
        .select(position_col, "average_coverage")
        .sort(position_col)
    )

def plot_enrichment_over_COX1(coverage_df: pl.DataFrame, fold_enrichment_df: pl.DataFrame, output_dir: str, min_average_sample_fraction: float = 0.8) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 4), width_ratios=[3, 1.4])
    axs = axs.flatten()
    enrichment_data = []
    position_col = "aligned_position" if "aligned_position" in coverage_df.columns else "position"
    for sample in coverage_df['sample'].unique():
        plot_lineplot_to_ax(coverage_df.filter(pl.col("sample") == sample), x=position_col, y="coverage", label=sample, ax=axs[0], color='tab:grey', alpha=0.2, legend=False)
    # plot average line in red
    avg_coverage_by_position = _average_coverage_at_well_represented_positions(coverage_df, position_col, min_average_sample_fraction)
    avg_coverage_by_position_pd = avg_coverage_by_position.to_pandas()
    axs[0].plot(avg_coverage_by_position_pd[position_col], avg_coverage_by_position_pd["average_coverage"], color="black", alpha=1, linewidth=2, label="Average")
    _add_cox1_boundary_lines(axs[0], coverage_df, position_col)
    # log scale y axis
    axs[0].set_yscale("log")
    _format_log_axis(axs[0], "y")
    _style_ax(axs[0])
    # label axes
    axs[0].set_xlabel("aligned position on mt genome" if position_col == "aligned_position" else "position on mt genome")
    axs[0].set_ylabel("read depth")
    # have no legend 
    axs[0].legend().set_visible(False)
    # make legend 3 cols
    # handles, labels = axs[0].get_legend_handles_labels()
    # new_labels = [shorted_genus(label.replace("_", " ")) for label in labels]
    # axs[0].legend(handles, new_labels, title="Sample", loc="upper right", bbox_to_anchor=(1.3, 1), ncol=3)
    # now draw scatter and line of best fit of median coverage in range vs median coverage out of range
    # make a fold enrichment column
    plot_scatter_with_line_of_best_fit_to_ax(fold_enrichment_df, y="mean_in_cox1", x="mean_out_cox1", ax=axs[1])
    # log scale axes
    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    _format_log_axis(axs[1], "x")
    _format_log_axis(axs[1], "y")
    _style_ax(axs[1])
    # label axes
    axs[1].set_xlabel("mean background depth")
    axs[1].set_ylabel("mean enriched depth")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/enrichment_over_COX1.png")

def plot_enrichment_over_COX1_with_BFV_culture(coverage_df: pl.DataFrame, BFV_results_path: str, output_dir: str, min_average_sample_fraction: float = 0.8) -> None:
    # parse BFV data
    BFV_results = pl.read_csv(BFV_results_path)
    # plot coverage across COX1
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 4), width_ratios=[1.4, 3])
    axs = axs.flatten()
    position_col = "aligned_position" if "aligned_position" in coverage_df.columns else "position"
    for sample in coverage_df['sample'].unique():
        plot_lineplot_to_ax(coverage_df.filter(pl.col("sample") == sample), x=position_col, y="coverage", label=sample, ax=axs[1], color='tab:grey', alpha=0.2, legend=False)
    # plot average line in red
    avg_coverage_by_position = _average_coverage_at_well_represented_positions(coverage_df, position_col, min_average_sample_fraction)
    avg_coverage_by_position_pd = avg_coverage_by_position.to_pandas()
    axs[1].plot(avg_coverage_by_position_pd[position_col], avg_coverage_by_position_pd["average_coverage"], color="black", alpha=1, linewidth=2, label="Average")
    _add_cox1_boundary_lines(axs[1], coverage_df, position_col)
    # log scale y axis
    axs[1].set_yscale("log")
    _format_log_axis(axs[1], "y")
    _style_ax(axs[1])
    # label axes
    axs[1].set_xlabel("aligned position on mt genome" if position_col == "aligned_position" else "position on mt genome")
    axs[1].set_ylabel("read depth")
    # have no legend 
    axs[1].legend().set_visible(False)
    # plot BFV scatter
    plot_scatter_with_line_of_best_fit_to_ax(BFV_results, y="n_reads_dedup", x="ct", ax=axs[0], logx=False, logy=True)
    # _format_log_axis(axs[1], "x")
    _format_log_axis(axs[0], "y")
    _style_ax(axs[0])
    # label axes
    axs[0].set_xlabel("cycle threshold in original\nsample")
    axs[0].set_ylabel("unique read count")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/enrichment_over_COX1_plus_BFV.png")
