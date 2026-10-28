import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from ..virus.plot import virus_coverage_plot_paneled
from matplotlib.ticker import FuncFormatter, NullLocator
from matplotlib.text import Text
from importlib import import_module
import pickle
from .virus_scatters_and_heatmap import plot_nc2_curve_panels

abundance_validation = import_module("analysis.abundance-validation.plot")

A4_WIDTH_INCHES = 8.27

def style_coverage_yaxis(ax: plt.Axes) -> None:
    # log10(coverage + 1) includes zero; ticks remain in read-depth units.
    ax.set_yscale("function", functions=(lambda depth: np.log10(1 + depth),
                                        lambda log_depth: np.power(10.0, log_depth) - 1))
    ax.set_ylim(bottom=-0.1)
    ax.set_yticks([0, 100, 10_000])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.yaxis.set_minor_locator(NullLocator())


def plot_mitochondrial_coverage(mito_df: pl.DataFrame, fig: plt.Figure, subplot_spec: plt.SubplotSpec, zoom_end: int = None, dashed_boundaries: bool = False) -> plt.Axes:
    mitochondrial_grid = subplot_spec.subgridspec(
        1,
        2,
        width_ratios=[0.025, 0.975],
        wspace=0.35,
        hspace=0.05
    )
    mitochondrial_bracket_ax = fig.add_subplot(mitochondrial_grid[0, 0])
    mitochondrial_ax = fig.add_subplot(mitochondrial_grid[0, 1])

    position_col = "aligned_position" if "aligned_position" in mito_df.columns else "position"
    mito_df = mito_df.with_columns(
        pl.col(position_col).cast(pl.Int64),
        pl.col("coverage").cast(pl.Float64),
    )

    total_samples = mito_df.select(pl.col("sample").n_unique()).item()
    minimum_samples = max(1, int(np.ceil(total_samples * 0.8)))
    # Median in log10(depth + 1) space; invert for the axis, which reapplies
    # log10(depth + 1). The displayed height is exactly the log10 median.
    median_coverage_df = mito_df.group_by(position_col).agg(
        (pl.lit(10.0).pow((pl.col("coverage") + 1).log(10).median()) - 1).alias("median_coverage"),
        pl.col("sample").n_unique().alias("n_samples"),
    )
    all_positions_df = pl.DataFrame(
        {
            position_col: range(
                mito_df[position_col].min(),
                mito_df[position_col].max() + 1,
            )
        }
    )
    median_coverage_df = (
        all_positions_df
        .join(median_coverage_df, on=position_col, how="left")
        .with_columns(
            pl.when(pl.col("n_samples") >= minimum_samples)
            .then(pl.col("median_coverage"))
            .otherwise(None)
            .alias("median_coverage")
        )
        .select(position_col, "median_coverage")
        .sort(position_col)
    )

    for sample in mito_df.select("sample").unique().to_series().to_list():
        sample_coverage = (
            mito_df
            .filter(pl.col("sample") == sample)
            .select(position_col, "coverage")
            .sort(position_col)
            .to_pandas()
        )
        mitochondrial_ax.plot(
            sample_coverage[position_col],
            sample_coverage["coverage"],
            color="#BDBDBD",
            alpha=0.5,
            linewidth=0.4,
        )

    median_coverage = median_coverage_df.to_pandas()
    mitochondrial_ax.plot(
        median_coverage[position_col],
        median_coverage["median_coverage"],
        color="black",
        alpha=1,
        linewidth=1.2,
    )

    cox1_boundaries = (
        mito_df
        .filter(pl.col("is_cox1"))
        .select(
            pl.col(position_col).min().alias("start"),
            pl.col(position_col).max().alias("end"),
        )
        .row(0, named=True)
    )
    if dashed_boundaries:
        for boundary in [cox1_boundaries["start"], cox1_boundaries["end"]]:
            mitochondrial_ax.axvline(boundary, color="#A62635", linestyle="--", linewidth=1)
    else:
        mitochondrial_ax.axvspan(
            cox1_boundaries["start"], cox1_boundaries["end"],
            color="#A62635", alpha=0.3, zorder=20,
        )
    # annotate the cox1 region
    mitochondrial_ax.text(
        (cox1_boundaries["start"] + cox1_boundaries["end"]) / 2,
        1.07,
        "COX1",
        color="#A62635",
        fontsize=9,
        fontstyle="italic",
        ha="center",
        va="top",
        transform=mitochondrial_ax.get_xaxis_transform(),
        zorder=5,
    )

    max_position = int(mito_df.select(pl.col(position_col).max()).item())
    max_position = zoom_end if zoom_end is not None else max_position
    tick_interval = 1_500 if zoom_end is not None else 3_000
    mitochondrial_ax.set_xlim(0, max_position)
    mitochondrial_ax.set_xticks(np.arange(0, max_position + 1, tick_interval))
    mitochondrial_ax.set_xlabel("")
    mitochondrial_ax.set_ylabel("")
    mitochondrial_ax.tick_params(axis="both", which="both", width=0.6, length=3, color="#555555")
    mitochondrial_ax.tick_params(axis="x", labelsize=9)
    mitochondrial_ax.tick_params(axis="y", which="major", labelsize=8, pad=1)
    # make the x axis span to the final tick mark
    largest_tick = mitochondrial_ax.get_xticks()[-1]
    style_coverage_yaxis(mitochondrial_ax)
    mitochondrial_ax.set_yticks([0, 10, 100, 1_000, 10_000, 100_000])
    mitochondrial_ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    mitochondrial_ax.spines["top"].set_visible(False)
    mitochondrial_ax.spines["right"].set_visible(False)
    mitochondrial_ax.spines["left"].set_position(("outward", 5))
    mitochondrial_ax.spines["bottom"].set_bounds(0, largest_tick)
    mitochondrial_ax.spines["bottom"].set_position(("outward", 5))

    bracket_position = mitochondrial_bracket_ax.get_position()
    mitochondrial_bracket_ax.set_position([
        bracket_position.x0 - 0.040, bracket_position.y0,
        bracket_position.width, bracket_position.height,
    ])
    label_x = -0.6
    mitochondrial_bracket_ax.text(
        label_x,
        0.52,
        "Mosquito\nmitochondrial\ngenome",
        ha="right",
        va="bottom",
        multialignment="right",
        fontsize=10,
        transform=mitochondrial_bracket_ax.transAxes,
    )

    mitochondrial_bracket_ax.text(
        label_x,
        0.48,
        f"(n={total_samples})",
        ha="right",
        va="top",
        fontsize=9,
        fontstyle="italic",
        transform=mitochondrial_bracket_ax.transAxes,
    )

    mitochondrial_bracket_ax.set_xticks([])
    mitochondrial_bracket_ax.set_yticks([])
    mitochondrial_bracket_ax.patch.set_alpha(0)
    mitochondrial_bracket_ax.spines["right"].set_visible(False)
    for side in ["left", "top", "bottom"]:
        mitochondrial_bracket_ax.spines[side].set_linewidth(0.8)
        mitochondrial_bracket_ax.spines[side].set_color("#555555")
    return mitochondrial_ax


def plot_enrichment_and_estimation(mitochondria_coverage_file: str, virus_coverage_file: str, heatmap_data_file: str, output_file: str) -> None:
    fig = plt.figure(
        figsize=(A4_WIDTH_INCHES, 6.5),
        constrained_layout=False,
    )
    outer_grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[0.54, 0.46],
        height_ratios=[0.32, 0.68],
        left=0.1,
        right=0.9,
        top=0.97,
        bottom=0.17,
        wspace=0.45,
        hspace=0.12,
        )
    mitochondrial_spec = outer_grid[0, 0]
    virus_spec = outer_grid[1, 0]
    abundance_spec = outer_grid[:, 1]

    # parse files
    mito_df = pl.read_csv(mitochondria_coverage_file, separator=",")
    virus_df = pl.read_csv(virus_coverage_file, separator=",")
    with open(heatmap_data_file, "rb") as handle:
        heatmap_data = pickle.load(handle)

    # plot the mitochondrial coverage
    plot_mitochondrial_coverage(mito_df, fig, mitochondrial_spec)

    # plot the virus coverage
    virus_coverage_plot_paneled(virus_df, output_file, fig=fig, subplot_spec=virus_spec)

    # plot the abundance validation heatmap
    abundance_validation.plot_control_abundance_heatmap(
        heatmap_data,
        simulated_mixture_panel="heatmap",
        cmap="Blues",
        fig=fig,
        subplot_spec=abundance_spec,
    )

    # add a shared y axis label on the left side of the figure
    fig.text(
        0.02,
        0.55,
        "Read depth",
        rotation="vertical",
        fontsize=10,
        ha="center",
        va="center",
    )

    # save the figure
    fig.savefig(output_file, dpi=300)

def plot_virus_coverage_paneled(virus_df: pl.DataFrame, fig: plt.Figure, subplot_spec: plt.SubplotSpec) -> list[plt.Axes]:
    virus_acronyms = {
        "stratfordvirus": "STRV",
        "rossrivervirus": "RRV",
        "murrayvalleyencephalitisvirus": "MVEV",
        "kunjinvirus": "KUN",
        "kokoberavirus": "KOK",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "barmahforestvirus": "BFV",
        "sindbisvirus": "SINDV",
        "umbrevirus": "LSV",
    }
    virus_df = virus_df.with_columns(
        pl.col("virus").replace(virus_acronyms).alias("virus_acronym"),
        pl.col("position").cast(pl.Int64),
        pl.col("coverage").cast(pl.Float64),
    )
    sample_sizes = virus_df.group_by("virus_acronym").agg(
        pl.col("sample").n_unique().alias("n_samples"),
    ).sort(["n_samples", "virus_acronym"], descending=[True, False])
    # Take the median of log10(depth + 1), including explicit zeros.
    # Invert for the axis, which reapplies log10(depth + 1), so the displayed
    # height is exactly the median of the transformed depths.
    median_df = virus_df.group_by(["virus_acronym", "position"]).agg(
        (pl.lit(10.0).pow((pl.col("coverage") + 1).log(10).median()) - 1).alias("median_coverage"),
    ).sort(["virus_acronym", "position"])

    virus_grid = subplot_spec.subgridspec(
        nrows=sample_sizes.height, ncols=2,
        width_ratios=[0.025, 0.975], wspace=0.35, hspace=0.4,
    )
    max_position = virus_df["position"].max()
    shared_ax = None
    coverage_axes = []
    for row, (virus, n_samples) in enumerate(sample_sizes.iter_rows()):
        bracket_ax = fig.add_subplot(virus_grid[row, 0])
        bracket_position = bracket_ax.get_position()
        bracket_ax.set_position([
            bracket_position.x0 - 0.040, bracket_position.y0,
            bracket_position.width, bracket_position.height,
        ])
        ax = fig.add_subplot(virus_grid[row, 1], sharex=shared_ax)
        coverage_axes.append(ax)
        if shared_ax is None:
            shared_ax = ax
        coverage_df = virus_df.filter(pl.col("virus_acronym") == virus)
        for sample_df in coverage_df.partition_by("sample"):
            sample_df = sample_df.sort("position")
            ax.plot(
                sample_df["position"], sample_df["coverage"],
                color="#BDBDBD", alpha=0.5, linewidth=0.4,
            )
        median_coverage = median_df.filter(pl.col("virus_acronym") == virus)
        ax.plot(
            median_coverage["position"], median_coverage["median_coverage"],
            color="black", linewidth=1.2,
        )
        genome_length = coverage_df["position"].max()
        for boundary in [0, genome_length]:
            ax.axvline(boundary, color="#A62635", linestyle="--", linewidth=1, clip_on=False)

        style_coverage_yaxis(ax)
        ax.set_xlim(0, max_position)
        ax.set_xticks(np.arange(0, max_position + 1, 3_000))
        ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_position(("outward", 5))
        ax.spines["bottom"].set_bounds(0, genome_length)
        ax.spines["bottom"].set_color("black")
        ax.spines["bottom"].set_linewidth(0.5)
        ax.tick_params(axis="both", which="both", width=0.6, length=3, color="#555555")
        ax.tick_params(axis="y", which="major", labelsize=8, pad=1)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        if row == sample_sizes.height - 1:
            ax.tick_params(axis="x", bottom=True, labelbottom=True, labelsize=9)
            ax.set_xlabel("Genome position (bp)", fontsize=10)

        label = f"*{virus}" if virus in ["SINDV", "LSV"] else virus
        bracket_ax.text(-0.35, 0.68, label, ha="right", va="center", fontsize=9,
                        transform=bracket_ax.transAxes)
        bracket_ax.text(-0.3, 0.22, f"(n={n_samples})", ha="right", va="center",
                        fontsize=9, fontstyle="italic", transform=bracket_ax.transAxes)
        bracket_ax.set_xticks([])
        bracket_ax.set_yticks([])
        bracket_ax.spines["right"].set_visible(False)
        for side in ["left", "top", "bottom"]:
            bracket_ax.spines[side].set_linewidth(0.8)
            bracket_ax.spines[side].set_color("#555555")
    return coverage_axes


def plot_enrichment(mitochondria_coverage_file: str, virus_coverage_file: str, output_file: str) -> None:
    fig = plt.figure(
        figsize=(A4_WIDTH_INCHES/2, 6.5),
        constrained_layout=False,
    )
    outer_grid = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[0.32, 0.68],
        left=0.30,
        right=0.96,
        top=0.94,
        bottom=0.09,
        hspace=0.18,
    )
    mitochondrial_spec = outer_grid[0, 0]
    virus_spec = outer_grid[1, 0]

    mito_df = pl.read_csv(mitochondria_coverage_file, separator=",")
    virus_df = pl.read_csv(virus_coverage_file, separator="\t")
    plot_mitochondrial_coverage(
        mito_df, fig, mitochondrial_spec, zoom_end=4_500, dashed_boundaries=True,
    )
    plot_virus_coverage_paneled(virus_df, fig, virus_spec)
    fig.text(0.025, 0.55, "Read depth", rotation="vertical", fontsize=10,
             ha="center", va="center")
    fig.savefig(output_file, dpi=600)
    fig.savefig(output_file.replace(".png", ".pdf"), dpi=600)
    plt.close(fig)


def plot_enrichment_and_nc2(
    mitochondria_coverage_file: str,
    virus_coverage_file: str,
    metadata_file: str,
    virus_nc2_file: str,
    output_file: str,
    fit_bin_medians: bool = True,
) -> None:

    fig = plt.figure(figsize=(A4_WIDTH_INCHES, 6.5), constrained_layout=False)
    fig_spec = fig.add_gridspec(
        nrows=2, ncols=2, width_ratios=[0.6, 0.4], height_ratios=[0.3, 0.7],
        left=0.18, right=0.96, top=0.94, bottom=0.09,
        wspace=0.3, hspace=0.3,
    )
    mos_mito_spec = fig_spec[0, 0]
    virus_spec = fig_spec[1, :]
    nc2_spec = fig_spec[0, 1]

    # plot mosquito mito cov
    mito_df = pl.read_csv(mitochondria_coverage_file, separator=",")
    virus_df = pl.read_csv(virus_coverage_file, separator="\t")
    mitochondrial_ax = plot_mitochondrial_coverage(
        mito_df, fig, mos_mito_spec,
        zoom_end=4_500, dashed_boundaries=True,
    )
    mitochondrial_ax.set_xlabel("Genome position (bp)", fontsize=9)

    # plot virus cov
    virus_axes = plot_virus_coverage_paneled(virus_df, fig, virus_spec)
    # Both rows start at the same data-axis edge, despite different grid widths.
    coverage_left = mitochondrial_ax.get_position().x0
    for ax in virus_axes:
        bounds = ax.get_position()
        ax.set_position([
            coverage_left, bounds.y0, bounds.x1 - coverage_left, bounds.height,
        ])
    fig.text(
        0.025, 0.55, "Read depth", rotation="vertical", fontsize=10,
        ha="center", va="center",
    )
    virus_axes[-1].set_xlabel("Genome position (bp)", fontsize=9)

    # plot nc2 curve
    scatter_ax, _ = plot_nc2_curve_panels(
        metadata_file, virus_nc2_file, fig,
        nc2_spec, None, output_file,
        fit_bin_medians=fit_bin_medians,
    )
    scatter_ax.set_xlabel("Cycle threshold in original sample", fontsize=9)
    scatter_ax.set_ylabel("Virus genome coverage", fontsize=9)
    # Increase every label, annotation and tick by one point for this figure.
    fig.canvas.draw()
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize() + 1)

    fig.savefig(output_file, dpi=600)
    fig.savefig(output_file.replace(".png", ".pdf"), dpi=600)
    plt.close(fig)

def quick_plot_of_nc2_of_cox1_across_traps(nc2_file: str, metadata_file: str, mosquito_abundance_results_file: str, output_file: str) -> None:
    fig, ax = plt.subplots(figsize=(3, 3.2), constrained_layout=True)
    # parse files
    nc2_df = pl.read_csv(nc2_file, separator="\t")
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    mosquito_abundance_df = pl.read_csv(mosquito_abundance_results_file)
    # filter nc2 to just metadata samples
    nc2_df = nc2_df.join(metadata_df.select("sample"), left_on="sampleid", right_on="sample", how="inner")
    # if sample in mosquito_abundance_df, add column where "called" is True, else False
    nc2_df = nc2_df.with_columns(
        pl.when(pl.col("sampleid").is_in(mosquito_abundance_df.select("sample").to_series()))
        .then(True)
        .otherwise(False)
        .alias("mosquito_abundance_called")
    )
    # print totals in each group
    print(nc2_df['mosquito_abundance_called'].value_counts())
    # print median and range of mosquitocox1 nc2 in each group
    for group in [True, False]:
        group_df = nc2_df.filter(pl.col("mosquito_abundance_called") == group)
        median = group_df.select(pl.col("mosquitocox1").median()).item()
        min_value = group_df.select(pl.col("mosquitocox1").min()).item()
        max_value = group_df.select(pl.col("mosquitocox1").max()).item()
        print(f"Group {group}: median={median}, range=({min_value}, {max_value})")
    # if mosquitocox1 nc2 value is NaN, replace with 0
    nc2_df = nc2_df.with_columns(
        pl.when(pl.col("mosquitocox1").is_null())
        .then(0)
        .otherwise(pl.col("mosquitocox1"))
        .alias("mosquitocox1")
    )
    # plot box plots of mosquitocox1 nc2 by mosquito_abundance_called, with swarmplot overlayed
    plot_data = nc2_df.to_pandas()
    sns.boxplot(
        data=plot_data,
        x="mosquito_abundance_called",
        y="mosquitocox1",
        ax=ax,
        order=[False, True],
        palette=["#A6CEE3", "#1F78B4"],
        width=0.3,
        linewidth=1,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.5},
    )
    sns.swarmplot(
        data=plot_data,
        x="mosquito_abundance_called",
        y="mosquitocox1",
        ax=ax,
        order=[False, True],
        color=".25",
        size=5,
        alpha=0.6,
    )
    ax.set_xticks([0, 1], labels=["Not called", "Called"])
    ax.set_xlabel("Mosquito abundance", fontsize=10, labelpad=8)
    ax.set_ylabel("Proportion of COX1 covered by\n≥2 reads", fontsize=10, labelpad=8)
    ax.tick_params(axis="both", labelsize=9)
    ax.tick_params(axis="x", length=0)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.9", linewidth=0.6)
    ax.xaxis.grid(False)
    sns.despine(ax=ax)
    fig.savefig(output_file, dpi=600, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
