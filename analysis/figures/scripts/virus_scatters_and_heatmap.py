import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from ..virus.plot import prepare_virus_overview_scatter_data, style_unique_read_count_axis, plot_line_of_best_fit_to_ax
from ..virus.process import combine_ct_dedup_and_nc2_data
from importlib import import_module
import pickle
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pathlib import Path
from ..common.mosquito_ids import shorted_genus
from scipy.optimize import curve_fit
from ..common.regression import fit_linear_regression, log10_reads, inverse_log10_reads
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


abundance_validation = import_module("analysis.abundance-validation.plot")

A4_WIDTH_INCHES = 8.27

def plot_virus_scatter(metadata_file: str, nc2_file: str, dedup_file: str, fig, spec, stats_output: str, legend_ax) -> dict:
    ax = fig.add_subplot(spec)
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
    nc2_groups = [(0, 0.05), (0.05, 0.50), (0.50, 0.95), (0.95, 1)]
    nc2_labels = ["<5%", "5-50%", "50-95%", "≥95%"]
    # genome_coverage_groups = [(0, 0.3), (0.3, 1)]
    # genome_coverage_labels = ["< 30%", '≥ 30%']
    nc2_colours = [
            "#FDD49E",  # light orange
            "#FC8D59",  # orange
            "#E34A33",  # red-orange
            "#B30000",  # dark red
        ]
    # make a combined df
    df = combine_ct_dedup_and_nc2_data(metadata_file, nc2_file, dedup_file, virus_acronyms, nc2_groups, nc2_labels)
    (
        completeness_data,
        completeness_specs,
        stack_order,
        virus_order,
        virus_marker_map,
        completeness_palette,
    ) = prepare_virus_overview_scatter_data(
        df,
        nc2_groups,
        nc2_labels,
        nc2_colours,
    )
    # Missing viral read/coverage records represent non-detections in this panel.
    completeness_data = completeness_data.with_columns(
        pl.col("nc2", "dedup_reads").fill_null(0)
    )
    nc2_cmap = sns.blend_palette(nc2_colours, as_cmap=True)
    sns.scatterplot(
        data=completeness_data.to_pandas(),
        x="ct_value",
        y="dedup_reads",
        hue="nc2",
        hue_norm=(0, 1),
        # hue_order=stack_order,
        style="virus",
        style_order=virus_order,
        markers=virus_marker_map,
        # palette=completeness_palette,
        palette=nc2_cmap,
        s=34,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.25,
        legend=False,
        ax=ax,
    )
    # Match log(reads + 1) while labelling ticks in original read counts.
    ax.set_yscale("function", functions=(lambda values: np.log10(np.maximum(values, -0.999999) + 1), inverse_log10_reads))
    ax.set_ylim(-0.15, max(1, completeness_data["dedup_reads"].max()) * 1.2)
    ticks = [0, 10, 100, 1_000, 10_000, 100_000, 1_000_000]
    ticks = [value for value in ticks if value <= ax.get_ylim()[1]]
    ax.set_yticks(ticks, labels=[f"{value:,}" for value in ticks])
    ax.minorticks_off()
    ax.invert_xaxis()
    # plot line of best fit onto scatter
    plot_line_of_best_fit_to_ax(
        completeness_data,
        x="ct_value",
        y="dedup_reads",
        ax=ax,
        logx=False,
        logy=True,
        logy_offset=1,
        stats_output_file=stats_output,
    )
    ax.set_xlabel("Cycle threshold in original sample")
    ax.set_ylabel("Unique read count")
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # completeness_legend = legend_ax.legend(
    #     handles=[
    #         Patch(facecolor=spec["color"], edgecolor="none", label=spec["label"])
    #         for spec in completeness_specs
    #     ],
    #     title="Genome completeness",
    #     loc="upper left",
    #     bbox_to_anchor=(0.01, 0.935),
    #     bbox_transform=fig.transFigure,
    #     ncol=4,
    #     frameon=False,
    #     borderaxespad=0.0,
    #     fontsize=8,
    #     title_fontsize=10,
    #     columnspacing=1.0,
    #     handletextpad=0.4,
    #     alignment="left",
    # )
    legend_ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker=virus_marker_map[virus],
                linestyle="None",
                color="black",
                markerfacecolor="black",
                markersize=7,
                label=virus_acronyms.get(virus, virus),
            )
            for virus in virus_order
        ],
        title="Virus",
        loc="upper left",
        bbox_to_anchor=(0.01, 0.995),
        bbox_transform=fig.transFigure,
        frameon=False,
        borderaxespad=0.0,
        fontsize=8,
        title_fontsize=10,
        columnspacing=1.0,
        handletextpad=0.4,
        ncol=4,
        alignment="left",
    )
    # legend_ax.add_artist(completeness_legend)
    # add a colourbar for the nc2 values under the Virus legend
    sm = plt.cm.ScalarMappable(cmap=nc2_cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    # Explicit figure coordinates keep the bar below the virus legend.
    # [left, bottom, width, height]; edit the title separately below.
    cbar_ax = fig.add_axes([0.018, 0.907, 0.36, 0.015])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar_ax.set_title("Genome coverage", fontsize=10, loc="left", pad=5)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1], labels=["0", "0.25", "0.5", "0.75", "1"])
    cbar.outline.set_visible(False)
    cbar_ax.tick_params(labelsize=8, length=2, pad=2)
    
    return virus_marker_map

def plot_mosquito_scatter(mosquito_df: pl.DataFrame, fig, spec, legend_ax, virus_marker_map: dict) -> None:
    ax = fig.add_subplot(spec)
    species_order = sorted(mosquito_df["species"].unique().to_list())
    available_markers = ["o", "o"]
    species_marker_map = dict(zip(species_order, available_markers))
    colourblind_colours = ["white","black"]
    species_colour_map = dict(zip(species_order, colourblind_colours))
    sns.scatterplot(
        data=mosquito_df,
        x="surveillance_proportion",
        y="estimated_abundance",
        hue="species",
        style="species",
        style_order=species_order,
        markers=species_marker_map,
        hue_order=species_order,
        palette=species_colour_map,
        edgecolor="black",
        s=34,
        alpha=0.82,
        linewidth=0.5,
        clip_on=False,
        legend=False,
        ax=ax,
    )

    ax.set(xlim=(-0.03, 1.03), ylim=(-0.03, 1.03),
           xticks=np.linspace(0, 1, 6), yticks=np.linspace(0, 1, 6))
    ax.plot([0, 1], [0, 1], color="#333333", linestyle="--", linewidth=1, alpha=0.55, zorder=0)
    ax.set_xlabel("Morphological estimated\nabundance")
    ax.set_ylabel("COSMOS estimated\nabundance")
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend_ax.legend(
        handles=[
            Line2D([0], [0], marker=species_marker_map[species], linestyle="None",
                   markerfacecolor=species_colour_map[species],
                   markeredgecolor='black',
                   markersize=7, label=shorted_genus(species))
            for species in species_order
        ],
        title="Mosquito species",
        alignment="right",
        loc="upper right",
        bbox_to_anchor=(0.99, 0.995),
        bbox_transform=fig.transFigure,
        ncol=1,
        frameon=False,
        borderaxespad=0.0,
        prop={"size": 9, "style": "italic"},
        title_fontsize=10,
        columnspacing=1.0,
        handletextpad=0.4,
    )

def plot_virus_scatters_and_heatmap(metadata_file: str, virus_nc2_file: str, virus_dedup_file: str, mosquito_file: str, heatmap_data_file: str, output_file: str) -> None:
    # draw figure
    fig = plt.figure(
        figsize=(A4_WIDTH_INCHES*0.7, 9.5),
        constrained_layout=False,
    )
    outer_grid = fig.add_gridspec(
        nrows=2,
        ncols=4,
        width_ratios=[0.18, 0.32, 0.32, 0.18],
        height_ratios=[0.3, 0.7],
        left=0.2,
        right=0.95,
        top=0.97,
        bottom=0.04,
        wspace=0.0,
        hspace=0.65,
        )
    # Place both scatter plots across the full top row.
    scatter_grid = outer_grid[0, :].subgridspec(
        nrows=2,
        ncols=2,
        height_ratios=[0.85, 1.8],
        wspace=0.525,
        hspace=0.35,
    )
    virus_legend_ax = fig.add_subplot(scatter_grid[0, 0])
    virus_legend_ax.axis("off")
    virus_spec = scatter_grid[1, 0]
    mosquito_legend_ax = fig.add_subplot(scatter_grid[0, 1])
    mosquito_legend_ax.axis("off")
    mosquito_spec = scatter_grid[1, 1]
    abundance_spec = outer_grid[1, 1:3]

    # parse files
    mosquito_df = pl.read_csv(mosquito_file, separator='\t')
    with open(heatmap_data_file, "rb") as handle:
        heatmap_data = pickle.load(handle)

    # plot virus nc2 vs ct scatter
    stats_output_file = "results/publication_figures/tables/nc2_v_ct_line.txt"
    virus_marker_map = plot_virus_scatter(metadata_file, virus_nc2_file, virus_dedup_file, fig, virus_spec, stats_output_file, virus_legend_ax)

    # plot mosquito program vs our estimates scatter
    plot_mosquito_scatter(mosquito_df, fig, mosquito_spec, mosquito_legend_ax, virus_marker_map)

    # plot the abundance validation heatmap
    _, heatmap_axes = abundance_validation.plot_control_abundance_heatmap(
        heatmap_data,
        simulated_mixture_panel="heatmap",
        cmap="Blues",
        fig=fig,
        subplot_spec=abundance_spec,
        show_group_brackets=False,
    )
    # Leave room on the right for species labels angled upwards to the right.
    right_edge = mosquito_spec.get_position(fig).x1 - 0.08
    shift = right_edge - heatmap_axes[0, 0].get_position().x1
    for ax in heatmap_axes[:, 0]:
        bounds = ax.get_position()
        ax.set_position([bounds.x0 + shift, bounds.y0, bounds.width, bounds.height])
    for ax in heatmap_axes[:, 1]:
        ax.set_visible(False)
    # Move the genus strip and shared species labels above the first panel.
    species_ax = heatmap_axes[0, 0]
    genus_ax = heatmap_axes[-1, 0]
    species_bounds = species_ax.get_position()
    genus_bounds = genus_ax.get_position()
    genus_ax.set_position([
        species_bounds.x0, species_bounds.y1 + 0.025,
        species_bounds.width, genus_bounds.height,
    ])
    genus_ax.xaxis.tick_top()
    genus_ax.xaxis.set_label_position("top")
    genus_ax.tick_params(axis="x", labeltop=True, labelbottom=False, length=0, pad=2)
    for label in genus_ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("left")
        label.set_rotation_mode("anchor")
    genus_ax.set_xlabel("COSMOS assigned species", fontsize=11,
                        labelpad=8, fontweight="bold")

    # Position each vertical label beside its panel's row tick labels.
    heatmap_label_x = 0.16
    species_ax.set_ylabel("Morphological\nidentification", fontsize=12,
                          rotation=90, ha="center", va="center")
    species_ax.yaxis.set_label_coords(
        heatmap_label_x, (species_bounds.y0 + species_bounds.y1) / 2,
        transform=fig.transFigure,
    )
    dilution_ax = heatmap_axes[-2, 0]
    bounds = dilution_ax.get_position()
    label_y = (bounds.y0 + bounds.y1) / 2
    dilution_ax.set_ylabel("Variable input", fontsize=12, rotation=90,
                           ha="center", va="center")
    dilution_label_x = heatmap_label_x
    dilution_ax.yaxis.set_label_coords(dilution_label_x, label_y, transform=fig.transFigure)
    for fraction_down, species_label in [
        (0.25, "Cx. siti &\nCx. pipi"),
        (0.75, "Ae. vigi &\nCx. annu"),
    ]:
        fig.text(
            dilution_label_x + 0.13, bounds.y1 - fraction_down * bounds.height,
            species_label, fontsize=10, fontstyle="italic", rotation=0, ha="center", va="center",
        )

    # save figure
    fig.savefig(output_file, dpi=600)
    fig.savefig(output_file.replace(".png", ".pdf"), dpi=600)

def supplementary_nc2_curve_figure(
    metadata_file: str,
    virus_nc2_file: str,
    output_file: str,
    fit_bin_medians: bool = False,
) -> None:
    fig = plt.figure(
        figsize=(A4_WIDTH_INCHES, 3.5),
        constrained_layout=False,
    )
    outer_grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.13,
        right=0.9,
        top=0.88,
        bottom=0.25,
        wspace=0.35,
    )
    scatter_spec = outer_grid[0, 0]
    boxplot_spec = outer_grid[0, 1]

    try:
        plot_nc2_curve_panels(
            metadata_file, virus_nc2_file, fig, scatter_spec, boxplot_spec,
            output_file, fit_bin_medians=fit_bin_medians,
        )
        fig.savefig(output_file, dpi=300)
    finally:
        plt.close(fig)


def plot_nc2_curve_panels(
    metadata_file: str,
    virus_nc2_file: str,
    fig,
    scatter_spec,
    boxplot_spec,
    output_file: str,
    fit_bin_medians: bool = False,
):
    """Draw NC2 panels and sidecars; pass boxplot_spec=None for the curve alone."""
    # parse files
    metadata = pl.read_csv(metadata_file, separator='\t')
    nc2 = pl.read_csv(virus_nc2_file, separator='\t')

    # for every positive ct in a sample across a list of viruses, make a table with nc2 value
    nc2_with_metadata = nc2.join(metadata, left_on="sampleid", right_on="sample", how="right")
    test_viruses = {
        "barmahforestvirus": "BFV",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "kokoberavirus": "KOK",
        "murrayvalleyencephalitisvirus": "MVEV",
        "rossrivervirus": "RRV",
        "stratfordvirus": "STRV",
    }  
    positive_cts = pl.concat([
        nc2_with_metadata.select(
            "sample",
            pl.lit(virus).alias("virus"),
            pl.lit(acronym).alias("virus_acronym"),
            pl.col(f"ct_{acronym}")
            .cast(pl.Float64, strict=False)
            .alias("ct_value"),
            pl.col(virus).cast(pl.Float64).fill_null(0).alias("nc2"),
        )
        for virus, acronym in test_viruses.items()
    ]).filter((pl.col("ct_value") > 0) & pl.col("ct_value").is_finite())

    if positive_cts.is_empty():
        raise ValueError("No positive Ct values found for the selected viruses.")

    # Right-closed bins match the supplied table: (0, 16], (16, 19], etc.
    breaks = [16, 19, 22, 25, 28, 31, 34, 37]
    labels = [f"({low}, {high}]" for low, high in zip([0] + breaks[:-1], breaks)]
    positive_cts = positive_cts.with_columns(
        pl.col("ct_value").cut(breaks, labels=labels + ["(37, ∞)"])
        .cast(pl.String).alias("ct_category")
    )
    if positive_cts["ct_value"].max() > 37:
        labels.append("(37, ∞)")
    summary = pl.DataFrame({"ct_category": labels}).join(
        positive_cts.group_by("ct_category").agg(
            pl.len().alias("count"),
            pl.col("ct_value").median().alias("median_ct"),
            pl.col("nc2").median().alias("median_nc2"),
        ),
        on="ct_category", how="left",
    ).with_columns(pl.col("count").fill_null(0))

    scatter_ax = fig.add_subplot(scatter_spec)
    boxplot_ax = (fig.add_subplot(boxplot_spec, sharey=scatter_ax)
                  if boxplot_spec is not None else None)
    scatter_ax.scatter(positive_cts["ct_value"], positive_cts["nc2"],
                       s=12, color="0.55", alpha=0.45, edgecolors="none")
    populated = summary.filter(pl.col("count") > 0)
    x = populated["median_ct"].to_numpy()
    y = populated["median_nc2"].to_numpy()
    if len(x) < 2:
        raise ValueError("At least two populated Ct bins are needed for the curve.")
    # scatter_ax.plot(x, y, "o-", color="#B30000", markersize=4, linewidth=1.5)

    # add EC50 curve fit to scatter plot
    if fit_bin_medians:
        fit_x = populated["median_ct"].to_numpy()
        fit_y = populated["median_nc2"].to_numpy()
    else:
        data = positive_cts.select(["ct_value", "nc2"]).drop_nulls().to_pandas()
        fit_x = data["ct_value"].to_numpy()
        fit_y = data["nc2"].to_numpy()
    initial = [
        0.0,                       # bottom
        1.0,                       # top
        np.median(fit_x),           # EC50
        1.0,                   # slope
    ]

    parameters, covariance = curve_fit(
        decreasing_logistic,
        fit_x,
        fit_y,
        p0=initial,
        bounds=(
            [0, 0, fit_x.min(), 0.01],
            [1, 1, fit_x.max(), 100],
        ),
        maxfev=20_000,
    )

    bottom, top, midpoint, slope = parameters
    ct50 = midpoint + slope * np.log(
        (top - bottom) / (0.5 - bottom) - 1
    )
    ct95 = (
        midpoint + slope * np.log((top - bottom) / (0.95 - bottom) - 1)
        if bottom < 0.95 < top else None
    )
    rng = np.random.default_rng(42)
    ct50_values = []
    bootstrap_parameters = []
    for _ in range(1000):
        indices = rng.choice(len(fit_x), len(fit_x), replace=True)
        try:
            boot_parameters, _ = curve_fit(
                decreasing_logistic,
                fit_x[indices], fit_y[indices], p0=parameters,
                bounds=([0, 0, fit_x.min(), 0.01], [1, 1, fit_x.max(), 100]),
                maxfev=20_000,
            )
            boot_bottom, boot_top, boot_midpoint, boot_slope = boot_parameters
            bootstrap_parameters.append(boot_parameters)
            if boot_bottom < 0.5 < boot_top:
                boot_ct50 = boot_midpoint + boot_slope * np.log(
                    (boot_top - boot_bottom) / (0.5 - boot_bottom) - 1
                )
                if np.isfinite(boot_ct50):
                    ct50_values.append(boot_ct50)
        except (RuntimeError, ValueError):
            pass
    if len(ct50_values) < 20:
        raise ValueError("Too few successful bootstrap fits to estimate a Ct50 CI.")
    ci_low, ci_high = np.percentile(ct50_values, [2.5, 97.5])
    print(f"Ct50 = {ct50:.2f} (95% CI {ci_low:.2f}–{ci_high:.2f})")

    ct_95_values = []
    for boot_parameters in bootstrap_parameters:
        boot_bottom, boot_top, boot_midpoint, boot_slope = boot_parameters
        if boot_bottom < 0.95 < boot_top:
            try:
                boot_ct95 = boot_midpoint + boot_slope * np.log(
                    (boot_top - boot_bottom) / (0.95 - boot_bottom) - 1
                )
                if np.isfinite(boot_ct95):
                    ct_95_values.append(boot_ct95)
            except ValueError:
                pass
    if len(ct_95_values) >= 20:
        ct95_ci_low, ct95_ci_high = np.percentile(ct_95_values, [2.5, 97.5])
        print(f"Ct95 = {ct95:.2f} (95% CI {ct95_ci_low:.2f}–{ct95_ci_high:.2f})")
    else:
        ct95_ci_low, ct95_ci_high = None, None
        print("Too few successful bootstrap fits to estimate a Ct95 CI.")

    ct_grid = np.linspace(10, 40, 400)
    fit = decreasing_logistic(ct_grid, *parameters)
    bootstrap_curves = np.asarray([
        decreasing_logistic(ct_grid, *boot_parameters)
        for boot_parameters in bootstrap_parameters
    ])
    curve_ci_low, curve_ci_high = np.percentile(
        bootstrap_curves, [2.5, 97.5], axis=0
    )

    scatter_ax.plot(
        ct_grid,
        fit,
        color="black",
        linewidth=2,
        label=f"Logistic fit, Ct EC50 = {ct50:.1f}",
    )
    # scatter_ax.fill_between(
    #     ct_grid,
    #     curve_ci_low,
    #     curve_ci_high,
    #     color="firebrick",
    #     alpha=0.18,
    #     linewidth=0,
    #     label="95% bootstrap CI",
    # )
    scatter_ax.axhline(0.95, color="0.4", linestyle="--", linewidth=0.8, zorder=0)
    scatter_ax.axvline(ct50, color="black", linestyle=":")
    scatter_ax.annotate(
        rf"Ct$_{{50}}$ ≈ {ct50:.1f}",
        xy=(1, 0.5), xycoords=scatter_ax.get_yaxis_transform(),
        xytext=(-4, 4), textcoords="offset points",
        fontsize=8, ha="right", va="bottom",
    )
    if ct95 is not None:
        scatter_ax.axvline(ct95, color="black", linestyle=":")
        scatter_ax.annotate(
            rf"Ct$_{{95}}$ ≈ {ct95:.1f}",
            xy=(1, 0.95), xycoords=scatter_ax.get_yaxis_transform(),
            xytext=(-4, -4), textcoords="offset points",
            fontsize=8, ha="right", va="top",
        )
    crossing_text = f"Ct EC50 = {ct50:.2f} (95% CI {ci_low:.2f}–{ci_high:.2f})"
    scatter_ax.invert_xaxis()
    scatter_ax.set(xlabel="Cycle threshold (Ct)",
                   ylabel="Genome coverage\n(proportion ≥2×)")

    if boxplot_ax is not None:
        sns.boxplot(data=positive_cts.to_pandas(), x="ct_category", y="nc2",
                    order=labels, color="0.35", width=0.65, linewidth=0.8,
                    boxprops={"facecolor": "none"}, showfliers=False, ax=boxplot_ax)
        jitter_rng = np.random.default_rng(42)
        for i, category in enumerate(labels):
            values = positive_cts.filter(pl.col("ct_category") == category)["nc2"].to_numpy()
            boxplot_ax.scatter(i + jitter_rng.uniform(-0.2, 0.2, len(values)), values,
                               s=14, facecolors="none", edgecolors="0.35", alpha=0.75,
                               linewidths=0.6, zorder=3)
        # Display square brackets as requested; bin assignment remains right-closed.
        boxplot_ax.set_xticks(range(len(labels)), [label.replace("(", "[").replace(")", "]") for label in labels])
        for i, count in enumerate(summary["count"]):
            boxplot_ax.text(
                i, 1.03, f"n={count}", ha="left", va="bottom", fontsize=8,
                rotation=45, rotation_mode="default",
            )
        boxplot_ax.set(xlabel="Ct category", ylabel="")
        boxplot_ax.tick_params(axis="x", rotation=45, labelsize=8)
        boxplot_ax.invert_xaxis()
        for label in boxplot_ax.get_xticklabels():
            label.set_ha("right")
    for ax in ([scatter_ax, boxplot_ax] if boxplot_ax is not None else [scatter_ax]):
        ax.axhline(0.5, color="0.4", linestyle="--", linewidth=0.8, zorder=0)
        ax.set_ylim(-0.04, 1.12)
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.grid(axis="y", alpha=0.2)
        sns.despine(ax=ax, trim=True, left=False, bottom=False)

    scatter_ax.set_xlim(40, 10)

    # save figure
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    summary.write_csv(Path(output_file).with_suffix(".summary.tsv"), separator="\t")
    pl.DataFrame({"median_ct": x, "median_nc2": y}).write_csv(
        Path(output_file).with_suffix(".curve.tsv"), separator="\t"
    )
    Path(output_file).with_suffix(".threshold.txt").write_text(
        f"NC2 input: {virus_nc2_file}\n{crossing_text}\n"
        f"Logistic parameters:\n"
        f"  bottom = {bottom:.6f} (SE {np.sqrt(covariance[0, 0]):.6f})\n"
        f"  top = {top:.6f} (SE {np.sqrt(covariance[1, 1]):.6f})\n"
        f"  midpoint parameter = {midpoint:.6f} (SE {np.sqrt(covariance[2, 2]):.6f})\n"
        f"  Ct at 50% coverage = {ct50:.6f}\n"
        f"  Ct at 95% coverage = {ct95 if ct95 is not None else 'not reached by fitted curve'}\n"
        f"  Ct50 95% CI = {ci_low:.6f}–{ci_high:.6f}\n"
        f"  slope = {slope:.6f} (SE {np.sqrt(covariance[3, 3]):.6f})\n"
        f"Fit data: {'Ct-bin medians' if fit_bin_medians else 'individual sample–virus pairs'}.\n"
        "Pooled positive sample–virus pairs; missing NC2 values are zero.\n"
        "Curve: joined median coverage points at observed median Ct in each bin.\n"
        "Ct bins remain right-closed (lower, upper], despite square-bracket display labels.\n"
        "This is a descriptive estimate, not a diagnostic cutoff.\n"
    )
    return scatter_ax, boxplot_ax

def decreasing_logistic(ct, bottom, top, ec50, slope):
    return bottom + (top - bottom) / (
        1 + np.exp((ct - ec50) / slope)
    )

def association_between_abundance_and_log_read_count(data: pl.DataFrame, abundance_col: str, read_count_col: str, perform_log=True) -> tuple:
    filtered = data.select(abundance_col, read_count_col).drop_nulls()
    reads = filtered[read_count_col].to_numpy()
    response = log10_reads(reads) if perform_log else reads
    fit = fit_linear_regression(filtered[abundance_col].to_numpy(), response)
    return fit.slope, fit.intercept, fit.rvalue, fit.pvalue, fit.stderr, fit.n


def ISV_heatmap_figure(
    heatmap_values,
    metadata_df,
    sample_col="sample",
    longitude_col="longitude",
    figsize=(A4_WIDTH_INCHES, 6),
    output_path=None,
    ):
    # Remove samples with no signal across any selected virus
    heatmap_values = heatmap_values.loc[
        :,
        heatmap_values.max(axis=0) > 0
    ]
    # Get west -> east sample order from metadata
    sample_order = (
        metadata_df
        .select([sample_col, longitude_col])
        .drop_nulls(longitude_col)
        .unique(subset=[sample_col])
        .sort(longitude_col)
        [sample_col]
        .to_list()
    )
    # Retain only samples present in heatmap
    sample_order = [
        sample
        for sample in sample_order
        if sample in heatmap_values.columns
    ]
    heatmap_values = heatmap_values[sample_order]
    # Cluster virus rows only
    if heatmap_values.shape[0] > 1:
        row_distances = pdist(
            heatmap_values.values,
            metric="euclidean"
        )

        row_linkage = linkage(
            row_distances,
            method="average"
        )
        row_order = leaves_list(row_linkage)
        heatmap_values = heatmap_values.iloc[row_order, :]
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        heatmap_values,
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        xticklabels=False,
        yticklabels=True,
        linewidths=0,
        cbar_kws={
            "label": "Proportion of genome covered ≥2×"
        },
        cbar=False,
        ax=ax
    )
    ax.set_xlabel("Surveillance samples (west → east)")
    ax.set_ylabel("Virus")
    ax.tick_params(
        axis="y",
        labelrotation=0,
        labelsize=10
    )
    fig.tight_layout()
    if output_path is not None:
        fig.savefig(
            output_path,
            dpi=300,
            bbox_inches="tight"
        )
    return fig, ax, heatmap_values

def wuhan6_supp_figure():
    ## PARSE IN ALL THE DATA
    COSMOS1 = f"/Volumes/research-data/PRJ-CSM/DATA/2023-03-17_COSMOS1/castanet-insectvirs-output/COSMOS1_depth.csv"
    COSMOS2 = f"/Volumes/research-data/PRJ-CSM/DATA/2023-06-22_COSMOS2/castanet-insectvirs-output/COSMOS2_depth.csv"
    COSMOS3 = f"/Volumes/research-data/PRJ-CSM/DATA/2023-11-01_COSMOS3/castanet-insectvirs-output/COSMOS3_depth.csv"
    COSMOS4 = "/Volumes/research-data/PRJ-CSM/DATA/2024-05-17_COSMOS4/castanet-insectvirs-output/COSMOS4_depth_with_clin.csv"
    # Read the per-sample depth summaries for all four batches.
    df1 = pl.read_csv(COSMOS1, separator=",")
    df2 = pl.read_csv(COSMOS2, separator=",")
    df3 = pl.read_csv(COSMOS3, separator=",")
    df4 = pl.read_csv(COSMOS4, separator=",")
    # Batch exports have different columns; align by name and fill absent fields with null.
    df = pl.concat([df1, df2, df3, df4], how="diagonal_relaxed")
    # filter samples to analysis set 
    metadata = "/Users/afor0966/Documents/code/COSMOS-analysis/data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    metadata_df = pl.read_csv(metadata, separator="\t")
    df = df.join(metadata_df, left_on="sampleid", right_on="sample", how="inner", validate="m:1")
    # print total number of probetypes with more than 0.5 in prop_npos_cov2 column in at least one sample
    probetypes = df.filter(pl.col("prop_npos_cov2") > 0.5).select("probetype").unique().to_series().to_list()
    print(f"Total number of probetypes with >0.5 prop_npos_cov2 in at least one sample: {len(probetypes)}")
    # print top 5 most common probetypes with >0.5 prop_npos_cov2 in at least one sample
    top5_probetypes = (
        df
        .filter(pl.col("prop_npos_cov2") > 0.5)
        .group_by("probetype")
        .agg(
            pl.col("sampleid").n_unique().alias("count")
        )
        .sort("count", descending=True)
        .head(5)
    )
    print(f"Top 5 most common probetypes with >0.5 prop_npos_cov2 in at least one sample:\n{top5_probetypes}")

    ## ISV HEATMAP
    selected_probetypes = (
        df
        .filter(pl.col("prop_npos_cov2") > 0.5)
        .select("probetype")
        .unique()
        .to_series()
        .to_list()
    )
    # Retain all coverage values for those probetypes
    heatmap_data = (
        df
        .filter(pl.col("probetype").is_in(selected_probetypes))
        .group_by(["probetype", "sampleid"])
        .agg(
            pl.col("prop_npos_cov2")
            .max()
            .alias("prop_npos_cov2")
        )
    )
    # Pivot to viruses x samples
    heatmap_pivot = (
        heatmap_data
        .pivot(
            index="probetype",
            on="sampleid",
            values="prop_npos_cov2"
        )
        .fill_null(0)
    )
    heatmap_values = (
        heatmap_pivot
        .to_pandas()
        .set_index("probetype")
    )
    ISV_heatmap_figure(
        heatmap_values,
        metadata_df,
        sample_col="sample",
        longitude_col="longitude",
        output_path=(
            "/Users/afor0966/Documents/code/COSMOS-analysis/"
            "results/publication_figures/plots/ISVs_heatmap.png"
        )
    )

    ## WU6 ANALYSIS
    abundance_data = "/Users/afor0966/Documents/code/COSMOS-analysis/results/publication_figures/data/combined_abundance_data.tsv"
    abundance_df = pl.read_csv(abundance_data)
    # filter df to where probetype is Wuhan-Mosquito-Virus
    wu6_df = df.filter(pl.col("probetype") == "Wuhan-Mosquito-Virus")
    # Keep every analysis sample with an abundance estimate, including no Wuhan detection.
    analysis_samples = metadata_df.select("sample").join(
        abundance_df, on="sample", how="inner", validate="1:1"
    ).rename({"sample": "sampleid"})
    joined = analysis_samples.join(
        wu6_df.select("sampleid", "n_reads_dedup", "prop_npos_cov2"),
        on="sampleid", how="left", validate="1:1",
    ).with_columns(
        pl.col("n_reads_dedup", "prop_npos_cov2").fill_null(0)
    )
    plot_columns = ["Culex australicus", "n_reads_dedup", "prop_npos_cov2"]
    # scatter plot of prop_npos_cov2 vs estimated_abundance
    fig, ax = plt.subplots(figsize=(4, 3.5))
    colour_map = sns.color_palette("flare", as_cmap=True)
    # replace null values in prop_npos_cov2 and n_reads_dedup with 0
    joined = joined.with_columns(
        pl.col("n_reads_dedup", "prop_npos_cov2").fill_null(0)
    )
    scatter = sns.scatterplot(
        data=joined.to_pandas(),
        y="n_reads_dedup",
        x="Culex australicus",
        ax=ax,
        hue="prop_npos_cov2",
        hue_norm=(0, 1),
        palette=colour_map,
        edgecolor="black",
        s=30,
        alpha=0.8,
        legend=False,
    )
    rendered_points = sum(
        len(collection.get_offsets())
        for collection in scatter.collections
    )
    assert rendered_points == len(joined), "Mismatch between scatter points and data rows."
    scalar_mappable = plt.cm.ScalarMappable(
        norm=plt.Normalize(vmin=0, vmax=1), cmap=colour_map
    )
    scalar_mappable.set_array([])
    colourbar_ax = fig.add_axes([0.28, 0.9, 0.44, 0.025])
    colourbar = fig.colorbar(
        scalar_mappable, cax=colourbar_ax, orientation="horizontal"
    )
    colourbar_ax.set_title(
        "Genome coverage (proportion ≥2×)", fontsize=8, loc="left", pad=3
    )
    colourbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
    colourbar.outline.set_visible(False)
    colourbar_ax.tick_params(labelsize=7, length=2, pad=2)
    # Fit log(1 + reads) so every sample, including zero counts, contributes.
    x = joined["Culex australicus"].to_numpy()
    log_reads = log10_reads(joined["n_reads_dedup"].to_numpy())
    fit = fit_linear_regression(x, log_reads)
    n = fit.n
    x_vals, fitted_log_reads, lower, upper = fit.predict()
    # Pointwise mean confidence interval on the log10(reads + 1) scale.
    ax.fill_between(
        x_vals, inverse_log10_reads(lower), inverse_log10_reads(upper),
        color="black", alpha=0.15, linewidth=0, label="95% confidence band", zorder=0,
    )
    y_vals = inverse_log10_reads(fitted_log_reads)
    ax.plot(x_vals, y_vals, color="black", linestyle="-", linewidth=1.5, label="Linear regression")
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    ax.set_ylabel("Unique read count of\nWuhan Mosquito Virus", fontsize=10)
    ax.set_xlabel(r"Estimated abundance of $\mathit{Culex\ australicus}$", fontsize=10)
    # change fontsoze of tck labels
    ax.tick_params(axis="both", which="major", labelsize=8)
    # Match the model's transformation, retaining a genuine zero tick.
    ax.set_yscale("function", functions=(lambda values: np.log10(np.maximum(values, -0.999999) + 1), inverse_log10_reads))
    ax.set_ylim(-0.15, max(1100, joined["n_reads_dedup"].max() * 1.1, inverse_log10_reads(upper).max() * 1.1))
    ax.set_yticks([0, 1, 10, 100, 1000], labels=["0", "1", "10", "100", "1,000"])
    ax.minorticks_off()
    sns.despine(ax=ax, trim=True, left=False, bottom=False)
    ax.grid(axis="both", alpha=0.2, linewidth=0.8)
    # make spines black and have ticks
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", length=4, width=0.8, color="black")
    fig.tight_layout(rect=[0, 0, 1, 0.87])
    fig.savefig("/Users/afor0966/Documents/code/COSMOS-analysis/results/publication_figures/plots/wuhan6_supp_figure.png", dpi=300)
    # save stats into a text file
    rsquared = fit.rvalue ** 2
    with open("/Users/afor0966/Documents/code/COSMOS-analysis/results/publication_figures/plots/wuhan6_supp_figure_stats.txt", "w") as f:
        f.write(f"Log10(read_count + 1) linear regression: slope={fit.slope:.4f}, intercept={fit.intercept:.4f}, r_value={fit.rvalue:.4f}, p_value={fit.pvalue:.4e}, std_err={fit.stderr:.4f}, n={n}\n")
        f.write(f"R-squared: {rsquared:.4f}\n")
        f.write("response_transform\tlog10(read_count + 1)\n")
        for name, value in fit.summary().items():
            f.write(f"{name}\t{value}\n")
    f.close()


def bfv_supp_figure(
    input_file="data/BFVculture-dilutions/bfv_expt_COSMOS3.csv",
    output_file="results/publication_figures/plots/bfv_supp_figure.png",
):
    """Standalone BFV dilution figure with OLS and Student-t mean intervals."""
    data = pl.read_csv(input_file)
    fit = fit_linear_regression(data["ct"].to_numpy(), log10_reads(data["n_reads_dedup"].to_numpy()))
    x, mean, lower, upper = fit.predict()
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with sns.axes_style("ticks"), sns.plotting_context("paper", font_scale=1.1):
        fig, ax = plt.subplots(figsize=(4.2, 3.8))
        ax.fill_between(x, inverse_log10_reads(lower), inverse_log10_reads(upper),
                        color="black", alpha=0.15, linewidth=0, label="95% confidence band")
        ax.plot(x, inverse_log10_reads(mean), color="black", linewidth=1.5, label="OLS fit")
        ax.scatter(fit.x, inverse_log10_reads(fit.y), s=42, color="black",
                   edgecolor="white", linewidth=0.6, zorder=3)
        ax.set_yscale("function", functions=(lambda values: np.log10(np.maximum(values, -0.999999) + 1), inverse_log10_reads))
        ax.set_ylim(0, max(data["n_reads_dedup"].max(), inverse_log10_reads(upper).max()) * 1.5)
        ax.set_yticks([0, 10, 100, 1_000, 10_000], labels=["0", "10", "100", "1,000", "10,000"])
        ax.set_xticks([20, 22, 24, 26, 28, 30])
        ax.minorticks_off()
        ax.set(xlabel="Cycle threshold (Ct)", ylabel="Unique read count")
        ax.grid(alpha=0.2, linewidth=0.7)
        sns.despine(ax=ax)
        fig.tight_layout()
        fig.savefig(output, dpi=300)
        fig.savefig(output.with_suffix(".pdf"))
        plt.close(fig)
    summary = dict(input_file=input_file, response_transform="log10(read_count + 1)",
                   formula=f"log10(n_reads_dedup + 1) = {fit.intercept:.8g} + {fit.slope:.8g} * ct",
                   **fit.summary())
    output.with_suffix(".stats.tsv").write_text("".join(f"{key}\t{value}\n" for key, value in summary.items()))
    output.with_suffix(".caption.txt").write_text(
        "Barmah Forest virus (BFV) dilution series. Points show deduplicated read counts "
        "against cycle threshold (Ct) for four samples. The line is an ordinary least "
        "squares regression of log10(read count + 1) on Ct, including an intercept. "
        "Shading indicates the pointwise 95% confidence interval for the fitted mean "
        "on the transformed scale, calculated using Student-t critical values with "
        "2 residual degrees of freedom. Predictions and interval endpoints are "
        "back-transformed; the vertical axis uses log10(read count + 1) spacing "
        "with tick labels in original read counts.\n"
    )
