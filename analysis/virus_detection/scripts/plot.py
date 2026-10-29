# plot virus results from castanet

import polars as pl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
from analysis.common.regression import fit_linear_regression, log10_reads, inverse_log10_reads
from analysis.virus_detection.scripts.process import generate_matrix
import matplotlib.gridspec as gridspec

def filter_pcr_positive(data: pl.DataFrame) -> pl.DataFrame:
    return data.filter(pl.col("ct_value").cast(pl.Float64, strict=False).fill_null(0) > 0)


def filter_dedup_read_positive(data: pl.DataFrame) -> pl.DataFrame:
    return data.filter(pl.col("dedup_reads").cast(pl.Float64, strict=False).fill_null(0) > 0)


def plot_boxplot_with_hue(data: pl.DataFrame, x: str, y: str, hue: str, title: str, ax: plt.Axes) -> None:
    sns.boxplot(data=data.to_pandas(), x=x, y=y, hue=hue, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)

def plot_scatter_on_ax(data: pl.DataFrame, x: str, y: str, title: str, ax: plt.Axes, hue: str = None, style: str = None) -> None:
    sns.scatterplot(data=data.to_pandas(), x=x, y=y, hue=hue, style=style, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)

def plot_line_of_best_fit_to_ax(data: pl.DataFrame, x: str, y: str, ax: plt.Axes, logx: bool=True, logy: bool=True, stats_output_file: str | None = None, logy_offset: float = 1) -> None:
    plot_data = data.select([x, y]).drop_nulls().to_pandas()
    if logx:
        plot_data = plot_data[plot_data[x] > 0]
    x_values = plot_data[x].to_numpy(dtype=float)
    y_values = plot_data[y].to_numpy(dtype=float)
    fit_x = np.log10(x_values) if logx else x_values
    if logy:
        if logy_offset == 1:
            fit_y = log10_reads(y_values)
        else:
            raise ValueError("Read-count regression requires log10(reads + 1).")
    else:
        fit_y = y_values
    fit = fit_linear_regression(fit_x, fit_y)
    fit_x_line, fit_y_line, lower_y_line, upper_y_line = fit.predict()
    x_line = 10 ** fit_x_line if logx else fit_x_line
    if logy:
        fit_y_line, lower_y_line, upper_y_line = map(
            inverse_log10_reads, (fit_y_line, lower_y_line, upper_y_line)
        )
    ax.plot(x_line, fit_y_line, color="tab:grey", linewidth=2)
    ax.fill_between(x_line, lower_y_line, upper_y_line, color="tab:grey", alpha=0.2, linewidth=0)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if stats_output_file is not None:
        fit_x_label = f"log10({x})" if logx else x
        fit_y_label = f"log10({y} + 1)" if logy else y
        summary = fit.summary()
        summary["formula"] = f"{fit_y_label} = {fit.intercept:.8g} + {fit.slope:.8g} * {fit_x_label}"
        summary["response_transform"] = "log10(read_count + 1)" if logy else "identity"
        with open(stats_output_file, "w") as f:
            for name, value in summary.items():
                f.write(f"{name}\t{value}\n")


def plot_dedup_reads_vs_ct_value(data: pl.DataFrame, output_file: str, virus_acronyms: dict) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(9, 6))
    plot_scatter_on_ax(data, x="ct_value", y="dedup_reads", title="", ax=ax, hue="nc2", style="virus")
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _: f"{int(value):,}" if value >= 1 else f"{value:g}"))
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlabel("cycle threshold in original sample")
    ax.set_ylabel("unique read count")
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # change legend labels to virus acronyms and move to the outside top right of the plot
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [virus_acronyms.get(label, label) for label in labels]
    ax.legend(handles, new_labels, loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(output_file)

def assign_groups_to_values(value: float, groups: list) -> str:
    for i, (lower, upper) in enumerate(groups):
        if lower <= value < upper or (i == len(groups) - 1 and lower <= value <= upper):
            return f"{lower}-{upper}"
    return "unknown"

def plot_box_plot_of_dedup_reads_by_ct_value_coloured_by_nc2(
    data: pl.DataFrame,
    output_file: str,
    genome_coverage_labels: list[str],
    genome_coverage_colours: list[str],
) -> None:
    plot_data = data.with_columns(
        pl.when(pl.col("genome_coverage_category") == "missed")
        .then(pl.lit("missed"))
        .otherwise(pl.lit("detected"))
        .alias("nc2_group_binary")
    )
    category_order = ["missed", "detected"]
    colour_by_category = {
        "missed": genome_coverage_colours[0],
        "detected": "#FFFFFF",
    }
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(4, 6))
    sns.boxplot(
        data=plot_data.to_pandas(),
        x="nc2_group_binary",
        y="ct_value",
        # hue="nc2_group_binary",
        order=category_order,
        # hue_order=category_order,
        # palette=colour_by_category,
        color="white",
        ax=ax,
        showfliers=False,
        width=0.6,
        linewidth=1.2,
        dodge=False,
        legend=False,
    )
    sns.stripplot(
        data=plot_data.to_pandas(),
        x="nc2_group_binary",
        y="ct_value",
        order=category_order,
        color="tab:gray",
        ax=ax,
        dodge=False,
        jitter=0.22,
        size=10,
        alpha=0.75,
        linewidth=0.2,
        edgecolor="white",
    )
    ax.set_xlabel("virus genome\ncompleteness")
    ax.set_ylabel("virus cycle threshold in\noriginal sample")
    ax.grid(True, which="major", axis="y", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks(range(len(category_order)))
    ax.set_xticklabels(category_order, rotation=35, ha="right", rotation_mode="anchor")
    plt.tight_layout()
    plt.savefig(output_file)

def plot_stacked_bar_plot_for_each_virus_with_missed_incomplete_and_complete_categories(
    data: pl.DataFrame,
    output_file: str,
    nc2_groups: list,
    virus_acronyms: dict,
    nc2_labels: list[str],
    nc2_colours: list[str],
) -> None:
    stack_order = nc2_labels[::-1]
    if len(nc2_groups) != len(nc2_labels):
        raise ValueError("nc2_groups and nc2_labels must have the same length")
    # create a new column for completeness category based on nc2 value
    plot_data = data.with_columns(
        pl.col("nc2").fill_null(0).map_elements(
            lambda nc2: next(
                (
                    nc2_labels[index]
                    for index, (lower, upper) in enumerate(nc2_groups)
                    if lower <= nc2 < upper or (index == len(nc2_groups) - 1 and lower <= nc2 <= upper)
                ),
                "unknown",
            ),
            return_dtype=pl.Utf8,
        ).alias("completeness_category")
    )
    # count the number of samples in each category for each virus
    category_counts = plot_data.group_by(["virus", "completeness_category"]).agg(pl.len().alias("count"))
    # pivot the data to have categories as columns
    pivot_data = category_counts.pivot(index="virus", columns="completeness_category", values="count").fill_null(0)
    # plot a stacked bar plot
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    pivot_data_pd = pivot_data.to_pandas().set_index("virus").reindex(columns=stack_order, fill_value=0)
    sort_col = "complete" if "complete" in pivot_data_pd.columns else stack_order[0]
    pivot_data_pd = pivot_data_pd.sort_values(by=sort_col, ascending=False)
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    stack_colours = [
        nc2_colours[nc2_labels.index(label)]
        for label in stack_order
    ]
    pivot_data_pd.plot(kind="bar", stacked=True, color=stack_colours, ax=ax, legend=False, width=0.82)
    ax.set_xlabel("virus")
    ax.set_ylabel("number of samples")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels([virus_acronyms.get(virus, virus) for virus in pivot_data_pd.index], rotation=45, ha="right")
    legend_handles = [
        Patch(facecolor=nc2_colours[index], edgecolor="none", label=label)
        for index, label in enumerate(nc2_labels)
    ]
    handles, labels = ax.get_legend_handles_labels()
    # add ranges to legend labels
    new_labels = [f"{label} ({lower}-{upper})" if upper != 1000000000 else f"{label} ({lower}+)" for label, (lower, upper) in zip(labels, nc2_groups[::-1])]
    ax.legend(handles=legend_handles, labels=new_labels[::-1], frameon=False, loc="upper right", title="genome completeness", fontsize=10, title_fontsize=11, bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

def plot_virus_codetection_heatmap_discrete(codetection_matrix: pl.DataFrame, output_file: str, categories: dict) -> None:
    category_values = {}
    for category in categories:
        category_values[category] = categories[category]["value"]
    category_colours = {}
    for category in categories:
        category_colours[category] = categories[category]["colour"]
    virus_cols = [col for col in codetection_matrix.columns if col != "sample"]
    plot_matrix = codetection_matrix.select(
        [pl.col(col).replace(category_values).cast(pl.Int64) for col in virus_cols]
    ).to_numpy()
    # try clustering the rows and columns of the codetection matrix using seaborn clustermap
    clustered_matrix = sns.clustermap(
        plot_matrix,
        method="average",
        metric="euclidean",
        cmap=ListedColormap([category_colours[label] for label in category_values]),
        norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], len(category_colours)),
        figsize=(10, 10),
        xticklabels=virus_cols,
        yticklabels=codetection_matrix.get_column("sample").to_list(),
    )
    row_order = clustered_matrix.dendrogram_row.reordered_ind
    col_order = clustered_matrix.dendrogram_col.reordered_ind
    plot_matrix = plot_matrix[row_order, :][:, col_order]
    sample_labels = list(np.array(codetection_matrix.get_column("sample").to_list())[row_order])
    virus_cols = list(np.array(virus_cols)[col_order])
    end_viruses = [virus for virus in ["SINDV", "LSV"] if virus in virus_cols]
    ordered_virus_cols = [virus for virus in virus_cols if virus not in end_viruses] + end_viruses
    virus_order = [virus_cols.index(virus) for virus in ordered_virus_cols]
    plot_matrix = plot_matrix[:, virus_order]
    virus_cols = ordered_virus_cols
    plt.close(clustered_matrix.fig)
    n_samples = codetection_matrix.height
    # figure_width = max(4.5, min(18, 1.3 + (0.12 * n_samples)))
    # figure_height = max(6, 1.2 + (0.55 * len(virus_cols)))
    figure_width = 11.5
    figure_height = 2.5
    sns.set_theme(context="paper", style="white", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    cmap = ListedColormap([category_colours[label] for label in category_values])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(plot_matrix.T, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_yticks(range(len(virus_cols)))
    ax.set_yticklabels(virus_cols, rotation=0, ha="right", rotation_mode="anchor")
    if n_samples <= 80:
        ax.set_xticks(range(n_samples))
        ax.set_xticklabels(sample_labels, rotation=90, ha="center", fontsize=6)
    else:
        ax.set_xticks([])
    ax.set_ylabel("virus")
    ax.set_xlabel("sample")
    ax.set_yticks(np.arange(-0.5, len(virus_cols), 1), minor=True)
    ax.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend_handles = [
        Patch(facecolor=colour, edgecolor="none", label=label)
        for label, colour in category_colours.items()
        if label != "neither"
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), title="detection")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)

def plot_virus_codetection_heatmap_continuous_with_boxplot(
    raw_data: pl.DataFrame,
    value_matrix: pl.DataFrame,
    output_file: str,
    cmap: str,
    outline_by: pl.DataFrame = None,
    outline_colours: dict = None,
    boxplot_colours: dict = None,
) -> None:
    virus_cols = [col for col in value_matrix.columns if col != "sample"]
    plot_matrix = value_matrix.select(virus_cols).to_numpy()
    outline_matrix = None
    if outline_by is not None and outline_colours is not None:
        outline_matrix = outline_by.select(virus_cols).to_numpy()
    # custom order
    virus_counts = np.nansum(plot_matrix > 0, axis=0)
    virus_order = np.argsort(-virus_counts, kind="stable")
    plot_matrix = plot_matrix[:, virus_order]
    if outline_matrix is not None:
        outline_matrix = outline_matrix[:, virus_order]
    virus_cols = [virus_cols[i] for i in virus_order]
    sample_labels = value_matrix.get_column("sample").to_list()
    detection_matrix = np.nan_to_num(plot_matrix, nan=0.0) > 0
    detection_counts = detection_matrix.sum(axis=1)
    no_detection_rank = len(virus_cols)
    primary_virus_ranks = np.array(
        [
            int(np.flatnonzero(sample_detections)[0])
            if sample_detections.any()
            else no_detection_rank
            for sample_detections in detection_matrix
        ]
    )
    primary_ct_values = np.array(
        [
            plot_matrix[sample_index, primary_rank]
            if primary_rank < no_detection_rank
            else np.inf
            for sample_index, primary_rank in enumerate(primary_virus_ranks)
        ]
    )
    primary_ct_values = np.nan_to_num(primary_ct_values, nan=np.inf)
    row_order = sorted(
        range(plot_matrix.shape[0]),
        key=lambda sample_index: (
            primary_virus_ranks[sample_index],
            -detection_counts[sample_index],
            primary_ct_values[sample_index],
            sample_labels[sample_index],
        ),
    )
    plot_matrix = plot_matrix[row_order, :]
    if outline_matrix is not None:
        outline_matrix = outline_matrix[row_order, :]
    sample_labels = [sample_labels[index] for index in row_order]
    # Retain samples and viruses detected by either PCR (positive Ct) or
    # COSMOS (a genome-coverage category that has a visible outline).
    ordered_detection_matrix = np.nan_to_num(plot_matrix, nan=0.0) > 0
    if outline_matrix is not None:
        outline_detection_matrix = np.vectorize(
            lambda category: (
                category in outline_colours
                and outline_colours[category] != "None"
            ),
            otypes=[bool],
        )(outline_matrix)
        ordered_detection_matrix |= outline_detection_matrix
    sample_detection_mask = ordered_detection_matrix.any(axis=1)
    virus_detection_mask = ordered_detection_matrix.any(axis=0)
    plot_matrix = plot_matrix[sample_detection_mask, :][:, virus_detection_mask]
    if outline_matrix is not None:
        outline_matrix = outline_matrix[sample_detection_mask, :][:, virus_detection_mask]
    sample_labels = [
        label for label, keep in zip(sample_labels, sample_detection_mask) if keep
    ]
    virus_cols = [
        virus for virus, keep in zip(virus_cols, virus_detection_mask) if keep
    ]
    n_samples = plot_matrix.shape[0]
    # figure_width = max(4.5, min(18, 1.3 + (0.12 * n_samples)))
    # figure_height = max(6, 1.2 + (0.55 * len(virus_cols)))
    figure_width = 10
    figure_height = 2.5
    sns.set_theme(context="paper", style="white", font_scale=1.0)
    fig, axs = plt.subplots(
        ncols=4,
        nrows=1,
        figsize=(figure_width, figure_height),
        width_ratios=[0.77, 0.02, 0.02, 0.19],
        gridspec_kw={"wspace": 0.04},
    )
    ax = axs[0]
    plot_matrix = plot_matrix.T
    if outline_matrix is not None:
        outline_matrix = outline_matrix.T
    cmap_obj = plt.get_cmap(cmap).reversed()
    positive_ct_values = plot_matrix[np.isfinite(plot_matrix) & (plot_matrix > 0)]
    ct_min = 5
    ct_max = 38
    ct_ticks = np.arange(10, 36, 5)
    norm = Normalize(vmin=ct_min, vmax=ct_max)
    outline_x_inset = 0.12
    outline_y_inset = 0.02
    new_detection_colour = "#80FB82"
    for i in range(plot_matrix.shape[0]):
        for j in range(plot_matrix.shape[1]):
            value = plot_matrix[i, j]
            edgecolour = (
                outline_colours.get(outline_matrix[i, j], "None")
                if outline_matrix is not None
                else "None"
            )
            has_visible_outline = edgecolour != "None"
            if not np.isfinite(value) or value == 0:
                facecolour = new_detection_colour if has_visible_outline else "white"
            else:
                facecolour = cmap_obj(norm(value))
            ax.add_patch(
                plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    facecolor=facecolour,
                    edgecolor="none",
                    linewidth=0,
                    antialiased=False,
                )
            )
            if has_visible_outline:
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5 + outline_x_inset, i - 0.5 + outline_y_inset),
                        1 - 2 * outline_x_inset,
                        1 - 2 * outline_y_inset,
                        fill=False,
                        edgecolor=edgecolour,
                        linewidth=0.8,
                        joinstyle="miter",
                        antialiased=True,
                        zorder=2,
                    )
                )
    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(len(virus_cols) - 0.5, -0.5)
    ax.set_aspect("auto")
    ax.set_yticks(range(len(virus_cols)))
    ax.set_yticklabels(virus_cols, rotation=0, ha="right", rotation_mode="anchor")
    if n_samples <= 80:
        ax.set_xticks(range(n_samples))
        ax.set_xticklabels(sample_labels, rotation=90, ha="center", fontsize=6)
    else:
        ax.set_xticks([])
    ax.set_ylabel("Virus")
    ax.set_title("Sample")
    ax.legend(
        handles=[
            # Patch(
            #     facecolor="white",
            #     edgecolor="black",
            #     linewidth=0.5,
            #     label=">5% genome coverage",
            # ),
            Patch(
                facecolor=new_detection_colour,
                edgecolor="black",
                linewidth=0,
                label="new detection",
            )
        ],
        loc="upper left",
        bbox_to_anchor=(0, 0),
        frameon=False,
        fontsize=7,
        handlelength=1.0,
        handleheight=0.8,
        borderaxespad=0,
        ncols=2
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#B8B8B8")
        spine.set_linewidth(0.01)
    # plot a cmap colourbar legend to the right of the plot
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj), cax=axs[1])
    cbar.outline.set_visible(False)
    cbar.set_ticks(ct_ticks)
    axs[1].set_title("Cycle\nthreshold", fontsize=10)
    axs[1].tick_params(axis="y", labelsize=8)
    for spine in axs[1].spines.values():
        spine.set_visible(True)
        spine.set_color("#B8B8B8")
        spine.set_linewidth(0.01)
    axs[2].set_axis_off()
    box_ax = axs[3]
    boxplot_order = list(boxplot_colours.keys())[::-1]
    hue_order = list(boxplot_colours.keys())[::-1]
    boxplot_data = filter_pcr_positive(raw_data)
    stripplot_point_count = boxplot_data.filter(
        pl.col("genome_coverage_category").is_in(boxplot_order)
    ).height
    heatmap_detection_count = positive_ct_values.size
    assert stripplot_point_count == heatmap_detection_count, (
        f"Stripplot contains {stripplot_point_count} points, but the heatmap "
        f"contains {heatmap_detection_count} positive-Ct cells"
    )
    sns.boxplot(
        data=boxplot_data.to_pandas(),
        x="genome_coverage_category",
        y="ct_value",
        hue="genome_coverage_category",
        order=boxplot_order,
        hue_order=hue_order,
        palette={category: "white" for category in boxplot_order},
        showfliers=False,
        width=0.62,
        linewidth=1.1,
        dodge=False,
        legend=False,
        ax=box_ax,
    )
    box_edge_colours = [boxplot_colours[category] for category in boxplot_order]
    for patch, edgecolour in zip(box_ax.patches, box_edge_colours):
        patch.set_facecolor("white")
        patch.set_edgecolor(edgecolour)
        patch.set_linewidth(1.1)
    for line in box_ax.lines:
        x_values = np.asarray(line.get_xdata(), dtype=float)
        if x_values.size and np.isfinite(x_values).any():
            category_index = int(round(np.nanmean(x_values)))
            if 0 <= category_index < len(box_edge_colours):
                line.set_color(box_edge_colours[category_index])
    box_ax.set_ylim(ct_min, ct_max)
    box_ax.set_xlabel("Genome completeness", fontsize=10)
    box_ax.set_ylabel("")
    box_ax.tick_params(
        axis="y",
        labelleft=False,
        left=True,
        length=3,
        width=0.6,
        color="#666666",
    )
    box_ax.grid(True, which="major", axis="y", color="#B8B8B8", linewidth=0.9, alpha=0.55)
    for spine in box_ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#B8B8B8")
        spine.set_linewidth(0.01)
    sns.stripplot(
        data=boxplot_data.to_pandas(),
        x="genome_coverage_category",
        y="ct_value",
        order=boxplot_order,
        hue_order=hue_order,
        ax=box_ax,
        color="#444444",
        dodge=False,
        jitter=0.18,
        size=3.6,
        alpha=0.6,
        linewidth=0,
    )
    box_ax.set_xticks(range(len(boxplot_order)))
    box_ax.set_xticklabels(
        boxplot_order,
        rotation=35,
        ha="right",
        rotation_mode="anchor",
        fontsize=8,
    )
    box_ax.tick_params(which="both", axis="x", bottom=True, length=3, width=0.6, pad=2)
    box_ax.set_yticks(ct_ticks)
    fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.95, top=0.85, bottom=0.3, hspace=0, wspace=0.04)
    fig.savefig(output_file.replace(".png", ".pdf"), dpi=600)
    plt.close(fig)

def prepare_virus_overview_scatter_data(
    combined_df: pl.DataFrame,
    nc2_groups: list,
    nc2_labels: list[str],
    nc2_colours: list[str],
) -> tuple[pl.DataFrame, list[dict], list[str], list[str], dict, dict]:
    plot_data = filter_pcr_positive(combined_df)
    completeness_specs = []
    for group, label, colour in zip(nc2_groups, nc2_labels, nc2_colours):
        completeness_specs.append(
            {
                "label": label,
                "range": f"{group[0]}-{group[1]}",
                "color": colour
            }
        )
    stack_order = nc2_labels[::-1]
    virus_order = sorted(plot_data.select(pl.col("virus")).unique().to_series().to_list())
    virus_markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "*", "h", "H", "8"]
    virus_marker_map = {virus: virus_markers[index % len(virus_markers)] for index, virus in enumerate(virus_order)}
    completeness_palette = {spec["label"]: spec["color"] for spec in completeness_specs}
    completeness_data = plot_data.with_columns(
        pl.col("nc2").fill_null(0).map_elements(
            lambda nc2: next(
                (
                    completeness_specs[index]["label"]
                    for index, (lower, upper) in enumerate(nc2_groups)
                    if lower <= nc2 < upper or (index == len(nc2_groups) - 1 and lower <= nc2 <= upper)
                ),
                "unknown",
            ),
            return_dtype=pl.Utf8,
        ).alias("completeness_category"),
        pl.when(pl.col("dedup_reads").cast(pl.Float64, strict=False).fill_null(0) > 0)
        .then(pl.col("dedup_reads").cast(pl.Float64, strict=False))
        .otherwise(pl.lit(0.0))
        .alias("dedup_reads_for_plot"),
    )
    return (
        completeness_data,
        completeness_specs,
        stack_order,
        virus_order,
        virus_marker_map,
        completeness_palette,
    )


def style_unique_read_count_axis(ax: plt.Axes) -> None:
    # Use the same transformation as the model, with genuine zero-count ticks.
    ax.set_yscale("function", functions=(lambda values: np.log10(np.maximum(values, -0.999999) + 1), inverse_log10_reads))
    ax.set_ylim(bottom=0)
    ax.set_yticks([0, 10, 100, 1_000, 10_000, 100_000, 1_000_000])
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.yaxis.set_minor_locator(mticker.NullLocator())


def plot_virus_overview_scatter(
    combined_df: pl.DataFrame,
    output_file: str,
    virus_acronyms: dict,
    nc2_groups: list,
    nc2_labels: list[str],
    nc2_colours: list[str],
) -> None:
    (
        completeness_data,
        completeness_specs,
        stack_order,
        virus_order,
        virus_marker_map,
        completeness_palette,
    ) = prepare_virus_overview_scatter_data(
        combined_df,
        nc2_groups,
        nc2_labels,
        nc2_colours,
    )
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    sns.scatterplot(
        data=completeness_data.to_pandas(),
        x="ct_value",
        y="dedup_reads_for_plot",
        hue="completeness_category",
        hue_order=stack_order,
        style="virus",
        style_order=virus_order,
        markers=virus_marker_map,
        palette=completeness_palette,
        s=34,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.25,
        legend=False,
        ax=ax,
    )
    style_unique_read_count_axis(ax)
    plot_line_of_best_fit_to_ax(
        completeness_data,
        x="ct_value",
        y="dedup_reads",
        ax=ax,
        logx=False,
        logy=True,
        stats_output_file=output_file.replace(".png", ".line_of_best_fit_stats.tsv"),
    )
    ax.set_xlabel("cycle threshold in original sample")
    ax.set_ylabel("unique read count")
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    completeness_legend = ax.legend(
        handles=[
            Patch(facecolor=spec["color"], edgecolor="none", label=f'{spec["label"]} ({spec["range"]})')
            for spec in completeness_specs
        ],
        title="genome completeness",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        borderaxespad=0.0,
        fontsize=10,
        title_fontsize=11,
    )
    ax.legend(
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
        title="virus",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.55),
        frameon=False,
        borderaxespad=0.0,
        fontsize=10,
        title_fontsize=11,
        ncol=2,
    )
    ax.add_artist(completeness_legend)
    fig.subplots_adjust(left=0.2, right=0.75, top=0.96, bottom=0.14)
    fig.savefig(output_file, dpi=300)#, bbox_inches="tight")
    plt.close(fig)


def plot_virus_results_overview(combined_df: pl.DataFrame, output_file: str, virus_acronyms: dict,
                                inset_data: pl.DataFrame,
                                nc2_groups: list = [(0, 0.2), (0.2, 0.95), (0.95, 1)], nc2_labels=["missing", "incomplete", "complete"], nc2_colours=["#D9D9D2", "#B7C7A0", "#4E6E5D"],
                                ) -> None:
    (
        completeness_data,
        completeness_specs,
        stack_order,
        virus_order,
        virus_marker_map,
        completeness_palette,
    ) = prepare_virus_overview_scatter_data(
        combined_df,
        nc2_groups,
        nc2_labels,
        nc2_colours,
    )
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig = plt.figure(figsize=(11.5, 9.2))
    grid = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[1.5, 1.1], width_ratios=[1, 2], hspace=0.3, wspace=0.4)

    scatter_ax = fig.add_subplot(grid[0, :])
    sns.scatterplot(
        data=completeness_data.to_pandas(),
        x="ct_value",
        y="dedup_reads_for_plot",
        hue="completeness_category",
        hue_order=stack_order,
        style="virus",
        style_order=virus_order,
        markers=virus_marker_map,
        palette=completeness_palette,
        s=34,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.25,
        legend=False,
        ax=scatter_ax,
    )
    style_unique_read_count_axis(scatter_ax)

    # plot line of best fit onto scatter
    plot_line_of_best_fit_to_ax(
        completeness_data,
        x="ct_value",
        y="dedup_reads",
        ax=scatter_ax,
        logx=False,
        logy=True,
        stats_output_file=output_file.replace(".png", ".line_of_best_fit_stats.tsv"),
    )

    scatter_ax.set_xlabel("cycle threshold in original sample")
    scatter_ax.set_ylabel("unique read count")
    scatter_ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    scatter_ax.spines["top"].set_visible(False)
    scatter_ax.spines["right"].set_visible(False)
    scatter_xlim = scatter_ax.get_xlim()
    scatter_ylim = scatter_ax.get_ylim()

    inset_ax = scatter_ax.inset_axes([0.08, 0.08, 0.34, 0.34])
    inset_ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    inset_ax.spines["top"].set_visible(False)
    inset_ax.spines["right"].set_visible(False)
    sns.scatterplot(
        data=inset_data.to_pandas(),
        x="ct",
        y="n_reads_dedup",
        ax=inset_ax,
        s=30,
        alpha=0.8,
        color='black',
        edgecolor="white",
        linewidth=0.25,
        legend=False,
    )
    style_unique_read_count_axis(inset_ax)
    plot_line_of_best_fit_to_ax(
        inset_data,
        x="ct",
        y="n_reads_dedup",
        ax=inset_ax,
        logx=False,
        logy=True,
        stats_output_file=output_file.replace(".png", ".inset_line_of_best_fit_stats.tsv"),
    )
    inset_ax.set_xlabel("")
    inset_ax.set_ylabel("")
    inset_ax.set_xticks([20, 25, 30])
    inset_ax.tick_params(axis="both", which="major", labelsize=8, pad=-4)
    # inset_ax.set_xlim(scatter_xlim)
    # inset_ax.set_ylim(scatter_ylim)

    completeness_legend_handles = [
        Patch(facecolor=spec["color"], edgecolor="none", label=f'{spec["label"]} ({spec["range"]})')
        for spec in completeness_specs
    ]
    completeness_legend = scatter_ax.legend(
        handles=completeness_legend_handles,
        title="genome completeness",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        borderaxespad=0.0,
        fontsize=12,
        title_fontsize=14,
    )

    virus_legend_handles = [
        Line2D([0], [0], marker=virus_marker_map[virus], linestyle="None", color="black", markerfacecolor="black", markersize=7, label=virus_acronyms.get(virus, virus))
        for virus in virus_order
    ]
    virus_legend = scatter_ax.legend(
        handles=virus_legend_handles,
        title="virus",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.58),
        frameon=False,
        borderaxespad=0.0,
        fontsize=12,
        title_fontsize=14,
        ncol=2,
    )
    scatter_ax.add_artist(completeness_legend)

    box_ax = fig.add_subplot(grid[1, 0])
    sns.boxplot(
        data=completeness_data.to_pandas(),
        x="completeness_category",
        y="ct_value",
        hue="completeness_category",
        order=stack_order,
        hue_order=stack_order,
        palette={spec["label"]: spec["color"] for spec in completeness_specs},
        ax=box_ax,
        showfliers=False,
        width=0.62,
        linewidth=1.1,
        dodge=False,
        legend=False,
    )
    sns.stripplot(
        data=completeness_data.to_pandas(),
        x="completeness_category",
        y="ct_value",
        order=stack_order,
        ax=box_ax,
        color="#444444",
        dodge=False,
        jitter=0.18,
        size=3.6,
        alpha=0.6,
        linewidth=0,
    )
    box_ax.set_xlabel("virus genome completeness")
    box_ax.set_ylabel("virus cycle threshold in\noriginal sample")
    box_ax.grid(True, which="major", axis="y", linewidth=0.8, alpha=0.22)
    box_ax.spines["top"].set_visible(False)
    box_ax.spines["right"].set_visible(False)
    box_ax.set_xticks(range(len(stack_order)))
    box_ax.set_xticklabels(stack_order, rotation=35, ha="right", rotation_mode="anchor")
    box_ax.tick_params(which="both", axis="x", pad=-5)

    bar_data = (filter_dedup_read_positive(combined_df)
    .filter(pl.col("nc2").fill_null(0) >= 0.05) # take only COSMOS detections
    .with_columns(
        pl.col("nc2").fill_null(0).map_elements(
            lambda nc2: next(
                (
                    completeness_specs[index]["label"]
                    for index, (lower, upper) in enumerate(nc2_groups)
                    if lower <= nc2 < upper or (index == len(nc2_groups) - 1 and lower <= nc2 <= upper)
                ),
                "unknown",
            ),
            return_dtype=pl.Utf8,
        ).alias("completeness_category")
    ))
    bar_counts = bar_data.group_by(["virus", "completeness_category"]).agg(pl.len().alias("count"))
    bar_pivot = bar_counts.pivot(index="virus", columns="completeness_category", values="count").fill_null(0)
    bar_pivot_pd = bar_pivot.to_pandas().set_index("virus").reindex(columns=stack_order, fill_value=0)
    bar_pivot_pd = bar_pivot_pd.sort_values(by=stack_order[0], ascending=False)
    # put SINDV and LSV at end
    virus_order = [virus for virus in bar_pivot_pd.index.to_list() if virus_acronyms.get(virus) not in ["SINDV", "LSV"]] + [virus for virus in bar_pivot_pd.index.to_list() if virus_acronyms.get(virus) in ["SINDV", "LSV"]]
    bar_pivot_pd = bar_pivot_pd.reindex(virus_order)

    bar_ax = fig.add_subplot(grid[1, 1])
    bar_pivot_pd.plot(kind="bar", stacked=True, color=[spec["color"] for spec in completeness_specs[::-1]], ax=bar_ax, legend=False, width=0.82)
    bar_ax.set_xlabel("virus")
    bar_ax.set_ylabel("number of samples")
    bar_ax.grid(False)
    bar_ax.spines["top"].set_visible(False)
    bar_ax.spines["right"].set_visible(False)
    bar_ax.set_xticks(range(len(bar_pivot_pd.index)))
    bar_ax.set_xticklabels([virus_acronyms.get(virus, virus) for virus in bar_pivot_pd.index], rotation=35, ha="right", rotation_mode="anchor")
    bar_ax.tick_params(which="both", axis="x", pad=-5)

    # let's add a little asterisks to the bars where SINDV and LSV to highlight them as additional discoveries
    for index, virus in enumerate(bar_pivot_pd.index):
        if virus_acronyms.get(virus) in ["SINDV", "LSV"]:
            total_count = bar_pivot_pd.loc[virus].sum()
            bar_ax.text(index, total_count + 0.1, "*", ha="center", va="bottom", fontsize=30)

    fig.subplots_adjust(left=0.18, right=0.74, top=0.97, bottom=0.18)
    plt.savefig(output_file)

def plot_virus_results(
    combined_df: pl.DataFrame,
    output_dir: str,
    virus_acronyms: dict,
    inset_results_file: str,
    genome_coverage_groups: list,
    genome_coverage_labels: list[str],
    genome_coverage_colours: list[str],
) -> None:
    plot_dedup_reads_vs_ct_value(combined_df, output_dir + "/dedup_reads_vs_ct_value.png", virus_acronyms)
    plot_virus_overview_scatter(
        combined_df,
        output_dir + "/virus_ct_unique_reads_scatter.png",
        virus_acronyms,
        genome_coverage_groups,
        genome_coverage_labels,
        genome_coverage_colours,
    )
    plot_box_plot_of_dedup_reads_by_ct_value_coloured_by_nc2(
        filter_pcr_positive(combined_df),
        output_dir + "/box_plot_dedup_reads_by_ct_value_coloured_by_nc2.png",
        genome_coverage_labels,
        genome_coverage_colours,
    )
    plot_stacked_bar_plot_for_each_virus_with_missed_incomplete_and_complete_categories(
        filter_dedup_read_positive(combined_df),
        output_dir + "/stacked_bar_plot_by_virus.png",
        genome_coverage_groups,
        virus_acronyms,
        genome_coverage_labels,
        genome_coverage_colours,
    )
    inset_results = pl.read_csv(inset_results_file)
    plot_virus_results_overview(combined_df, 
                                output_dir + "/virus_results_overview.png", 
                                virus_acronyms,
                                inset_data=inset_results,
                                nc2_groups=genome_coverage_groups,
                                nc2_labels=genome_coverage_labels,
                                nc2_colours=genome_coverage_colours)

def draw_colourbar_by_metadata( # from abundance-geospatial module
    ax: plt.Axes,
    metadata: pl.DataFrame,
    col: str,
    cmap=None,
    reference_ax: plt.Axes | None = None,
    label_map: dict | None = None,
    palette: dict | None = None,
) -> dict:
    AXIS_LABEL_SIZE = 11
    TICK_LABEL_SIZE = 8
    ANNOTATION_LABEL_SIZE = 9
    LEGEND_LABEL_SIZE = 10
    values = metadata[col].to_list()
    if palette is not None:
        value_to_color = {value: to_rgba(color) for value, color in palette.items()}
    else:
        unique_values = metadata.select(pl.col(col).drop_nulls().unique()).to_series().to_list()
        denom = max(len(unique_values), 1)
        value_to_color = {value: cmap(index / denom) for index, value in enumerate(unique_values)}
    value_to_color[None] = to_rgba("#FFFFFF")
    colors = [value_to_color.get(value, (1.0, 1.0, 1.0, 1.0)) for value in values]
    color_array = np.array(colors).reshape(1, len(colors), 4)
    ax.imshow(color_array, aspect="auto", interpolation="nearest", extent=[0, len(colors), 0, 1])

    longest_run_by_value = {}
    run_value = None
    run_start = 0
    for index, value in enumerate(values + [None]):
        if value != run_value:
            if run_value is not None:
                run_len = index - run_start
                previous = longest_run_by_value.get(run_value)
                if previous is None or run_len > previous[1]:
                    longest_run_by_value[run_value] = (run_start, run_len)
            run_value = value
            run_start = index

    def text_color(rgba):
        r, g, b = rgba[:3]
        luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
        return "black" if luminance > 0.6 else "white"

    min_segment_width = max(3, int(0.03 * len(values)))
    for value, (start, run_len) in longest_run_by_value.items():
        if run_len < min_segment_width:
            continue
        label = label_map.get(value, value) if label_map is not None else value
        ax.text(
            start + run_len / 2,
            0.5,
            str(label),
            ha="center",
            va="center",
                fontsize=ANNOTATION_LABEL_SIZE,
            color=text_color(value_to_color.get(value, (1.0, 1.0, 1.0, 1.0))),
            clip_on=True,
        )

    if reference_ax is not None:
        ax.set_xlim(reference_ax.get_xlim())
    ax.set_yticks([])
    ax.set_xticks([])
    sns.despine(ax=ax, left=True, bottom=True, top=True, right=True)
    return value_to_color


def plot_virus_codetection_with_boxplot_with_cbars(
    data: pl.DataFrame,
    metadata: pl.DataFrame,
    plot_feature: str,
    cbar_features: list[str],
    boxplot_x_feature: str,
    boxplot_y_feature: str,
    output_path: str,
    boxplot_order: list[str],
    virus_acronyms: dict,
    cmap="flare",
    outline_by=None,
) -> None:
    # where virus is sindbis/SINV, change nc2 = 0 to NaN
    # data = data.with_columns(
    #     pl.when((pl.col("virus") == "sindbisvirus") & (pl.col("nc2") <  0.01)).then(pl.lit(np.nan)).otherwise(pl.col("nc2")).alias("nc2")
    # )

    # going to apply a blanket filter of nc2 > 0.05 for all viruses
    data = data.with_columns(
        pl.when(pl.col("nc2") < 0.05).then(pl.lit(np.nan)).otherwise(pl.col("nc2")).alias("nc2")
    )

    # palettes from geospatial module
    COASTAL_PALETTE = {
        1: "#8EBFE8",
        0: "#8FD19E",
        None: "#D9D9D9",
    }
    COASTAL_LABEL_MAP = {
        1: "coastal",
        0: "inland",
        None: "unknown",
    }

    URBAN_RURAL_PALETTE = {
        "Urban": "#EFBA74",
        "Met_Fringe": "#E38975",
        "Rural": "#7A1F5C",
    }
    URBAN_RURAL_LABEL_MAP = {
        "Urban": "urban",
        "Met_Fringe": "met fringe",
        "Rural": "rural",
    }

    # now let's order cols by cbar features in metadata
    custom_sort_order = {
        "coastal": [0, 1],
        "urban_rural": ["Rural", "Met_Fringe", "Urban"],
    }
    for feature, order in custom_sort_order.items():
        if feature in cbar_features:
            order_map = {value: index for index, value in enumerate(order)}
            metadata = metadata.with_columns(
                pl.col(feature)
                .replace(order_map)
                .cast(pl.Int64)
                .alias(f"{feature}_sort")
            )
    sample_order = metadata.sort(by=[f"{feature}_sort" for feature in cbar_features]).get_column("sample").to_list()
    sorted_metadata = metadata.filter(pl.col("sample").is_in(sample_order)).sort(by=[f"{feature}_sort" for feature in cbar_features])
    # generate matrix of plot feature for heatmap
    plot_matrix_df = generate_matrix(data, virus_acronyms, plot_feature, nan_value=np.nan).filter(pl.col("sample").is_in(sample_order))
    plot_matrix_df = plot_matrix_df.sort("sample").with_columns(
        pl.col("sample").replace({sample: index for index, sample in enumerate(sample_order)}).cast(pl.Int64).alias("sample_order")
    ).sort("sample_order")
    virus_order = plot_matrix_df.drop(["sample", "sample_order"]).columns
    plot_matrix = plot_matrix_df.drop(["sample", "sample_order"]).to_numpy().T
    # order rows by virus detection counts
    detection_counts = np.nansum(plot_matrix > 0, axis=1)
    row_order = np.argsort(-detection_counts, kind="stable")
    plot_matrix = plot_matrix[row_order, :]
    virus_order = [virus_order[index] for index in row_order]
    if outline_by is not None:
        outline_matrix_df = generate_matrix(data, virus_acronyms, outline_by, nan_value=np.nan).filter(pl.col("sample").is_in(sample_order))
        outline_matrix_df = outline_matrix_df.sort("sample").with_columns(
            pl.col("sample").replace({sample: index for index, sample in enumerate(sample_order)}).cast(pl.Int64).alias("sample_order")
        ).sort("sample_order")
        outline_matrix = outline_matrix_df.drop(["sample", "sample_order"]).to_numpy().T # implement later if want to use
        outline_matrix = outline_matrix[row_order, :]
    # remove KUN and put SINDV and LSV at the end
    row_order = [index for index, virus in enumerate(virus_order) if virus not in ["KUN", "SINDV", "LSV"]] + [index for index, virus in enumerate(virus_order) if virus in ["SINDV", "LSV"]]
    plot_matrix = plot_matrix[row_order, :]
    if outline_by is not None:
        outline_matrix = outline_matrix[row_order, :]
    virus_order = [virus_order[index] for index in row_order]

    # let's build figure - heatmap with cbars below and boxplot to the right. Boxplot covers all rows
    fig = plt.figure(figsize=(10, 3 + 0.5 * len(cbar_features)))
    grid = fig.add_gridspec(
        nrows=4,
        ncols=2,
        height_ratios=[0.08, 1, 0.15, 0.15],
        width_ratios=[0.85, 0.15],
        hspace=0.35,
        wspace=0.15,
    )

    heatmap_ax = fig.add_subplot(grid[0:2, 0])
    heatmap_cbar_ax = fig.add_subplot(grid[0, 1])
    boxplot_ax = fig.add_subplot(grid[1:4, 1])
    cbar_axes = [
        fig.add_subplot(grid[2, 0]),
        fig.add_subplot(grid[3, 0]),
    ]

    # plot heatmap
    print(plot_matrix.shape)
    masked_matrix = np.ma.masked_invalid(plot_matrix)
    heatmap_cmap = plt.get_cmap(cmap)
    heatmap_cmap.set_bad(color="#FCF6E9FF") # set NaN values to light grey
    heatmap_ax.imshow(masked_matrix, aspect="auto", interpolation="nearest", cmap=heatmap_cmap, vmin=0, vmax=1)
    heatmap_ax.set_yticks(range(len(plot_matrix)))
    ytick_labels = [i + "*" if i in ["SINDV", "LSV"] else i for i in virus_order]
    heatmap_ax.set_yticklabels(ytick_labels, rotation=0, ha="right", rotation_mode="anchor")
    heatmap_ax.set_xticks([])
    heatmap_ax.set_ylabel("Virus")
    # remove spines
    for spine in heatmap_ax.spines.values():
        spine.set_visible(False)

    # draw colorbar for heatmap
    cbar_norm = plt.Normalize(vmin=0, vmax=1)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=cbar_norm, cmap=heatmap_cmap), cax=heatmap_cbar_ax, orientation="horizontal")
    heatmap_cbar_ax.set_title("Proportion of genome\nrecovered", fontsize=8)
    # remove spines
    for spine in heatmap_cbar_ax.spines.values():
        spine.set_visible(False)
    # remove ticks
    heatmap_cbar_ax.tick_params(axis="x", which="both", bottom=False, top=True, labelbottom=False, labeltop=True)
    heatmap_cbar_ax.set_xticks([0.2, 0.4, 0.6, 0.8])
    heatmap_cbar_ax.set_xticklabels(["", "", "", ""], fontsize=8)
    # annotate either end with min and max values
    heatmap_cbar_ax.text(-0.1, 0.5, "0", transform=heatmap_cbar_ax.transAxes, ha="left", va="center", fontsize=8)
    heatmap_cbar_ax.text(1.1, 0.5, "1", transform=heatmap_cbar_ax.transAxes, ha="right", va="center", fontsize=8)
    # annotate a tiny triangle at 0.3 below
    heatmap_cbar_ax.scatter(0.3, -0.35, marker="^", s=20, color="black", transform=heatmap_cbar_ax.transAxes, clip_on=False)

    # plot cbars for each feature
    for i, feature in enumerate(cbar_features):
        if feature == "coastal":
            palette = COASTAL_PALETTE
            label_map = COASTAL_LABEL_MAP
        elif feature == "urban_rural":
            palette = URBAN_RURAL_PALETTE
            label_map = URBAN_RURAL_LABEL_MAP
        else:
            raise ValueError(f"Unknown cbar feature: {feature}")
        draw_colourbar_by_metadata(
            ax=cbar_axes[i],
            metadata=sorted_metadata,
            col=feature,
            label_map=label_map,
            palette=palette
        )
    # # plot cbars for each feature
    # for i, feature in enumerate(cbar_features):
    #     cbar_ax = cbar_axes[i]
    #     feature_values = metadata.filter(pl.col("sample").is_in(sample_order)).sort("sample").with_columns(
    #         pl.col("sample").replace({sample: index for index, sample in enumerate(sample_order)}).cast(pl.Int64).alias("sample_order")
    #     ).sort("sample_order").select(feature).to_numpy().flatten()
    #     unique_values = np.unique(feature_values)
    #     if feature == "coastal":
    #         palette = COASTAL_PALETTE
    #     elif feature == "urban_rural":
    #         palette = URBAN_RURAL_PALETTE
    #     else:
    #         palette = {val: plt.cm.tab20(i) for i, val in enumerate(unique_values)}
    #     cbar_matrix = np.array([[to_rgba(palette[val]) for val in feature_values]])
    #     cbar_ax.imshow(cbar_matrix, aspect="auto")
    #     cbar_ax.set_yticks([])
    #     cbar_ax.set_xticks([])
    #     # remove spines
    #     for spine in cbar_ax.spines.values():
    #         spine.set_visible(False)
        
    # plot boxplot
    boxplot_data = filter_pcr_positive(data).with_columns(
        pl.col("nc2").fill_null(0)
    ).to_pandas()
    print(boxplot_data.shape)
    # print value counts for each category
    print(boxplot_data.groupby(boxplot_x_feature)[boxplot_y_feature].describe())
    sns.boxplot(
        data=boxplot_data,
        x=boxplot_x_feature,
        y=boxplot_y_feature,
        order=boxplot_order,
        showfliers=False,
        color="white",
        width=0.52,
        linewidth=1.1,
        dodge=False,
        legend=False,
        ax=boxplot_ax,
    )
    sns.stripplot(
        data=boxplot_data,
        x=boxplot_x_feature,
        y=boxplot_y_feature,
        order=boxplot_order,
        hue="nc2",
        ax=boxplot_ax,
        palette=cmap,
        # color="#444444",
        dodge=False,
        jitter=0.18,
        size=5,
        alpha=0.6,
        linewidth=0,
        legend=False,
    )
    boxplot_ax.set_xlabel("Genome recovery")
    boxplot_ax.set_ylabel("Cycle threshold")
    boxplot_ax.grid(True, which="major", axis="y", linewidth=0.8, alpha=0.22)
    # remove spines
    for spine in boxplot_ax.spines.values():
        spine.set_visible(False)
    xticks = boxplot_ax.get_xticks()
    xtick_labels = boxplot_ax.get_xticklabels()
    boxplot_ax.set_xticks(xticks)
    boxplot_ax.set_xticklabels(xtick_labels, fontsize=8)

    # save fig
    fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.95, top=0.85, bottom=0.2)
    fig.savefig(output_path, dpi=300)

def virus_coverage_plot_combined(df: pl.DataFrame, output_file: str) -> None:
    # to make relatve genome position, divide position by max position for each virus
    df = df.with_columns(
        max_position=pl.col("position").max().over("virus_acronym"),
    ).with_columns(
        relative_position=pl.col("position") / pl.col("max_position"),
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    # manually draw using mean, lower and upper shadding as lineplot very slow
    sns.lineplot(
        data=df.to_pandas(),
        x="relative_position",
        y="median_depth",
        hue="virus_acronym",
        palette="tab10",
        ax=ax
    )
    # add shading for lower and upper depth
    for virus in df.select("virus_acronym").unique().to_series().to_list():
        virus_df = df.filter(pl.col("virus_acronym") == virus).to_pandas()
        ax.fill_between(
            virus_df["relative_position"],
            virus_df["depth_q25"],
            virus_df["depth_q75"],
            alpha=0.2
        )
    ax.set_xlabel("Relative genome position")
    ax.set_ylabel("Median depth of coverage")
    style_unique_read_count_axis(ax)
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

def style_virus_coverage_yaxis(ax: plt.Axes) -> None:
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.4)
    ylims = ax.get_ylim()
    # select ticks that are within the limits
    ytick_values = [1,  100, 10_000]
    # ytick_values = [value for value in ytick_values if ylims[0] <= value <= ylims[1]]
    ax.set_yticks(ytick_values)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(
            lambda value, _: "0"
            if np.isclose(value, 0.5)
            else f"{int(value):,}"
        )
    )
    ax.yaxis.set_minor_locator(mticker.NullLocator())

def virus_coverage_plot_paneled(df: pl.DataFrame, output_file: str, fig: plt.Figure = None, subplot_spec: gridspec.GridSpec = None) -> None:
    standalone = fig is None
    # to make relatve genome position, divide position by max position for each virus
    # df = df.with_columns(
    #     max_position=pl.col("position").max().over("virus_acronym"),
    # ).with_columns(
    #     relative_position=pl.col("position") / pl.col("max_position"),
    # )
    virus_list = df.select("virus_acronym").unique().to_series().to_list()
    # order by sample size 
    sample_sizes = {}
    for virus in virus_list:
        sample_sizes[virus] = df.filter(pl.col("virus_acronym") == virus).select("n_samples").to_series().to_list()[0]
    virus_list = sorted(virus_list, key=lambda x: sample_sizes[x], reverse=True)
    n_viruses = len(virus_list)
    grid_kwargs = {
        "height_ratios": [1] * n_viruses + [0.25],
        "hspace": 0.4,
        "width_ratios": [0.025, 0.975],
        "wspace": 0.35,
    }
    if standalone:
        fig, axes = plt.subplots(
            nrows=n_viruses + 1,
            ncols=2,
            figsize=(6, 0.4 * n_viruses + 0.6),
            sharex="col",
            gridspec_kw=grid_kwargs,
        )
    else:
        if subplot_spec is None:
            raise ValueError("subplot_spec is required when drawing into an existing figure")
        inner_grid = subplot_spec.subgridspec(
            nrows=n_viruses + 1,
            ncols=2,
            **grid_kwargs,
        )
        axes = np.empty((n_viruses + 1, 2), dtype=object)
        shared_coverage_ax = None
        for row in range(n_viruses + 1):
            axes[row, 0] = fig.add_subplot(inner_grid[row, 0])
            axes[row, 1] = fig.add_subplot(
                inner_grid[row, 1],
                sharex=shared_coverage_ax,
            )
            if shared_coverage_ax is None:
                shared_coverage_ax = axes[row, 1]
    virus_axes = axes[:-1, 1]
    virus_label_axes = axes[:-1, 0]
    genome_position_ax = axes[-1, 1]
    unused_bottom_left_ax = axes[-1, 0]
    for i, virus in enumerate(virus_list):
        virus_df = df.filter(pl.col("virus_acronym") == virus).to_pandas()
        ax = virus_axes[i]
        sns.lineplot(
            data=virus_df,
            x="position",
            y="median_depth",
            color="black",
            ax=ax,
            linewidth=1.2
        )
        ax.fill_between(
            virus_df["position"],
            virus_df["depth_q25"],
            virus_df["depth_q75"],
            alpha=0.2,
            color="dimgray"
        )
        ax.set_ylabel(f"")
        style_virus_coverage_yaxis(ax)
        ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # for each virus, make the bottom axis stop at the max position for that virus without changing the xlim of the other axes
        bottom_spine = ax.spines["bottom"]
        bottom_spine.set_visible(False)
        # offset the y axis and draw the genome-length baseline separately
        left_spine = ax.spines["left"]
        left_spine.set_position(("outward", 5))
        left_spine.set_clip_on(False)
        left_spine.set_zorder(10)
        left_spine.set_linewidth(0.8)
        left_spine.set_color("black")
        genome_line = Line2D(
            [0, virus_df["position"].max()],
            [0, 0],
            transform=ax.get_xaxis_transform(),
            color="black",
            linewidth=0.8,
            zorder=10,
            clip_on=False,
        )
        fig.add_artist(genome_line)
        # retain the genome-length spine, but put ticks and labels on the separate scale axis
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        # y axis fontsize to 6
        ax.tick_params(axis="y", which="major", labelsize=6, pad=1)
        # add square bracket with virus label to side
        bracket_ax = virus_label_axes[i]
        bracket_ax.spines["left"].set_visible(True)
        bracket_ax.spines["right"].set_visible(False)
        # label the opposite side of the bracket with the virus name
        sample_size = sample_sizes[virus]
        label = f"*{virus}" if virus in ["SINDV", "LSV"] else virus
        bracket_ax.text(
            -0.35,
            0.62,
            label,
            ha="right",
            va="center",
            fontsize=9,
            transform=bracket_ax.transAxes,
        )
        bracket_ax.text(
            -0.3,
            0.30,
            f"(n={sample_size})",
            ha="right",
            va="center",
            fontsize=7,
            fontstyle="italic",
            transform=bracket_ax.transAxes,
        )
        bracket_ax.set_xticks([])
        bracket_ax.set_yticks([])

        ax.set_xlabel("")

    max_genome_position = df.select(
        pl.col("position").max()
    ).item()

    last_genome_tick = 12_000
    genome_tick_interval = 2_000
    # Show data through its actual final position
    genome_position_ax.set_xlim(0, max_genome_position)
    # Only draw ticks through 12,000
    genome_position_ax.set_xticks(
        np.arange(
            0,
            last_genome_tick + 1,
            genome_tick_interval,
        )
    )

    genome_position_ax.set_ylim(0, 1)
    genome_position_ax.set_xlabel("Genome position", fontsize=9)
    genome_position_ax.set_yticks([])
    genome_position_ax.grid(False)
    genome_position_ax.spines["top"].set_visible(False)
    genome_position_ax.spines["left"].set_visible(False)
    genome_position_ax.spines["right"].set_visible(False)
    genome_position_ax.spines["bottom"].set_bounds(0, last_genome_tick)
    genome_position_ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True, labelsize=8)
    unused_bottom_left_ax.set_axis_off()
    final_virus_position = virus_axes[-1].get_position()
    genome_axis_position = genome_position_ax.get_position()
    gap = 0.005

    genome_position_ax.set_position([
        final_virus_position.x0,
        final_virus_position.y0 - gap - genome_axis_position.height,
        final_virus_position.width,
        genome_axis_position.height,
    ])
    if standalone:
        fig.subplots_adjust(left=0.25, right=0.9, top=0.95, bottom=0.12, hspace=0.4, wspace=0.05)
        # add a ylabel for the whole figure
        fig.text(
            0.15,
            0.55,
            "Median read\ndepth",
            va="center",
            ha="right",
            rotation=0,
            fontsize=9,
            multialignment="right",
        )
        fig.savefig(output_file, dpi=300)
    else:
        return fig, axes
