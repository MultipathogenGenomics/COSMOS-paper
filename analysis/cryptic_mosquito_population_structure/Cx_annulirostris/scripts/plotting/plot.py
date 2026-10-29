# Plot cryptic subpopulation maps.

from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Patch, Wedge, Circle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import geopandas as gpd
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from analysis.common.regression import fit_linear_regression
from scipy.cluster import hierarchy
from Bio import Phylo
from ete3 import Tree
from io import StringIO

SYDNEY_COAST_LON_LAT = (151.28, -33.89)

def plot_sequence_dendrogram(tree_data: dict, output_file: str) -> None:
    metadata = tree_data["sequence_metadata"].set_index("sequence_id")
    labels = [metadata.loc[sequence_id, "leaf_label"] for sequence_id in tree_data["sequence_ids"]]
    figure_height = max(6, len(labels) * 0.1)
    fig, ax = plt.subplots(figsize=(9, figure_height))
    dendrogram_data = hierarchy.dendrogram(
        tree_data["linkage_matrix"],
        labels=labels,
        orientation="left",
        leaf_font_size=6,
        above_threshold_color="black",
        link_color_func=lambda _: "black",
        ax=ax,
    )
    group_by_sequence = metadata["group"].to_dict()
    leaf_groups = [
        group_by_sequence[tree_data["sequence_ids"][leaf_index]]
        for leaf_index in dendrogram_data["leaves"]
    ]
    group_positions = {}
    for leaf_index, group in enumerate(leaf_groups):
        group_positions.setdefault(group, []).append(5 + leaf_index * 10)
    group_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for group, positions in group_positions.items():
        ax.plot(
            [1.32, 1.32],
            [min(positions) - 3, max(positions) + 3],
            color="black",
            linewidth=1,
            transform=group_transform,
            clip_on=False,
        )
        ax.text(
            1.35,
            (min(positions) + max(positions)) / 2,
            group,
            fontsize=10,
            fontweight="bold",
            va="center",
            ha="left",
            transform=group_transform,
            clip_on=False,
        )
    if tree_data.get("pid_threshold") is not None:
        threshold = 1 - tree_data["pid_threshold"]
        ax.axvline(
            threshold,
            color="black",
            linestyle="--",
            linewidth=0.8,
            label=f"{tree_data['pid_threshold']:.2f} identity threshold",
        )
        ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.set_xlabel("proportion sequence difference")
    ax.set_ylabel("")
    ax.grid(False)
    sns.despine(ax=ax, top=True, right=True, left=True, bottom=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)

def draw_pie(ax: plt.Axes, x: float, y: float, ratios: list[float], radius: float, colours: list[str]) -> None:
    theta = 90
    for ratio, colour in zip(ratios, colours):
        next_theta = theta + 360 * ratio
        wedge = Wedge((x, y), radius, theta, next_theta, facecolor=colour, edgecolor="none")
        ax.add_patch(wedge)
        theta = next_theta

def draw_gradient_point(ax: plt.Axes, x: float, y: float, radius: float, colour: str) -> None:
    circle = Circle((x, y), radius, facecolor=colour, edgecolor="none")
    ax.add_patch(circle)

def draw_points_as_pies(ax: plt.Axes, points, ratio_columns: list[str], colours: list[str], radius: float) -> None:
    for _, row in points.iterrows():
        ratios = [row[column] for column in ratio_columns]
        draw_pie(ax, row.geometry.x, row.geometry.y, ratios, radius, colours)

def draw_points_as_gradients(ax: plt.Axes, points, ratio_columns: list[str], colours: list[str], radius: float) -> None:
    custom_cmap = sns.blend_palette(colours, as_cmap=True)
    for _, row in points.iterrows():
        ratios = [row[column] for column in ratio_columns]
        if len(ratios) > 2:
            raise ValueError("ERROR: gradient can't be used for more than two ratios.")
        # get a colour between two in list at that ratio
        colour = custom_cmap(ratios[0])
        draw_gradient_point(ax, row.geometry.x, row.geometry.y, radius, colour)

def get_sydney_inset_limits(points, crs) -> tuple[float, float, float, float]:
    penrith_row = points[points["location2"] == "Penrith"].iloc[0]
    penrith_x = penrith_row.geometry.x
    penrith_y = penrith_row.geometry.y
    coast_point = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy([SYDNEY_COAST_LON_LAT[0]], [SYDNEY_COAST_LON_LAT[1]]),
        crs="EPSG:4326",
    ).to_crs(crs)
    coast_x = coast_point.geometry.iloc[0].x
    coast_y = coast_point.geometry.iloc[0].y

    padding_left = 0.15
    padding_right = 0.15
    padding_top = 0.25
    padding_bottom = 0.15
    xmin = penrith_x - padding_left
    xmax = coast_x + padding_right
    ymin = min(penrith_y, coast_y) - padding_bottom
    ymax = max(penrith_y, coast_y) + padding_top
    return xmin, xmax, ymin, ymax

def add_panel_legend(ax: plt.Axes, title: str, labels: list[str], colours: list[str], to_add_to_legend: dict[str, str] | None = None) -> None:
    handles = [Patch(facecolor=colour, edgecolor="none", label=label) for label, colour in zip(labels, colours)]
    if to_add_to_legend:
        for label, colour in to_add_to_legend.items():
            handles.append(Patch(facecolor=colour, edgecolor="none", label=label))
    legend = ax.legend(handles=handles, title=title, loc="center right", bbox_to_anchor=(-0.03, 0.4), frameon=False, fontsize=9, title_fontsize=9)
    legend.get_title().set_fontstyle("italic")
    legend.get_title().set_ha("left")
    legend._legend_box.align = "left"

def plot_subpopulation_map_panel(fig, ax: plt.Axes, map_boundary, panel: dict, inset_bounds: list[float], to_add_to_legend: dict[str, str] | None = None, draw_type = "pie", map_xlim: tuple[float, float] | None = None) -> None:
    points = panel["points"]
    colours = panel["colours"]
    gradient_colours = panel.get("gradient_colours", colours)
    ratio_columns = panel["ratio_columns"]

    map_boundary.boundary.plot(ax=ax, color="black", linewidth=0.5)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if draw_type == "gradient":
        draw_points_as_gradients(ax, points, ratio_columns, gradient_colours, radius=0.2)
    else:
        draw_points_as_pies(ax, points, ratio_columns, colours, radius=0.2)
    if map_xlim is None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(map_xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_axis_off()

    ax_inset = fig.add_axes(inset_bounds)
    ax_inset.set_box_aspect(1)
    map_boundary.boundary.plot(ax=ax_inset, color="black", linewidth=0.5)
    if draw_type == "gradient":
        draw_points_as_gradients(ax_inset, points, ratio_columns, gradient_colours, radius=0.04)
    else:
        draw_points_as_pies(ax_inset, points, ratio_columns, colours, radius=0.04)
    xmin, xmax, ymin, ymax = get_sydney_inset_limits(points, points.crs)
    ax_inset.set_xlim(xmin, xmax)
    ax_inset.set_ylim(ymin, ymax)
    ax_inset.set_aspect("equal")
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    ax_inset.set_title("Sydney metropolitan\narea:", fontsize=8)
    mark_inset(ax, ax_inset, loc1=2, loc2=3, fc="none", ec="black", lw=0.8)
    add_panel_legend(ax, panel["label"], panel["legend_labels"], colours, to_add_to_legend=to_add_to_legend)
    # if panel has used gradient, let's plot a colour bar under the inset
    if draw_type == "gradient":
        labels = panel["legend_labels"]
        inset_pos = ax_inset.get_position()
        cbar_bounds = [
            inset_pos.x0,
            inset_pos.y0 - 0.035,
            inset_pos.width,
            0.018,
        ]
        cbar_ax = fig.add_axes(cbar_bounds)
        draw_gradient_cbar(cbar_ax, gradient_colours, title="Relative proportion", labels=labels)

def draw_gradient_cbar(ax: plt.Axes, colours: list[str], title="", labels=[]) -> None:
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    cmap = LinearSegmentedColormap.from_list("subpopulation_gradient", colours)
    ax.imshow(
        gradient,
        aspect="auto",
        cmap=cmap,
        extent=[0, 1, 0, 1],
    )
    ax.set_xticks([0.25, 0.5, 0.75])
    ax.set_xticklabels([0.25, 0.5, 0.75], fontsize=7)
    ax.tick_params(axis='both', which='major', pad=-7)
    ax.set_yticks([])
    ax.set_title(title, fontsize=7, pad=2)
    ax.grid(False)
    # annotate labels at either end
    ax.text(-0.05, 0.5, labels[1], ha='right', va='center', fontsize=8)
    ax.text(1.05, 0.5, labels[0], va='center', fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)

def plot_cryptic_subpopulations_overview(subpopulation_data: dict, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig = plt.figure(figsize=(8, 7.8))
    map_x_padding = -0.2  # increase/decrease this to move maps away from/towards the insets
    map_min_x, _, map_max_x, _ = subpopulation_data["map_boundary"].total_bounds
    map_xlim = (map_min_x, map_max_x + map_x_padding)

    panel_bounds = [
        [0.08, 0.49, 0.84, 0.40],
        [0.08, 0.10, 0.84, 0.40],
    ]
    inset_bounds = [
        [0.60, 0.64, 0.28, 0.14],
        [0.60, 0.25, 0.28, 0.14],
    ]
    for panel_bound, panel, inset_bound in zip(panel_bounds, subpopulation_data["panels"], inset_bounds):
        ax = fig.add_axes(panel_bound)
        plot_subpopulation_map_panel(fig, ax, subpopulation_data["map_boundary"], panel, inset_bound, to_add_to_legend=panel.get("to_add_to_legend"), draw_type=panel.get("draw_type"), map_xlim=map_xlim)

    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def plot_ae_notoscriptus_population_latitude_associations(panel: dict, output_plot_file: str, output_stats_file: str) -> None:
    points = panel["points"]
    association_configs = [
        {
            "column": "proportion_cluster1",
            "label": "population 1 vs populations 2/3",
            "ylabel": "population proportion",
            "colour": panel["colours"][0],
        },
        {
            "column": "proportion_cluster2_given_23",
            "label": "population 2 vs population 3",
            "ylabel": "population 2 proportion among populations 2/3",
            "colour": panel["colours"][1],
        },
    ]
    association_columns = [config["column"] for config in association_configs]
    plot_df = points[["sample", "latitude", "longitude"] + association_columns].copy()
    label_map = {config["column"]: config["label"] for config in association_configs}
    long_df = plot_df.melt(
        id_vars=["sample", "latitude", "longitude"],
        value_vars=association_columns,
        var_name="comparison_column",
        value_name="proportion",
    )
    long_df["comparison"] = long_df["comparison_column"].map(label_map)

    sns.set_theme(context="talk", style="whitegrid", font_scale=0.85)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9.2, 4.2), sharex=True, sharey=True)
    stats_rows = []
    for ax, config in zip(axes, association_configs):
        comparison = config["label"]
        panel_data = long_df[long_df["comparison"] == comparison].replace([np.inf, -np.inf], np.nan).dropna(subset=["latitude", "proportion"])
        fit = fit_linear_regression(panel_data["latitude"], panel_data["proportion"])
        spearman = stats.spearmanr(panel_data["latitude"], panel_data["proportion"])
        sns.scatterplot(
            data=panel_data,
            x="latitude",
            y="proportion",
            color=config["colour"],
            edgecolor="white",
            linewidth=0.4,
            s=42,
            alpha=0.85,
            ax=ax,
        )
        x_line, mean, lower, upper = fit.predict()
        ax.plot(x_line, mean, color="black", linewidth=1.6)
        ax.fill_between(x_line, lower, upper, color="black", alpha=0.15, linewidth=0)
        ax.set_title(comparison)
        ax.set_xlabel("latitude")
        ax.set_ylabel(config["ylabel"])
        ax.set_ylim(-0.04, 1.04)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(
            0.04,
            0.96,
            f"R$^2$ = {fit.rvalue ** 2:.2f}\nP {'< 0.0001' if fit.pvalue < 0.0001 else '= ' + format(fit.pvalue, '.3g')}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        stats_rows.extend(
            [
                {"comparison": comparison, "metric": "n", "value": len(panel_data)},
                {"comparison": comparison, "metric": "formula", "value": f"{comparison}_proportion = {fit.intercept:.8g} + {fit.slope:.8g} * latitude"},
                {"comparison": comparison, "metric": "intercept", "value": fit.intercept},
                {"comparison": comparison, "metric": "slope_per_latitude_degree", "value": fit.slope},
                {"comparison": comparison, "metric": "r_value", "value": fit.rvalue},
                {"comparison": comparison, "metric": "r_squared", "value": fit.rvalue ** 2},
                {"comparison": comparison, "metric": "p_value", "value": fit.pvalue},
                {"comparison": comparison, "metric": "std_error_slope", "value": fit.stderr},
                {"comparison": comparison, "metric": "spearman_rho", "value": spearman.statistic},
                {"comparison": comparison, "metric": "spearman_p_value", "value": spearman.pvalue},
            ]
        )
        existing_metrics = {row["metric"] for row in stats_rows if row["comparison"] == comparison}
        stats_rows.extend(
            {"comparison": comparison, "metric": key, "value": value}
            for key, value in fit.summary().items() if key not in existing_metrics
        )
    fig.tight_layout()
    fig.savefig(output_plot_file, dpi=300)
    plt.close(fig)
    pd.DataFrame(stats_rows).to_csv(output_stats_file, sep="\t", index=False)


def heatmap_Cx_annulirostris_subtype_kmer_detections(
    kmer_counts_df: pd.DataFrame,
    output_file: str,
    s_aus_proportions: pd.Series | None = None,
) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig = plt.figure(figsize=(5.5, 6))
    grid = fig.add_gridspec(
        nrows=1,
        ncols=5,
        width_ratios=[0.25, 1.6, 0.18, 0.25, 0.18],
        wspace=0.35,
    )
    prop_ax = fig.add_subplot(grid[0, 0])
    ax = fig.add_subplot(grid[0, 1])
    prop_cbar_ax = fig.add_subplot(grid[0, 2])
    kmer_cbar_ax = fig.add_subplot(grid[0, 4])

    if s_aus_proportions is None:
        s_aus_proportions = pd.Series(np.nan, index=kmer_counts_df.index)
    s_aus_proportions = s_aus_proportions.reindex(kmer_counts_df.index)
    sns.heatmap(
        s_aus_proportions.to_frame("S-AUS"),
        cmap="Greens",
        vmin=0,
        vmax=1,
        yticklabels=False,
        xticklabels=False,
        cbar=True,
        cbar_ax=prop_cbar_ax,
        cbar_kws={"label": ""},
        ax=prop_ax,
    )
    prop_ax.set_xlabel("")
    prop_ax.set_ylabel("")
    prop_cbar_ax.tick_params(labelsize=8)
    prop_cbar_ax.set_title("S-AUS\nproportion", fontsize=12)
    prop_ax.set_ylabel("sample")

    positive_kmer_counts = kmer_counts_df.where(kmer_counts_df > 0)
    max_kmer_count = positive_kmer_counts.max().max()
    sns.heatmap(
        kmer_counts_df,
        cmap=ListedColormap(["white"]),
        linewidths=0.5,
        linecolor="lightgray",
        cbar=False,
        ax=ax,
        mask=kmer_counts_df > 0,
        vmin=0,
        vmax=1,
    )
    kmer_heatmap = sns.heatmap(
        positive_kmer_counts,
        cmap=sns.color_palette("magma_r", as_cmap=True),
        norm=LogNorm(vmin=1, vmax=max(1, max_kmer_count)),
        linewidths=0.5,
        linecolor="lightgray",
        cbar=True,
        cbar_ax=kmer_cbar_ax,
        cbar_kws={"label": ""},
        ax=ax,
    )
    ax.set_xlabel("subtype-specific k-mer set")
    ax.set_ylabel("")
    ax.set_yticklabels([])

    subtype_kmer_counts = kmer_counts_df.attrs.get("subtype_kmer_counts", {})
    ax.set_xticklabels(
        [f"{subtype}\n(n = {subtype_kmer_counts.get(subtype, 0)})" for subtype in kmer_counts_df.columns],
        rotation=0,
        ha="center",
        fontsize=8,
    )
    for column_index in range(1, len(kmer_counts_df.columns)):
        ax.axvline(column_index, color="black", linewidth=1.5)
    ax.tick_params(axis="x", length=0)
    kmer_cbar_ticks = [tick for tick in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000] if tick <= max_kmer_count]
    if max_kmer_count not in kmer_cbar_ticks:
        kmer_cbar_ticks.append(max_kmer_count)
    kmer_cbar = kmer_heatmap.collections[-1].colorbar
    kmer_cbar.set_ticks(kmer_cbar_ticks)
    kmer_cbar.set_ticklabels([str(int(tick)) for tick in kmer_cbar_ticks])
    kmer_cbar_ax.tick_params(labelsize=8)
    kmer_cbar_ax.set_title("k-mer\ncounts", fontsize=12)

    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.12, top=0.9, wspace=0.35)
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_iqtree_tree(tree_file: str, output_file: str, label_metadata: str, subtype_colours: dict) -> None:
    label_metadata = pd.read_csv(label_metadata, sep="\t")

    label_dict = {}
    for _, row in label_metadata.iterrows():
        label_dict[row["ncbi_accession"]] = row["subtype"]

    def get_terminal_subtype(clade):
        ncbi_accession = clade.name.split("_")[0]
        return label_dict.get(ncbi_accession)

    def get_branch_colour(subtype):
        colour = subtype_colours.get(subtype, "black")
        if isinstance(colour, str) and colour.startswith("#") and len(colour) == 9:
            colour = colour[:7]
        return colour

    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    for clade in tree.find_clades(order="postorder"):
        subtypes = [
            get_terminal_subtype(terminal)
            for terminal in clade.get_terminals()
        ]
        if len(set(subtypes)) == 1 and subtypes[0] is not None:
            clade.color = get_branch_colour(subtypes[0])
        else:
            clade.color = "black"

    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None

    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda clade: None,
        branch_labels=lambda clade: None,
        show_confidence=False,
    )

    terminals = tree.get_terminals()
    y_by_clade = {clade: y for y, clade in enumerate(terminals, start=1)}
    first_root_child = tree.root.clades[0]
    child_ys = [y_by_clade[tip] for tip in first_root_child.get_terminals()]
    root_y = (min(child_ys) + max(child_ys)) / 2
    _, x_max = ax.get_xlim()
    stem_length = x_max * 0.04

    ax.hlines(
        y=root_y,
        xmin=-stem_length,
        xmax=0,
        color="black",
        linewidth=0.75,
        zorder=1,
    )

    ax.set_xlim(-stem_length * 1.2, x_max)

    for collection in ax.collections:
        collection.set_linewidth(0.75)

    for y, clade in enumerate(terminals, start=1):
        subtype = get_terminal_subtype(clade)
        colour = subtype_colours.get(subtype, "black")
        x = tree.distance(tree.root, clade)

        ax.scatter(
            x,
            y,
            s=20,
            color=colour,
            edgecolor="None",
            linewidth=0.2,
            zorder=5,
            clip_on=False,
        )

    subtypes_in_tree = []
    for clade in terminals:
        subtype = get_terminal_subtype(clade)
        if subtype is not None and subtype not in subtypes_in_tree:
            subtypes_in_tree.append(subtype)

    # for i, subtype in enumerate(subtypes_in_tree):
    #     text_colour = subtype_colours.get(subtype, "black")
    #     if subtype in ["PNG-1", "PNG-2", "SI"]:
    #         text_colour = "#666666"

    #     ax.text(
    #         1.02,
    #         0.9 - (i * 0.08),
    #         subtype,
    #         transform=ax.transAxes,
    #         color=text_colour,
    #         fontsize=10,
    #         va="center",
    #     )

    ax.grid(False)
    ax.set_axis_off()
    ax.set_xlabel("")
    ax.set_ylabel("")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cx_annulirostris_coastal_association(data: pd.DataFrame, output_plot_file: str, output_stats_file: str, *, metadata_file: str, input_files: list[str]) -> None:
    """Plot S-AUS proportion against coastal distance, with OLS mean CIs."""
    # log10 distance to coast 
    data["distance_to_coast_km_log10"] = np.log10(data["distance_to_coast_km"])
    fit = fit_linear_regression(data["distance_to_coast_km_log10"], data["avg_subtype_ratio"])
    x, mean, lower, upper = fit.predict()
    fig, ax = plt.subplots(figsize=(5, 4.2))
    # if diance to coast < 20, colour blue, if more than colour orange
    data["colour"] = np.where(data["distance_to_coast_km"] < 20, "#0072B2", "#F28E2B")
    ax.scatter(fit.x, fit.y, color=data["colour"], edgecolor="0.3", linewidth=0.4, s=35)
    ax.plot(x, mean, color="black", linewidth=1.6)
    ax.fill_between(x, lower, upper, color="black", alpha=0.15, linewidth=0)
    ax.set(xlabel="Distance from coast (km)", ylabel="S-AUS proportion")
    # ax.text(0.04, 0.96, f"R$^2$ = {fit.rvalue ** 2:.2f}\nP = {fit.pvalue:.3g}",
    #         transform=ax.transAxes, va="top")
    # label x axis in log10 scale, but show ticks in km
    ax.set_xticks(np.log10([0.1, 1, 5, 10, 20, 50, 100, 200, 500]))
    ax.set_xticklabels([0.1, 1, 5, 10, 20, 50, 100, 200, 500])
    ax.set_xlim(np.log10(0.05), np.log10(1000))
    print(min(data["distance_to_coast_km"]), max(data["distance_to_coast_km"]))
    sns.despine(ax=ax)
    ax.set_ylim(-0.04, 1.04)
    # add a vertical line at 20 km
    ax.axvline(np.log10(20), color="black", linestyle="--", linewidth=0.8)
    # annotate the regression formula, p val and r squared
    ax.text(0, 1.04, f"R$^2$ = {fit.rvalue ** 2:.2f}\nP = {fit.pvalue:.3g}\nS-AUS proportion = {fit.intercept:.3g} + {fit.slope:.3g} * log10(distance from coast)",
            transform=ax.transAxes, va="bottom", ha="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_plot_file, dpi=300)
    plt.close(fig)
    pearson = stats.pearsonr(fit.x, fit.y)
    spearman = stats.spearmanr(fit.x, fit.y)
    summary = dict(
        analysis="Cx_annulirostris_S_AUS_proportion_vs_distance_from_coast",
        input_files=";".join(input_files), metadata_file=metadata_file,
        filters="metadata samples only; min_count > 20; avg_subtype_ratio finite; distance_to_coast_km finite",
        response="S_AUS_proportion", predictor="distance_to_coast_km", response_transform="identity",
        formula=f"S_AUS_proportion = {fit.intercept:.8g} + {fit.slope:.8g} * distance_to_coast_km",
        **fit.summary(), slope_per_km=fit.slope, slope_per_100km=100 * fit.slope,
        pearson_r=pearson.statistic, pearson_p_value=pearson.pvalue,
        spearman_rho=spearman.statistic, spearman_p_value=spearman.pvalue,
        plot_file=output_plot_file,
    )
    with open(output_stats_file, "w") as handle:
        for name, value in summary.items():
            handle.write(f"{name}\t{value}\n")
