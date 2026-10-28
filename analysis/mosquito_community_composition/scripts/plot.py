# Plot abundance geospatial figures.

import numpy as np
import polars as pl
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, to_rgba
import seaborn as sns


COASTAL_PALETTE = {
    1: "#8EBFE8",
    0: "#8FD19E",
    None: "#D9D9D9",
}

URBAN_RURAL_PALETTE = {
    "Urban": "#EFBA74",
    "Met Fringe": "#E38975",
    "Rural": "#7A1F5C",
}

URBAN_RURAL_METADATA_PALETTE = {
    "Urban": URBAN_RURAL_PALETTE["Urban"],
    "Met_Fringe": URBAN_RURAL_PALETTE["Met Fringe"],
    "Rural": URBAN_RURAL_PALETTE["Rural"],
}

COUNCIL_TYPE_PALETTE = {
    "Regional Town/City": URBAN_RURAL_PALETTE["Urban"],
    "Metropolitan": URBAN_RURAL_PALETTE["Urban"],
    "Rural": URBAN_RURAL_PALETTE["Rural"],
    "Large Rural": URBAN_RURAL_PALETTE["Rural"],
    "Metropolitan Fringe": URBAN_RURAL_PALETTE["Met Fringe"],
    "Missing": "#D9D9D9",
}

AXIS_LABEL_SIZE = 11
TICK_LABEL_SIZE = 8
ANNOTATION_LABEL_SIZE = 9
LEGEND_LABEL_SIZE = 10


def style_ax(ax: plt.Axes, grid_axis: str = "both") -> None:
    ax.grid(True, which="major", axis=grid_axis, linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)


def draw_colourbar_by_metadata(
    ax: plt.Axes,
    metadata: pl.DataFrame,
    col: str,
    cmap=None,
    reference_ax: plt.Axes | None = None,
    label_map: dict | None = None,
    palette: dict | None = None,
) -> dict:
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


def plot_abundance_heatmap_to_ax(heatmap_data: dict, ax: plt.Axes) -> None:
    sns.heatmap(
        heatmap_data["heatmap_data"],
        cmap="Greys",
        cbar=False,
        yticklabels=heatmap_data["species_labels"],
        xticklabels=False,
        ax=ax,
    )
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", rotation=0, labelsize=TICK_LABEL_SIZE)
    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")
    ax.set_ylabel("Mosquito species", fontsize=AXIS_LABEL_SIZE)
    ax.set_xlabel("")


def plot_species_distribution_group_ax(heatmap_ax: plt.Axes, group_ax: plt.Axes, n_species: int) -> None:
    group_ax.set_xlim(0, 1)
    group_ax.set_ylim(0, 1)
    group_ax.set_facecolor("white")
    group_ax.set_xticks([])
    group_ax.set_yticks([])
    group_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for _, spine in group_ax.spines.items():
        spine.set_visible(False)

    top_species_count = 6
    boundary = 1 - (top_species_count / n_species)
    bracket_gap = 0.018
    x_start = 0.4
    x_end = 1.0
    bracket_kwargs = {"color": "#333333", "linewidth": 1.5, "clip_on": False, "solid_capstyle": "butt"}
    group_ax.plot([x_start, x_end], [1, 1], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.plot([x_end, x_end], [boundary + bracket_gap, 1], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.plot([x_start, x_end], [boundary + bracket_gap, boundary + bracket_gap], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.plot([x_start, x_end], [boundary - bracket_gap, boundary - bracket_gap], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.plot([x_end, x_end], [0, boundary - bracket_gap], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.plot([x_start, x_end], [0, 0], transform=group_ax.transAxes, **bracket_kwargs)
    group_ax.text(1.15, (1 + boundary + bracket_gap) / 2, "widely\ndistributed", ha="left", va="center", fontsize=ANNOTATION_LABEL_SIZE, rotation=0, transform=group_ax.transAxes)
    group_ax.text(1.15, (boundary - bracket_gap) / 2, "coastal", ha="left", va="center", fontsize=ANNOTATION_LABEL_SIZE, rotation=0, transform=group_ax.transAxes)


def plot_umap_to_ax(embedding_df: pl.DataFrame, ax: plt.Axes) -> None:
    plot_df = embedding_df.to_pandas()
    plot_df["Coastal_status"] = plot_df["Coastal"].map({1: "coastal", 0: "inland"}).fillna("missing")
    plot_df["urban_rural"] = plot_df["urban_rural"].str.replace("_", " ").fillna("Missing")
    urban_rural_order = [label for label in ["Urban", "Met Fringe", "Rural"] if label in set(plot_df["urban_rural"])]
    coastal_order = [label for label in ["coastal", "inland"] if label in set(plot_df["Coastal_status"])]
    sns.scatterplot(
        data=plot_df,
        x="UMAP1",
        y="UMAP2",
        hue="urban_rural",
        style="Coastal_status",
        palette=URBAN_RURAL_PALETTE,
        hue_order=urban_rural_order,
        style_order=coastal_order,
        markers={"coastal": "o", "inland": "s"},
        ax=ax,
        s=55,
        alpha=0.85,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.set_xlabel("UMAP1", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("UMAP2", fontsize=AXIS_LABEL_SIZE)
    style_ax(ax)
    legend = ax.get_legend()
    if legend is not None:
        handles, labels = ax.get_legend_handles_labels()
        filtered_handles = []
        filtered_labels = []
        for handle, label in zip(handles, labels):
            if label in ["urban_rural", "Coastal_status", "Missing", "missing"]:
                continue
            filtered_handles.append(handle)
            filtered_labels.append(label.lower())
        ax.legend(filtered_handles, filtered_labels, loc="upper right", bbox_to_anchor=(-0.28, 1.0), frameon=False, fontsize=LEGEND_LABEL_SIZE)


def plot_pca_to_ax(embedding_df: pl.DataFrame, ax: plt.Axes) -> None:
    plot_df = embedding_df.to_pandas()
    plot_df["Coastal_status"] = plot_df["Coastal"].map({1: "coastal", 0: "inland"}).fillna("missing")
    plot_df["urban_rural"] = plot_df["urban_rural"].str.replace("_", " ").fillna("Missing")
    urban_rural_order = [label for label in ["Urban", "Met Fringe", "Rural"] if label in set(plot_df["urban_rural"])]
    coastal_order = [label for label in ["coastal", "inland"] if label in set(plot_df["Coastal_status"])]
    sns.scatterplot(
        data=plot_df,
        x="PCA1",
        y="PCA2",
        hue="urban_rural",
        style="Coastal_status",
        palette=URBAN_RURAL_PALETTE,
        hue_order=urban_rural_order,
        style_order=coastal_order,
        markers={"coastal": "o", "inland": "s"},
        ax=ax,
        s=55,
        alpha=0.85,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.set_xlabel("PC1", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("PC2", fontsize=AXIS_LABEL_SIZE)
    style_ax(ax)
    legend = ax.get_legend()
    if legend is not None:
        handles, labels = ax.get_legend_handles_labels()
        filtered_handles = []
        filtered_labels = []
        for handle, label in zip(handles, labels):
            if label in ["urban_rural", "Coastal_status", "Missing", "missing"]:
                continue
            filtered_handles.append(handle)
            filtered_labels.append(label.lower())
        ax.legend(filtered_handles, filtered_labels, loc="upper right", bbox_to_anchor=(-0.28, 1.0), frameon=False, fontsize=LEGEND_LABEL_SIZE)


def plot_embedding_comparison(
    embedding_df: pl.DataFrame,
    colour_specs: list[dict],
    output_file: str,
) -> None:
    plot_df = embedding_df.to_pandas()
    fig, axs = plt.subplots(
        2,
        len(colour_specs),
        figsize=(3.2 * len(colour_specs), 6.4),
        sharex="row",
        sharey="row",
        constrained_layout=False,
    )
    coordinate_pairs = [("PCA1", "PCA2"), ("UMAP1", "UMAP2")]

    for column_index, spec in enumerate(colour_specs):
        column = spec["column"]
        for row_index, (x_column, y_column) in enumerate(coordinate_pairs):
            ax = axs[row_index, column_index]
            missing = plot_df[column].isna()
            ax.scatter(
                plot_df.loc[missing, x_column],
                plot_df.loc[missing, y_column],
                s=28,
                color="#D9D9D9",
                edgecolor="white",
                linewidth=0.25,
                alpha=0.8,
                label="missing",
            )

            if spec["type"] == "categorical":
                for value in spec["order"]:
                    selected = (~missing) & (plot_df[column] == value)
                    ax.scatter(
                        plot_df.loc[selected, x_column],
                        plot_df.loc[selected, y_column],
                        s=30,
                        color=spec["palette"][value],
                        edgecolor="white",
                        linewidth=0.25,
                        alpha=0.85,
                        label=spec.get("labels", {}).get(value, str(value)),
                    )
                if row_index == 1:
                    handles, labels = ax.get_legend_handles_labels()
                    if not missing.any():
                        handles = handles[1:]
                        labels = labels[1:]
                    ax.legend(handles, labels, bbox_to_anchor=(0.5, -0.2), loc="upper center", frameon=False, fontsize=LEGEND_LABEL_SIZE)
            else:
                norm = Normalize(vmin=spec.get("vmin"), vmax=spec.get("vmax"))
                selected = ~missing
                scatter = ax.scatter(
                    plot_df.loc[selected, x_column],
                    plot_df.loc[selected, y_column],
                    c=plot_df.loc[selected, column],
                    cmap=spec["cmap"],
                    norm=norm,
                    s=30,
                    edgecolor="white",
                    linewidth=0.25,
                    alpha=0.85,
                )
                if row_index == 1:
                    cax = ax.inset_axes([0.10, -0.3, 0.80, 0.035])
                    cbar = fig.colorbar(
                        scatter,
                        cax=cax,
                        orientation="horizontal",
                    )
                    cbar.set_label(spec["label"], fontsize=9)
                    cbar.ax.tick_params(labelsize=8)
            ax.set_xlabel("PC1" if row_index == 0 else "UMAP1", fontsize=9)
            if column_index == 0:
                ax.set_ylabel("PC2" if row_index == 0 else "UMAP2", fontsize=9)
            style_ax(ax)
        axs[0, column_index].set_title(spec["label"], fontsize=10)
    fig.subplots_adjust(
        left=0.07,
        right=0.98,
        top=0.92,
        bottom=0.25,
        wspace=0.20,
        hspace=0.25,
    )
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sdi_boxplot_to_ax(sdi_df: pl.DataFrame, ax: plt.Axes) -> None:
    council_order = (
        sdi_df.group_by("council_type")
        .agg(pl.col("Simpsons_Diversity_Index").mean())
        .sort("Simpsons_Diversity_Index", descending=True)["council_type"]
        .to_list()
    )
    sdi_pandas = sdi_df.to_pandas()
    sns.boxplot(
        data=sdi_pandas,
        x="council_type",
        y="Simpsons_Diversity_Index",
        ax=ax,
        order=council_order,
        palette=COUNCIL_TYPE_PALETTE,
        width=0.4,
        showfliers=False,
    )
    sns.stripplot(
        data=sdi_pandas,
        x="council_type",
        y="Simpsons_Diversity_Index",
        ax=ax,
        order=council_order,
        color="black",
        alpha=0.4,
        size=4,
        jitter=True,
    )
    ax.set_xlabel("Council type", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Diversity Index (1 - D)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylim(-0.05, 0.95)
    label_map = {
        "Regional Town/City": "Regional\nTown/City",
        "Metropolitan": "Met",
        "Rural": "Rural",
        "Large Rural": "Large\nRural",
        "Metropolitan Fringe": "Met\nFringe",
        "Missing": "Missing",
    }
    ax.set_xticklabels([label_map.get(label.get_text(), label.get_text().replace(" ", "\n")) for label in ax.get_xticklabels()], fontsize=TICK_LABEL_SIZE)
    ax.tick_params(axis="x", rotation=0, pad=0)
    # for tick_label in ax.get_xticklabels():
    #     tick_label.set_ha("right")
    style_ax(ax, grid_axis="y")


def plot_abundance_geospatial_overview(geospatial_data: dict, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig = plt.figure(figsize=(9, 7))
    grid = fig.add_gridspec(nrows=4, ncols=2, height_ratios=[1.5, 0.18, 0.18, 2.1], hspace=0.12, wspace=0.3)

    heatmap_grid = grid[0, :].subgridspec(1, 2, width_ratios=[1, 0.045], wspace=0)
    coastal_grid = grid[1, :].subgridspec(1, 2, width_ratios=[1, 0.045], wspace=0)
    urban_grid = grid[2, :].subgridspec(1, 2, width_ratios=[1, 0.045], wspace=0)
    heatmap_ax = fig.add_subplot(heatmap_grid[0, 0])
    species_group_ax = fig.add_subplot(heatmap_grid[0, 1])
    coastal_bar_ax = fig.add_subplot(coastal_grid[0, 0])
    coastal_blank_ax = fig.add_subplot(coastal_grid[0, 1])
    urban_bar_ax = fig.add_subplot(urban_grid[0, 0])
    urban_blank_ax = fig.add_subplot(urban_grid[0, 1])
    embedding_ax = fig.add_subplot(grid[3, 0])
    sdi_ax = fig.add_subplot(grid[3, 1])

    plot_abundance_heatmap_to_ax(geospatial_data["heatmap"], heatmap_ax)
    plot_species_distribution_group_ax(heatmap_ax, species_group_ax, len(geospatial_data["heatmap"]["species_labels"]))
    coastal_blank_ax.axis("off")
    urban_blank_ax.axis("off")
    metadata = geospatial_data["heatmap"]["metadata"]
    draw_colourbar_by_metadata(coastal_bar_ax, metadata, "Coastal", reference_ax=heatmap_ax, label_map={1: "coastal", 0: "inland"}, palette=COASTAL_PALETTE)
    draw_colourbar_by_metadata(urban_bar_ax, metadata, "urban_rural", reference_ax=heatmap_ax, label_map={"Met_Fringe": "met fringe", "Rural": "rural", "Urban": "urban"}, palette=URBAN_RURAL_METADATA_PALETTE)
    if geospatial_data["embedding_method"] == "PCA":
        plot_pca_to_ax(geospatial_data["embedding"], embedding_ax)
    else:
        plot_umap_to_ax(geospatial_data["embedding"], embedding_ax)
    plot_sdi_boxplot_to_ax(geospatial_data["sdi"], sdi_ax)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def plot_pca_loadings(loading_df: pl.DataFrame, output_file: str) -> None:
    if loading_df is None:
        return
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    pc_cols = [col for col in loading_df.columns if col.startswith("PC")]
    if "PC1" not in pc_cols:
        return
    species_order = loading_df.sort("PC1", descending=True)["species"].to_list()
    plot_df = loading_df.to_pandas()
    plot_df["species"] = plot_df["species"].map(lambda species: species.replace("Ochlerotatus ", "Aedes "))

    fig, axs = plt.subplots(nrows=len(pc_cols), ncols=1, figsize=(9, 2.2 * len(pc_cols)), sharex=True)
    if len(pc_cols) == 1:
        axs = [axs]
    for ax, pc_col in zip(axs, pc_cols):
        sns.barplot(
            data=plot_df,
            x="species",
            y=pc_col,
            order=[species.replace("Ochlerotatus ", "Aedes ") for species in species_order],
            ax=ax,
            color="#5F6368",
            edgecolor="white",
            linewidth=0.2,
        )
        ax.axhline(0, color="#333333", linewidth=0.8)
        ax.set_xlabel("")
        ax.set_ylabel(pc_col, fontsize=AXIS_LABEL_SIZE)
        style_ax(ax, grid_axis="y")
    axs[-1].tick_params(axis="x", rotation=45, labelsize=TICK_LABEL_SIZE)
    for tick_label in axs[-1].get_xticklabels():
        tick_label.set_ha("right")
        tick_label.set_fontstyle("italic")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def plot_pca_variance_explained(variance_df: pl.DataFrame, output_file: str) -> None:
    if variance_df is None:
        return
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    plot_df = variance_df.to_pandas()
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(
        plot_df["component_number"],
        plot_df["explained_variance_ratio"],
        marker="o",
        color="#5F6368",
        linewidth=2,
        label="per PC",
    )
    ax.plot(
        plot_df["component_number"],
        plot_df["cumulative_explained_variance_ratio"],
        marker="o",
        color="#2F7EA1",
        linewidth=2,
        label="cumulative",
    )
    ax.set_xticks(plot_df["component_number"])
    ax.set_xticklabels(plot_df["principal_component"])
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("principal component", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("variance explained", fontsize=AXIS_LABEL_SIZE)
    style_ax(ax)
    ax.legend(frameon=False, fontsize=TICK_LABEL_SIZE)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def plot_pc1_distance_to_coast_association(embedding_df: pl.DataFrame, output_file: str) -> None:
    if "PCA1" not in embedding_df.columns or "distance_to_coast_km" not in embedding_df.columns:
        return
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    plot_df = embedding_df.select(
        [
            "PCA1",
            "distance_to_coast_km",
            "Coastal",
            "urban_rural",
        ]
    ).drop_nulls().to_pandas()
    plot_df["Coastal_status"] = plot_df["Coastal"].map({1: "coastal", 0: "inland"}).fillna("missing")
    plot_df["urban_rural"] = plot_df["urban_rural"].str.replace("_", " ").fillna("Missing")
    urban_rural_order = [label for label in ["Urban", "Met Fringe", "Rural"] if label in set(plot_df["urban_rural"])]
    coastal_order = [label for label in ["coastal", "inland"] if label in set(plot_df["Coastal_status"])]

    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    sns.regplot(
        data=plot_df,
        x="distance_to_coast_km",
        y="PCA1",
        scatter=False,
        ax=ax,
        color="#333333",
        line_kws={"linewidth": 1.5},
    )
    sns.scatterplot(
        data=plot_df,
        x="distance_to_coast_km",
        y="PCA1",
        hue="urban_rural",
        style="Coastal_status",
        palette=URBAN_RURAL_PALETTE,
        hue_order=urban_rural_order,
        style_order=coastal_order,
        markers={"coastal": "o", "inland": "s"},
        ax=ax,
        s=55,
        alpha=0.85,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.set_xlabel("distance to coast (km)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("PC1", fontsize=AXIS_LABEL_SIZE)
    style_ax(ax)
    legend = ax.get_legend()
    if legend is not None:
        handles, labels = ax.get_legend_handles_labels()
        filtered_handles = []
        filtered_labels = []
        for handle, label in zip(handles, labels):
            if label in ["urban_rural", "Coastal_status", "Missing", "missing"]:
                continue
            filtered_handles.append(handle)
            filtered_labels.append(label.lower())
        ax.legend(filtered_handles, filtered_labels, loc="best", frameon=False, fontsize=LEGEND_LABEL_SIZE)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
