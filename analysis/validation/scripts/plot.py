# Plot abundance validation figures.

from decimal import Decimal, ROUND_HALF_UP

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import seaborn as sns
from matplotlib.colors import Normalize


GROUP_LABELS = [
    "Single species pools",
    r"Morphologically similar taxa: $\it{Cx.pipi/Cx.aust}$",
    "Laboratory mixtures",
    r"$\it{In}$ $\it{silico}$ mixtures",
    "Laboratory serial dilutions",
]

GROUP_ORDER = [
    "species",
    "morphology test",
    "laboratory mixture",
    "simulated mixture",
    "laboratory serial dilution",
]

BRACKET_GROUPS = [
    ("Identification", ["species", "morphology test"]),
    ("Quantification", ["laboratory mixture", "simulated mixture", "laboratory serial dilution"]),
]

GENUS_PALETTE = {
    "Aedes": "#4E6E5D",
    "Anopheles": "#B7C7A0",
    "Coquillettidia": "#D9A441",
    "Culex": "#6C8EBF",
    "Mansonia": "#8E6C8A",
}


def format_heatmap_row(values) -> list[str]:
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0)
    scaled_values = values * 100
    rounded_units = np.floor(scaled_values).astype(int)
    target_units = int(
        Decimal(str(values.sum() * 100)).quantize(Decimal("1"), rounding=ROUND_HALF_UP)
    )
    units_to_add = target_units - int(rounded_units.sum())
    if units_to_add > 0:
        remainders = scaled_values - rounded_units
        for index in np.argsort(-remainders, kind="stable")[:units_to_add]:
            rounded_units[index] += 1
    return [f"{units / 100:.2f}" for units in rounded_units]


def plot_heatmap_to_ax(data, labels, ax, cmap, cmap_obj, norm, vmin: float, vmax: float) -> None:
    sns.heatmap(
        data,
        ax=ax,
        cmap=cmap,
        cbar=False,
        vmin=vmin,
        vmax=vmax,
        yticklabels=labels,
    )
    for i in range(data.shape[0]):
        formatted_values = format_heatmap_row(data.iloc[i].to_numpy())
        for j in range(data.shape[1]):
            value = data.iloc[i, j]
            if value > 0.005:
                r, g, b, _ = cmap_obj(norm(value))
                luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
                text_color = "white" if luminance < 0.5 else "black"
                ax.text(j + 0.5, i + 0.5, formatted_values[j], ha="center", va="center", color=text_color, fontsize=7)

def plot_genus_bar_to_ax(column_metadata: list[dict], x_labels: list[str], ax) -> None:
    genera = [metadata["genus"] for metadata in column_metadata]
    unique_genera = list(dict.fromkeys(genera))
    genus_palette = {
        genus: GENUS_PALETTE.get(genus, "#A6A6A6")
        for genus in unique_genera
    }
    genus_index = {genus: index for index, genus in enumerate(unique_genera)}
    genus_values = [[genus_index[genus] for genus in genera]]
    ax.imshow(
        genus_values,
        aspect="auto",
        cmap=ListedColormap([genus_palette[genus] for genus in unique_genera]),
        extent=(0, len(genera), 0, 1),
        interpolation="nearest",
    )
    ax.grid(False)
    ax.set_xlim(0, len(genera))
    ax.set_xticks([index + 0.5 for index in range(len(genera))])
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontstyle="italic", fontsize=10)
    ax.set_yticks([])
    ax.tick_params(axis="x", pad=1, length=0)
    ax.tick_params(axis="y", left=False, labelleft=False)
    for _, spine in ax.spines.items():
        spine.set_visible(False)


def get_genus_boundaries(column_metadata: list[dict]) -> list[int]:
    boundaries = []
    genera = [metadata["genus"] for metadata in column_metadata]
    for index in range(1, len(genera)):
        if genera[index] != genera[index - 1]:
            boundaries.append(index)
    return boundaries


def add_genus_boundaries_to_ax(ax, boundaries: list[int], color: str = "#333333", linewidth: float = 0.8) -> None:
    for boundary in boundaries:
        ax.axvline(boundary, color=color, linewidth=linewidth, alpha=0.8)


def style_group_label_ax(ax, ylabel: str) -> None:
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    sns.despine(ax=ax, left=True, bottom=False, right=False, top=False)
    ax.set_ylabel(ylabel, rotation=0, ha="left", va="center", fontsize=11)
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_label_coords(1.25, 0.5)
    ax.yaxis.label.set_verticalalignment("center")


def add_panel_title(ax, title: str) -> None:
    ax.set_title(title, fontsize=11, loc="left", pad=3)


def move_genus_bar_close_to_final_panel(genus_ax, final_panel_ax, gap: float = 0.012) -> None:
    final_panel_position = final_panel_ax.get_position()
    genus_position = genus_ax.get_position()
    genus_ax.set_position([
        final_panel_position.x0,
        final_panel_position.y0 - gap - genus_position.height,
        final_panel_position.width,
        genus_position.height,
    ])


def style_bracket_ax(ax) -> None:
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for _, spine in ax.spines.items():
        spine.set_visible(False)


def add_group_brackets(fig, group_to_ax: dict, bracket_ax, bracket_groups: list[tuple[str, list[str]]]) -> None:
    bracket_position = bracket_ax.get_position()
    x0 = bracket_position.x0 + 0.2 * bracket_position.width
    x1 = bracket_position.x1
    for label, groups in bracket_groups:
        positions = [group_to_ax[group].get_position() for group in groups]
        y0 = min(position.y0 for position in positions)
        y1 = max(position.y1 for position in positions)
        y_mid = (y0 + y1) / 2
        for x_values, y_values in [
            ([x1, x1], [y0, y1]),
            ([x0, x1], [y0, y0]),
            ([x0, x1], [y1, y1]),
        ]:
            fig.add_artist(plt.Line2D(x_values, y_values, transform=fig.transFigure, color="black", linewidth=1))
        fig.text(x1 + 0.005, y_mid, label, rotation=0, va="center", ha="left", fontsize=11)


def style_heatmap_ax(ax, x_labels: list[str], show_x_labels: bool = False) -> None:
    ax.tick_params(axis="x", labelrotation=45, length=0)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.tick_params(axis="y", rotation=0, pad=-5)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.2)

    ax.set_yticklabels([label.get_text().strip() for label in ax.get_yticklabels()], fontsize=10)

    if show_x_labels:
        ax.set_xticklabels(x_labels, fontstyle="normal", fontsize=10)
    else:
        ax.tick_params(axis="x", labelbottom=False)


def plot_simulated_mixture_boxplot_to_ax(
    boxplot_data,
    column_metadata: list[dict],
    ax,
    show_x_labels: bool = False,
) -> None:
    n_species = len(column_metadata)
    positions = np.arange(n_species) + 0.5
    values_by_species = []
    for metadata in column_metadata:
        values = (
            boxplot_data
            .filter(boxplot_data["species"] == metadata["species"])
            .select("estimated_abundance")
            .to_series()
            .to_list()
        )
        values_by_species.append(values)

    ax.boxplot(
        values_by_species,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "#D9D9D2", "edgecolor": "#333333", "linewidth": 0.8},
        medianprops={"color": "#333333", "linewidth": 1.0},
        whiskerprops={"color": "#333333", "linewidth": 0.8},
        capprops={"color": "#333333", "linewidth": 0.8},
    )
    rng = np.random.default_rng(42)
    for position, values in zip(positions, values_by_species):
        if len(values) == 0:
            continue
        jitter = rng.uniform(-0.18, 0.18, len(values))
        ax.scatter(
            np.repeat(position, len(values)) + jitter,
            values,
            s=8,
            facecolors="none",
            edgecolors="#5F6368",
            alpha=0.25,
            linewidth=0.45,
            zorder=3,
        )

    ax.axhline(0.25, color="#333333", linestyle="--", linewidth=1.0, alpha=0.8)
    # add a shaded bar around the line to show 10% relative error
    ax.axhspan(0.225, 0.275, color="#889EBF", alpha=0.2, zorder=2, linewidth=1, edgecolor="#3A43E5")

    ax.set_xlim(0, n_species)
    ax.set_ylim(0.18, 0.32)
    ax.set_xticks(positions)
    if show_x_labels:
        ax.set_xticklabels(
            [metadata["species"] for metadata in column_metadata],
            rotation=45,
            ha="right",
            fontstyle="italic",
            fontsize=9,
        )
    else:
        ax.set_xticklabels([])
    ax.set_ylabel("estimated\nabundance", fontsize=9, rotation=0, labelpad=10, ha="right", va="center")
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=8, pad=-5)
    ax.grid(True, which="major", axis="y", linewidth=0.8, alpha=0.25)
    ax.grid(True, which="major", axis="x", linewidth=0, alpha=0)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.2)


def plot_kmer_error_scatter_to_ax(scatter_data, x_column: str, x_label: str, ax) -> None:
    data = scatter_data.to_pandas()
    right_label_threshold = data[x_column].quantile(0.9)
    ax.scatter(
        data[x_column],
        data["mean_signed_error"],
        s=24,
        facecolors="#6C8EBF",
        edgecolors="#333333",
        linewidth=0.35,
        alpha=0.75,
    )
    for _, row in data.iterrows():
        align_right = row[x_column] >= right_label_threshold
        genus, species = row["plot_label"].split(" ", 1)
        plot_label = f"{genus[0]}. {species}"
        ax.annotate(
            plot_label,
            (row[x_column], row["mean_signed_error"]),
            xytext=(-3 if align_right else 3, 2),
            textcoords="offset points",
            ha="right" if align_right else "left",
            fontsize=6,
            fontstyle="italic",
            alpha=0.9,
        )
    ax.axhline(0, color="#555555", linewidth=0.6, zorder=0)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("mean signed error in\nestimated abundance", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, which="major", linewidth=0.8, alpha=0.25)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.2)
    ylims = ax.get_ylim()
    ax.set_ylim(ylims[0], ylims[1] + (0.05 * (ylims[1] - ylims[0])))  # add a small amount of space above the top of the plot for labels


def plot_non_confident_assignment_error_scatter_to_ax(scatter_data, ax) -> None:
    data = scatter_data.to_pandas()
    right_label_threshold = data["non_confidently_assigned"].quantile(0.9)
    ax.scatter(
        data["non_confidently_assigned"],
        data["mean_signed_error"],
        s=24,
        facecolors="#6C8EBF",
        edgecolors="#333333",
        linewidth=0.35,
        alpha=0.75,
    )
    for _, row in data.iterrows():
        align_right = row["non_confidently_assigned"] >= right_label_threshold
        ax.annotate(
            row["plot_label"],
            (row["non_confidently_assigned"], row["mean_signed_error"]),
            xytext=(-3 if align_right else 3, 2),
            textcoords="offset points",
            ha="right" if align_right else "left",
            fontsize=6,
            fontstyle="italic",
            alpha=0.9,
        )
    ax.axhline(0, color="#555555", linewidth=0.6, zorder=0)
    ax.set_xlabel("proportion assigned to\nincorrect species", fontsize=9)
    ax.set_ylabel("mean signed error in\nestimated abundance", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, which="major", linewidth=0.8, alpha=0.25)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.2)
    ylims = ax.get_ylim()
    ax.set_ylim(ylims[0], ylims[1] + (0.05 * (ylims[1] - ylims[0])))


def plot_simulated_mixture_supplement(
    boxplot_data,
    scatter_data,
    column_metadata: list[dict],
    output_file: str,
) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(9, 4.5),
        gridspec_kw={"width_ratios": [2, 1], "wspace": 0.32},
    )
    plot_simulated_mixture_boxplot_to_ax(boxplot_data, column_metadata, axs[0], show_x_labels=True)
    axs[0].set_title("", fontsize=10, loc="left", pad=4)
    axs[0].set_xlabel("Species", fontsize=9)
    axs[0].set_ylabel("estimated abundance", fontsize=9, rotation=90, labelpad=10)
    plot_non_confident_assignment_error_scatter_to_ax(scatter_data, axs[1])
    axs[1].set_title("", fontsize=10, loc="left", pad=4)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.88, bottom=0.34, wspace=0.32)
    fig.savefig(output_file, dpi=300)


def plot_kmer_error_supplement(scatter_data, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, axs = plt.subplots(1, 3, figsize=(9, 3), gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.1})
    scatter_panels = [
        ("correct_weighted_proportion", "correct k-mer proportion\n(including shared)"),
        ("incorrect_only_weighted_proportion", "incorrect-only k-mer\nproportion"),
        ("absent_weighted_proportion", "absent k-mer\nproportion"),
    ]
    for ax, (x_column, x_label) in zip(axs, scatter_panels):
        plot_kmer_error_scatter_to_ax(scatter_data, x_column, x_label, ax)
    for ax in axs[1:]:
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelleft=False)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.94, bottom=0.22, wspace=0.32)
    fig.savefig(output_file, dpi=300)


def plot_in_silico_mixture_error_violins(violin_data, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    data = violin_data.to_pandas()
    fig, axs = plt.subplots(1, 2, figsize=(5.8, 3.2))
    panels = [
        ("rmse", "Mixture RMSE"),
        ("proportion_off_target", "Proportion off-target"),
    ]
    palette = {"KMA": "#6C8EBF", "kallisto": "#D9A441"}
    for ax, (column, title) in zip(axs, panels):
        sns.violinplot(
            data=data,
            x="tool",
            y=column,
            hue="tool",
            order=["KMA", "kallisto"],
            hue_order=["KMA", "kallisto"],
            palette=palette,
            cut=0,
            inner=None,
            linewidth=0.6,
            legend=False,
            ax=ax,
        )
        sns.boxplot(
            data=data,
            x="tool",
            y=column,
            order=["KMA", "kallisto"],
            width=0.22,
            showfliers=False,
            boxprops={"facecolor": "lightgrey", "edgecolor": "#333333", "linewidth": 0.7},
            whiskerprops={"color": "#333333", "linewidth": 0.7},
            capprops={"color": "#333333", "linewidth": 0.7},
            medianprops={"color": "#333333", "linewidth": 0.8},
            ax=ax,
        )
        sns.stripplot(
            data=data,
            x="tool",
            y=column,
            order=["KMA", "kallisto"],
            color="#222222",
            size=1.4,
            alpha=0.22,
            jitter=0.18,
            ax=ax,
        )
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("proportion", fontsize=9)
        ax.set_ylim(0, max(data[column].max() * 1.08, 0.01))
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True, axis="y", linewidth=0.8, alpha=0.25)
        ax.grid(False, axis="x")
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(0.2)
    # replace kallisto with Kallisto in axis labels
    axs[0].set_xticklabels(["KMA", "Kallisto"])
    axs[1].set_xticklabels(["KMA", "Kallisto"])
    fig.subplots_adjust(left=0.11, right=0.99, top=0.86, bottom=0.2, wspace=0.32)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)


def plot_control_abundance_heatmap(
    heatmap_data: dict,
    output_file: str | None = None,
    cmap: str = "Blues",
    vmin: float = 0,
    vmax: float = 1,
    simulated_mixture_panel: str = "boxplot",
    fig=None,
    subplot_spec=None,
    show_group_brackets: bool = True,
) -> None:
    standalone = fig is None
    if simulated_mixture_panel not in ["boxplot", "heatmap"]:
        raise ValueError("simulated_mixture_panel must be 'boxplot' or 'heatmap'")
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    cmap_obj = plt.get_cmap(cmap)
    norm = Normalize(vmin=vmin, vmax=vmax)

    column_metadata = heatmap_data["column_metadata"]
    x_labels = [metadata.get("display_species", metadata["species"]) for metadata in column_metadata]
    height_ratios = []
    for group in GROUP_ORDER:
        if group == "simulated mixture" and simulated_mixture_panel == "boxplot" and "simulated_mixture_boxplot_data" in heatmap_data:
            height_ratios.append(5)
        else:
            height_ratios.append(len(heatmap_data[group][1]))
    height_ratios = height_ratios + [0.28]
    n_panel_rows = len(GROUP_ORDER)
    n_rows = n_panel_rows + 1
    grid_kwargs = {
        "height_ratios": height_ratios,
        "hspace": 0.48,
        "wspace": 0.05,
        "width_ratios": [1, 0.05],
    }
    if standalone:
        fig = plt.figure(figsize=(7, 9.0))
        inner_grid = fig.add_gridspec(n_rows, 2, **grid_kwargs)
    else:
        if subplot_spec is None:
            raise ValueError("subplot_spec is required when drawing into an existing figure")
        inner_grid = subplot_spec.subgridspec(
            nrows=n_rows, ncols=2, **grid_kwargs,
        )
    axs = np.empty((n_rows, 2), dtype=object)
    shared_heatmap_ax = None
    for row in range(n_rows):
        grid_row = row
        axs[row, 0] = fig.add_subplot(
            inner_grid[grid_row, 0], sharex=shared_heatmap_ax,
        )
        if shared_heatmap_ax is None:
            shared_heatmap_ax = axs[row, 0]
        axs[row, 1] = fig.add_subplot(inner_grid[grid_row, 1])
    axs = np.atleast_2d(axs)

    heatmap_axes = [axs[row_index, 0] for row_index in range(n_panel_rows)]
    bracket_axes = [axs[row_index, 1] for row_index in range(n_panel_rows)]
    for group, ax in zip(GROUP_ORDER, heatmap_axes):
        if group == "simulated mixture" and simulated_mixture_panel == "boxplot" and "simulated_mixture_boxplot_data" in heatmap_data:
            plot_simulated_mixture_boxplot_to_ax(heatmap_data["simulated_mixture_boxplot_data"], column_metadata, ax)
        else:
            data, labels = heatmap_data[group]
            plot_heatmap_to_ax(data, labels, ax, cmap, cmap_obj, norm, vmin=vmin, vmax=vmax)

    for group, ax, title in zip(GROUP_ORDER, heatmap_axes, GROUP_LABELS):
        if group == "morphology test":
            heading, species = title.split(":", 1)
            heading_text = ax.set_title(heading + ":", fontsize=11, loc="left", pad=3)
            ax.annotate(
                species.strip(), xy=(1, 0), xycoords=heading_text,
                xytext=(3, 0), textcoords="offset points",
                fontsize=9, ha="left", va="bottom", annotation_clip=False,
            )
        else:
            add_panel_title(ax, title)

    for group, ax in zip(GROUP_ORDER, heatmap_axes):
        if group == "simulated mixture" and simulated_mixture_panel == "boxplot" and "simulated_mixture_boxplot_data" in heatmap_data:
            continue
        style_heatmap_ax(ax, x_labels, show_x_labels=False)

    for ax in bracket_axes:
        style_bracket_ax(ax)

    genus_ax = axs[-1, 0]
    style_bracket_ax(axs[-1, 1])
    plot_genus_bar_to_ax(column_metadata, x_labels, genus_ax)

    # add a x axis label along the bottom of the figure
    genus_ax.set_xlabel("Assigned species", fontsize=11, labelpad=2)

    if standalone:
        fig.tight_layout(h_pad=2)
        # add a small amount of whitespace vertically between subplots
        fig.subplots_adjust(left=0.28, right=0.82, top=0.96, bottom=0.15, hspace=0.28, wspace=0.02)
    move_genus_bar_close_to_final_panel(genus_ax, heatmap_axes[-1], gap=0.01)
    if show_group_brackets:
        add_group_brackets(fig, dict(zip(GROUP_ORDER, heatmap_axes)), bracket_axes[0], BRACKET_GROUPS)
    if standalone:
        if output_file is None:
            raise ValueError("output_file is required in standalone mode")
        fig.savefig(output_file, dpi=300)
    else:
        return fig, axs
