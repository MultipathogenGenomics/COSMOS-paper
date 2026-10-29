import itertools
import math

import polars as pl
from matplotlib import pyplot as plt
import seaborn as sns

from analysis.common.mosquito_ids import shorted_genus


AXIS_LABEL_SIZE = 11
TICK_LABEL_SIZE = 8


def style_ax(ax: plt.Axes) -> None:
    ax.grid(True, which="major", axis="both", linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)


def plot_surveillance_validation(validation_df: pl.DataFrame, output_file: str) -> None:
    species = validation_df.select("species").unique().sort("species").to_series().to_list()
    ncols = len(species)
    fig_width = max(4.0 * ncols, 4.0)
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(fig_width, 4.0), squeeze=False)
    axs = axs.flatten()

    plot_df = validation_df.to_pandas()
    for ax, species_name in zip(axs, species):
        species_df = plot_df[plot_df["species"] == species_name]
        sns.scatterplot(
            data=species_df,
            x="surveillance_proportion",
            y="estimated_abundance",
            ax=ax,
            color="#5F6368",
            s=34,
            alpha=0.7,
            edgecolor="white",
            linewidth=0.25,
        )
        ax.plot([0, 1], [0, 1], color="#333333", linestyle="--", linewidth=1, alpha=0.55)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("NSW Health proportion", fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel("COSMOS estimated abundance", fontsize=AXIS_LABEL_SIZE)
        ax.set_title(shorted_genus(species_name), fontsize=AXIS_LABEL_SIZE, fontstyle="italic")
        ax.text(
            0.04,
            0.94,
            f"n = {len(species_df)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=TICK_LABEL_SIZE,
        )
        style_ax(ax)

    for ax in axs[len(species):]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def get_top_abundance_species(abundance_df: pl.DataFrame, top_n: int) -> list[str]:
    species_cols = [col for col in abundance_df.columns if col not in ["sample", "batch", "non_confidently_assigned"]]
    species_sums = abundance_df.select([pl.col(col).sum().alias(col) for col in species_cols]).row(0, named=True)
    return sorted(species_sums, key=species_sums.get, reverse=True)[:top_n]


def plot_pairwise_species_abundance_scatter(
    abundance_df: pl.DataFrame,
    output_file: str,
    top_n: int = 10,
    ncols: int = 5,
) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    species = get_top_abundance_species(abundance_df, top_n)
    species_pairs = list(itertools.combinations(species, 2))
    nrows = math.ceil(len(species_pairs) / ncols)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 2.8, nrows * 2.5), squeeze=False)
    axs = axs.flatten()

    plot_df = abundance_df.select(["sample"] + species).to_pandas()
    label_map = {species_name: shorted_genus(species_name) for species_name in species}

    for ax, (x_species, y_species) in zip(axs, species_pairs):
        sns.scatterplot(
            data=plot_df,
            x=x_species,
            y=y_species,
            ax=ax,
            color="#5F6368",
            s=18,
            alpha=0.55,
            edgecolor="white",
            linewidth=0.15,
        )
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel(label_map[x_species], fontsize=9, fontstyle="italic")
        ax.set_ylabel(label_map[y_species], fontsize=9, fontstyle="italic")
        style_ax(ax)

    for ax in axs[len(species_pairs):]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")


def write_pairwise_species_abundance_description(
    output_file: str,
    plot_output_file: str,
    assignment_reference: str,
    abundance_filename: str,
    top_n: int,
    abundance_df: pl.DataFrame,
) -> None:
    species = get_top_abundance_species(abundance_df, top_n)
    with open(output_file, "w") as handle:
        handle.write(f"plot_output_file\t{plot_output_file}\n")
        handle.write(f"assignment_reference\t{assignment_reference}\n")
        handle.write(f"abundance_filename\t{abundance_filename}\n")
        handle.write(f"top_n\t{top_n}\n")
        handle.write(f"n_species\t{len(species)}\n")
        handle.write(f"species\t{'; '.join(species)}\n")
        handle.write(f"n_samples\t{abundance_df.height}\n")
