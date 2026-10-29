# make plots to assess results of in silico mixtures analysis

import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import ast

def plot_histogram_to_ax(data: pl.DataFrame, x: str, title: str, ax: plt.Axes, hue: str | None = None) -> None:
    plot_kwargs = {"x": x, "multiple": "stack", "bins": 20, "ax": ax}
    if hue is not None:
        plot_kwargs["hue"] = hue
    sns.histplot(data.to_pandas(), **plot_kwargs)
    ax.set_xlabel(x)
    style_ax(ax, grid_axis="y")

def plot_scatter_to_ax(data: pl.DataFrame, x: str, y: str, title: str, ax: plt.Axes, hue: str | None = None, style: str | None = None) -> None:
    plot_kwargs = {"data": data.to_pandas(), "x": x, "y": y, "ax": ax}
    if hue is not None:
        plot_kwargs["hue"] = hue
    if style is not None:
        plot_kwargs["style"] = style
    sns.scatterplot(**plot_kwargs)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    style_ax(ax)

def style_ax(ax: plt.Axes, grid_axis: str = "both") -> None:
    ax.grid(True, which="major", axis=grid_axis, linewidth=0.8, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def format_legend(ax: plt.Axes) -> None:
    legend = ax.get_legend()
    if legend is not None:
        legend.set_frame_on(False)

def plot_RMSE_hist_per_tool(data: pl.DataFrame, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(12, 6), sharex=True)
    plot_histogram_to_ax(data, x="kma_rmse", ax=ax[0], title="KMA RMSE")
    plot_histogram_to_ax(data, x="kallisto_rmse", ax=ax[1], title="Kallisto RMSE")
    ax[0].set_xlabel("KMA RMSE")
    ax[1].set_xlabel("kallisto RMSE")
    plt.tight_layout()
    plt.savefig(output_file)

def plot_propn_off_target_hist_per_tool(data: pl.DataFrame, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(12, 6), sharex=True)
    plot_histogram_to_ax(data, x="kma_proportion_off_target", ax=ax[0], title="KMA Proportion Off Target")
    plot_histogram_to_ax(data, x="kallisto_proportion_off_target", ax=ax[1], title="Kallisto Proportion Off Target")
    ax[0].set_xlabel("KMA proportion off target")
    ax[1].set_xlabel("kallisto proportion off target")
    plt.tight_layout()
    plt.savefig(output_file)

def plot_template_identity_and_coverage_hist_in_true_and_offtarget_per_tool(data: pl.DataFrame, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(12, 6))
    plot_histogram_to_ax(data, x="kma_template_identity", hue="is_true_target", ax=ax[0], title="KMA Template Identity")
    plot_histogram_to_ax(data, x="kma_template_coverage", hue="is_true_target", ax=ax[1], title="KMA Template Coverage")
    ax[0].set_xlabel("KMA template identity")
    ax[1].set_xlabel("KMA template coverage")
    for current_ax in ax:
        format_legend(current_ax)
    plt.tight_layout()
    plt.savefig(output_file)

def plot_scatter_template_identity_vs_coverage_coloured_by_true_target(data: pl.DataFrame, output_file: str) -> None:
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_scatter_to_ax(data, x="kma_identity_across_covered_length", y="kma_template_coverage", title="KMA Template Identity vs Coverage", ax=ax, hue="is_true_target", style="assigned_by_both_tools")
    ax.set_xlabel("KMA identity across covered length")
    ax.set_ylabel("KMA template coverage")
    format_legend(ax)
    # ax.set_xlim(95, 101)
    # ax.set_ylim(95, 101)
    plt.tight_layout()
    plt.savefig(output_file)

def plot_boxplots_of_abundance_estimate_per_species(error_data: pl.DataFrame, mixture_info_data: pl.DataFrame, output_file: str) -> None:
    # gather data
    n_components = len([i for i in mixture_info_data.columns if "proportion" in i])
    boxplot_species = []
    boxplot_estimated_proportions = []
    for mixture in error_data['mixture_id']:
        true_species = [mixture_info_data.filter(pl.col("mixture_id") == mixture)[f"sample_{i}"].item() for i in range(1, n_components + 1)]
        estimated_proportions = ast.literal_eval(error_data.filter(pl.col("mixture_id") == mixture).select("kma_estimated_proportions").item())
        for species, est_prop in zip(true_species, estimated_proportions):
            boxplot_species.append(species)
            boxplot_estimated_proportions.append(est_prop)
    boxplot_df = pl.DataFrame({"species": boxplot_species, "estimated_proportion": boxplot_estimated_proportions})
    # plot
    sns.set_theme(context="talk", style="whitegrid", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=boxplot_df.to_pandas(), x="species", y="estimated_proportion", ax=ax)
    ax.set_xlabel("species")
    ax.set_ylabel("estimated proportion")
    ax.axhline(0.25, color="k", linestyle="--", label="true proportion")
    ax.tick_params(axis='x', rotation=45)
    xtick_labels = ax.get_xticklabels()
    xtick_labels = [label.get_text().replace("_", " ") for label in xtick_labels]
    ax.set_xticklabels(xtick_labels, rotation=35, ha="right", rotation_mode="anchor")
    format_legend(ax)
    # add scatter points on top of boxplot
    sns.stripplot(data=boxplot_df.to_pandas(), x="species", y="estimated_proportion", ax=ax, color="k", alpha=0.5)
    # save
    plt.tight_layout()
    plt.savefig(output_file)

def plot_in_silico_mixtures_results(error_file: str, qunatification_file: str, mixtures_info_file: str, output_dir: str) -> None:
    error_data = pl.read_csv(error_file, separator="\t")
    quantification_data = pl.read_csv(qunatification_file, separator="\t")
    mixtures_info_data = pl.read_csv(mixtures_info_file, separator="\t")
    # plot error histograms
    plot_RMSE_hist_per_tool(error_data, f"{output_dir}/RMSE_histograms.png")
    plot_propn_off_target_hist_per_tool(error_data, f"{output_dir}/propn_off_target_histograms.png")
    sample_columns = [column for column in mixtures_info_data.columns if column.startswith("sample_")]
    true_targets = (
        mixtures_info_data.select(["mixture_id", pl.concat_list(sample_columns).alias("species")])
        .explode("species")
        .drop_nulls("species")
        .with_columns(pl.lit(True).alias("is_true_target"))
    )
    quantification_data = (
        quantification_data.join(true_targets, on=["mixture_id", "species"], how="left")
        .with_columns(pl.col("is_true_target").fill_null(False))
    )
    # add a col to qunatification data where if both kma and kallisto have more than 5 reads then is true, otherwise false
    quantification_data = quantification_data.with_columns(
        pl.when((pl.col("kma_read_count") > 5) & (pl.col("kallisto_estimated_counts") > 5))
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("assigned_by_both_tools")
    )
    # compute identity across covered length as (kma_template_identity / min(100, template_coverage)) * 100
    quantification_data = quantification_data.with_columns(
        (
        pl.col("kma_template_identity")
        / pl.min_horizontal(pl.lit(100.0), pl.col("kma_template_coverage")).clip(lower_bound=1e-9)
        * 100
        ).alias("kma_identity_across_covered_length")    
    )
    # if there is more than one target for a species in a mixture, keep the one with the highest template coverage and identity
    quantification_data = (
        quantification_data.sort(["mixture_id", "species", "kma_template_coverage", "kma_template_identity"], descending=[False, False, True, True])
        .group_by(["mixture_id", "species"])
        .head(1)
    )
    # print row where is true target but identity across covered length is less than 97.5
    print(quantification_data.filter((pl.col("is_true_target") == True) & (pl.col("kma_identity_across_covered_length") < 97.5)).select(["mixture_id", "species", "kma_template_identity", "kma_template_coverage", "kma_identity_across_covered_length"]))
    plot_template_identity_and_coverage_hist_in_true_and_offtarget_per_tool(quantification_data, f"{output_dir}/template_identity_and_coverage_histograms.png")
    plot_scatter_template_identity_vs_coverage_coloured_by_true_target(quantification_data, f"{output_dir}/template_identity_vs_coverage.png")
    # plot boxplots of estimated proportions per species
    plot_boxplots_of_abundance_estimate_per_species(error_data, mixtures_info_data, f"{output_dir}/estimated_proportions_per_species.png")
