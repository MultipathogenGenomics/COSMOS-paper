import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import os
import json

def parse_data():
    # parse files
    nc2_post_filt_true = "/Users/afor0966/Documents/code/COSMOS-analysis/data/castanetv9.2-output/COSMOS1-4_castanetv9.2_nc2.tsv"
    reads_dedup_post_filt_true = "/Users/afor0966/Documents/code/COSMOS-analysis/data/castanetv9.2-output/COSMOS1-4_castanetv9.2_reads_dedup.tsv"
    nc2_post_filt_false = "/Users/afor0966/Documents/code/COSMOS-analysis/data/castanetv9.3_output/COSMOS1-4_castanetv9.3_nc2.post_filt_false.tsv"
    reads_dedup_post_filt_false = "/Users/afor0966/Documents/code/COSMOS-analysis/data/castanetv9.3_output/COSMOS1-4_castanetv9.3_reads_dedup.post_filt_false.tsv"
    config_file = "/Users/afor0966/Documents/code/COSMOS-analysis/config/COSMOS1-4.local.json"
    metadata = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"

    with open(config_file, 'r') as f:
        config = json.load(f)

    nc2_post_filt_true_df = pl.read_csv(nc2_post_filt_true, separator="\t")
    nc2_post_filt_false_df = pl.read_csv(nc2_post_filt_false, separator="\t")
    reads_dedup_post_filt_true_df = pl.read_csv(reads_dedup_post_filt_true, separator="\t")
    reads_dedup_post_filt_false_df = pl.read_csv(reads_dedup_post_filt_false, separator="\t")
    metadata_df = pl.read_csv(metadata, separator="\t")

    # Include every metadata sample and plot missing counts as zeros.
    metadata_samples = metadata_df.select(pl.col("sample").alias("sampleid")).unique()

    def include_missing_samples(df):
        return metadata_samples.join(df, on="sampleid", how="left").with_columns(
            pl.col(col).cast(pl.Float64).fill_null(0)
            for col in df.columns if col != "sampleid"
        )

    nc2_post_filt_true_df = include_missing_samples(nc2_post_filt_true_df)
    nc2_post_filt_false_df = include_missing_samples(nc2_post_filt_false_df)
    reads_dedup_post_filt_true_df = include_missing_samples(reads_dedup_post_filt_true_df)
    reads_dedup_post_filt_false_df = include_missing_samples(reads_dedup_post_filt_false_df)

    # make sample to sampleset dict
    sample_to_sampleset = {}
    for sample_set in config['sample_sets']:
        for sample in config['sample_sets'][sample_set]['samples']:
            sample_to_sampleset[sample] = sample_set

    # add sample_set column to nc2 and reads_dedup dataframes
    nc2_post_filt_true_df = nc2_post_filt_true_df.with_columns(
        pl.col("sampleid").replace_strict(sample_to_sampleset).alias("sample_set")
    )
    nc2_post_filt_false_df = nc2_post_filt_false_df.with_columns(
        pl.col("sampleid").replace_strict(sample_to_sampleset).alias("sample_set")
    )
    reads_dedup_post_filt_true_df = reads_dedup_post_filt_true_df.with_columns(
        pl.col("sampleid").replace_strict(sample_to_sampleset).alias("sample_set")
    )
    reads_dedup_post_filt_false_df = reads_dedup_post_filt_false_df.with_columns(
        pl.col("sampleid").replace_strict(sample_to_sampleset).alias("sample_set")
    )

    # add postfilt true to column headers
    nc2_post_filt_true_df = nc2_post_filt_true_df.rename({col: f"{col}_post_filt_true" for col in nc2_post_filt_true_df.columns if col != "sampleid"})
    reads_dedup_post_filt_true_df = reads_dedup_post_filt_true_df.rename({col: f"{col}_post_filt_true" for col in reads_dedup_post_filt_true_df.columns if col != "sampleid"})
    nc2_post_filt_false_df = nc2_post_filt_false_df.rename({col: f"{col}_post_filt_false" for col in nc2_post_filt_false_df.columns if col != "sampleid"})
    reads_dedup_post_filt_false_df = reads_dedup_post_filt_false_df.rename({col: f"{col}_post_filt_false" for col in reads_dedup_post_filt_false_df.columns if col != "sampleid"})

    # log10(reads + 1) keeps zero-read samples at zero; leave sample labels alone.
    reads_dedup_post_filt_true_df = reads_dedup_post_filt_true_df.with_columns(
        (pl.col(pl.Float64) + 1).log10()
    )
    reads_dedup_post_filt_false_df = reads_dedup_post_filt_false_df.with_columns(
        (pl.col(pl.Float64) + 1).log10()
    )

    # combine nc2 and reads_dedup dataframes
    combined_nc2_df = nc2_post_filt_true_df.join(nc2_post_filt_false_df, on="sampleid", how="outer")
    combined_reads_dedup_df = reads_dedup_post_filt_true_df.join(reads_dedup_post_filt_false_df, on="sampleid", how="outer")

    # convert all columns to float except for sampleid
    combined_nc2_df = combined_nc2_df.with_columns([pl.col(col).cast(pl.Float64) for col in combined_nc2_df.columns if "sample" not in col])
    combined_reads_dedup_df = combined_reads_dedup_df.with_columns([pl.col(col).cast(pl.Float64) for col in combined_reads_dedup_df.columns if "sample" not in col])

    # print their shapes
    print(f"Combined NC2 dataframe shape: {combined_nc2_df.shape}")
    print(f"Combined Reads Dedup dataframe shape: {combined_reads_dedup_df.shape}")

    # Report samples where BFV NC2 crosses the 0.05 threshold after filtering.
    lower_bfv_samples = combined_nc2_df.filter(
        (pl.col("barmahforestvirus_post_filt_false") > 0.05)
        & (pl.col("barmahforestvirus_post_filt_true") < 0.05)
    ).select("sampleid").sort("sampleid")
    print(
        "Samples with BFV NC2 >0.05 before and <0.05 after filtering: "
        f"{lower_bfv_samples.height}"
    )
    print(lower_bfv_samples["sampleid"].to_list())
    lower_bfv_by_sample_set = (
        lower_bfv_samples
        .join(
            combined_nc2_df.select(["sampleid", "sample_set_post_filt_true"]),
            on="sampleid",
            how="left",
        )
        .group_by("sample_set_post_filt_true")
        .agg(pl.len().alias("sample_count"))
        .sort("sample_set_post_filt_true")
        .rename({"sample_set_post_filt_true": "sample_set"})
    )
    print("Lower BFV samples by sample set:")
    print(lower_bfv_by_sample_set)

    return combined_nc2_df, combined_reads_dedup_df

def make_scatters(combined_nc2_df, combined_reads_dedup_df, output_dir, targets_to_plot={"mosquitocox1": "COX1"}):
    point_size = 80
    # make scatter plots
    sns.set(style="whitegrid")
    sample_sets = sorted(set(
        combined_nc2_df["sample_set_post_filt_true"].drop_nulls().to_list()
        + combined_reads_dedup_df["sample_set_post_filt_true"].drop_nulls().to_list()
    ))
    palette = dict(zip(sample_sets, sns.color_palette(n_colors=len(sample_sets))))

    # Give both legends their own rows below the plots, allowing for wrapping.
    sample_columns = len(sample_sets)
    target_columns = len(targets_to_plot.keys())
    sample_height = 0.35 * (1 + (len(sample_sets) + sample_columns - 1) // sample_columns)
    target_height = 0.35 * (1 + (len(targets_to_plot) + target_columns - 1) // target_columns)
    fig = plt.figure(figsize=(10, 5 + sample_height + target_height), layout="constrained")
    grid = fig.add_gridspec(3, 2, height_ratios=[5, sample_height, target_height])
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    sample_legend_ax = fig.add_subplot(grid[1, :])
    target_legend_ax = fig.add_subplot(grid[2, :])
    sample_legend_ax.set_axis_off()
    target_legend_ax.set_axis_off()

    markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "*", "h", "H", "8"]

    # nc2 scatter plot
    for i, target in enumerate(targets_to_plot.keys()):
        sns.scatterplot(data=combined_nc2_df.to_pandas(), x=f"{target}_post_filt_false", y=f"{target}_post_filt_true", ax=axes[0], 
                        hue="sample_set_post_filt_true", hue_order=sample_sets, palette=palette, marker=markers[i % len(markers)],
                        legend=False, s=point_size)
    axes[0].set_title("Genome coverage (propn covered ≥ 2x)", fontsize='large')
    axes[0].set_xlabel("PostFilt=False")
    axes[0].set_ylabel("PostFilt=True")

    # reads_dedup scatter plot
    for i, target in enumerate(targets_to_plot.keys()):
        sns.scatterplot(data=combined_reads_dedup_df.to_pandas(), x=f"{target}_post_filt_false", y=f"{target}_post_filt_true", ax=axes[1], 
                        hue="sample_set_post_filt_true", hue_order=sample_sets, palette=palette, marker=markers[i % len(markers)], 
                        legend=False, s=point_size)
    axes[1].set_title("Log10(unique reads + 1)", fontsize='large')
    axes[1].set_xlabel("PostFilt=False")
    axes[1].set_ylabel("PostFilt=True")

    # Match x/y limits with equal padding and a full-length equality line.
    for ax in axes:
        upper = max(ax.dataLim.xmax, ax.dataLim.ymax, 1)
        padding = 0.05 * upper
        limits = (-padding, upper + padding)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.axline((0, 0), slope=1, color="black", linestyle="--", linewidth=1, zorder=0)
        sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("black")
        sns.despine(ax=ax, top=True, right=True, left=False, bottom=False, trim=True)

    # add shared legend
    # a coloured patch for each sample set
    handles = [plt.Line2D([0], [0], marker='o', linestyle='None', color=palette[sample_set], label=sample_set, markersize=10) for sample_set in sample_sets]
    sample_legend_ax.legend(handles=handles, title="Sample set", loc="center", ncol=sample_columns, frameon=False)
    # shape for each target
    handles = [plt.Line2D([0], [0], marker=markers[i % len(markers)], linestyle='None', color='0.3', label=targets_to_plot[target], markersize=10) for i, target in enumerate(targets_to_plot)]
    target_legend_ax.legend(handles=handles, title="Target", loc="center", ncol=target_columns, frameon=False)

    fig.savefig(os.path.join(output_dir, "index_hopping_by_batch.png"), bbox_inches="tight")
    plt.close()

def main():
    virus_acronyms = {
        "barmahforestvirus": "BFV",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "kokoberavirus": "KOK",
        "murrayvalleyencephalitisvirus": "MVEV",
        "rossrivervirus": "RRV",
        "stratfordvirus": "STRV",
    }
    output_dir = "results/publication_figures/plots"
    nc2_df, reads_dedup_df = parse_data()

    # print max nc2 and dedup read count for a sample where BFV is <0.05 or <100 in PostFilt is True
    bfv_nc2_df = nc2_df.filter(pl.col("barmahforestvirus_post_filt_true") < 0.05)
    bfv_reads_dedup_df = reads_dedup_df.filter(pl.col("barmahforestvirus_post_filt_true") < 100)
    max_bfv_nc2 = bfv_nc2_df.select(pl.max("barmahforestvirus_post_filt_false")).to_numpy()[0][0]
    max_bfv_reads_dedup = bfv_reads_dedup_df.select(pl.max("barmahforestvirus_post_filt_false")).to_numpy()[0][0]
    print('number of samples with BFV < 0.05 in PostFilt=True:', bfv_nc2_df['sampleid'].n_unique())
    print('number of samples with BFV < 100 in PostFilt=True:', bfv_reads_dedup_df['sampleid'].n_unique())

    make_scatters(nc2_df, reads_dedup_df, output_dir, targets_to_plot=virus_acronyms)

if __name__ == "__main__":
    main()
