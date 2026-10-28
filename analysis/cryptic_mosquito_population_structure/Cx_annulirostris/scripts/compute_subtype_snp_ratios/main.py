from analysis.Cx_annulirostris.compute_subtype_snp_ratios.process import compute_snp_ratios_for_subtype_across_samples, make_results_readable
import os
import polars as pl


def summarise_readable_results(readable_df: pl.DataFrame, dataset: str, comparison: str) -> pl.DataFrame:
    ratio_cols = [col for col in readable_df.columns if col.endswith("_ratio")]
    summary_exprs = [
        pl.len().alias("n_positions"),
        pl.col("total_counts_at_position").mean().alias("mean_total_depth"),
        pl.col("total_counts_at_position").min().alias("min_total_depth"),
        pl.col("total_counts_at_position").max().alias("max_total_depth"),
    ]
    summary_exprs.extend(
        pl.col(col).mean().alias(col.replace("_ratio", "_avg_ratio"))
        for col in ratio_cols
    )
    summary_df = (
        readable_df
        .group_by("sample")
        .agg(summary_exprs)
        .with_columns(
            pl.lit(dataset).alias("dataset"),
            pl.lit(comparison).alias("comparison"),
        )
    )
    return summary_df.select(
        ["dataset", "comparison", "sample"] +
        [col for col in summary_df.columns if col not in ["dataset", "comparison", "sample"]]
    )


def main():
    datasets = ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]
    snp_dir = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris"
    output_dir = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/subtype_snp_ratios"
    comparisons = {
        "AUS_vs_S-AUS": "Cx_annulirostris.SNPs_between_AUS_and_S-AUS.tsv",
        "AUS_S-AUS_PNG-1": "Cx_annulirostris.SNPs_between_AUS_S-AUS_and_PNG-1.tsv",
        "AUS_S-AUS_PNG-2": "Cx_annulirostris.SNPs_between_AUS_S-AUS_and_PNG-2.tsv",
        "AUS_S-AUS_SI": "Cx_annulirostris.SNPs_between_AUS_S-AUS_and_SI.tsv",
        "all_subtypes": "Cx_annulirostris.SNPs_between_subtype_and_other_fragments.tsv",
    }
    target_id = "mosquito_NC_054323.1-1430-2966-Culex-annulirostris-isolate-W01-mitochondrion"
    tabulate_alignments_executable = "/Users/afor0966/Documents/code/software/iorek/build/apps/tabulate-alignments"

    os.makedirs(output_dir, exist_ok=True)

    for dataset in datasets:
        bam_dir = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/{dataset}/align_kma_assigned/clustered_0.995_pid_trimmed/Cx_annulirostris"
        for comparison, snp_filename in comparisons.items():
            snp_file = os.path.join(snp_dir, snp_filename)
            output_path = os.path.join(
                output_dir,
                f"{dataset}.{comparison}.Cx_annulirostris_subtype_defining_snp_ratios_across_samples.tsv"
            )
            readable_output_path = output_path.replace(".tsv", ".readable.tsv")
            summary_output_path = output_path.replace(".tsv", ".summary.tsv")

            print(f"Running {dataset} {comparison}")
            results_df = compute_snp_ratios_for_subtype_across_samples(
                snp_file, bam_dir, target_id, tabulate_alignments_executable
            )
            results_df.write_csv(output_path, separator='\t')

            readable_df = make_results_readable(results_df, snp_file)
            readable_df.write_csv(readable_output_path, separator='\t')

            summary_df = summarise_readable_results(readable_df, dataset, comparison)
            summary_df.write_csv(summary_output_path, separator='\t')

if __name__ == "__main__":
    main()
