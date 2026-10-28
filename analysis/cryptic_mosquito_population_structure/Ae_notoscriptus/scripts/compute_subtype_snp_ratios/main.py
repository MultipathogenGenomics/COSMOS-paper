from analysis.Cx_annulirostris.compute_subtype_snp_ratios.process import compute_snp_ratios_for_subtype_across_samples, make_results_readable
import os

def main():
    dataset = "COSMOS4"
    snp_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_2_or_3.tsv"
    bam_dir = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/{dataset}/align_kma_assigned/clustered_0.995_pid_trimmed/Ae_notoscriptus"
    target_id = "KM676218.1_Aedes_notoscriptus_isolate_MC04_mitochondrion_complete_genome"
    tabulate_alignments_executable = "/Users/afor0966/Documents/code/software/iorek/build/apps/tabulate-alignments"
    output_path = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples.1_vs_2or3.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # analysis
    results_df = compute_snp_ratios_for_subtype_across_samples(snp_file, bam_dir, target_id, tabulate_alignments_executable)
    results_df.write_csv(output_path, separator='\t')
    # make readable version of results
    results_df = make_results_readable(results_df, snp_file)
    readable_output_path = output_path.replace(".tsv", ".readable.tsv")
    cols = ["sample", "position_in_reference", "total_counts_at_position", "1_subtype_bases", "1_subtype_count", "1_non_subtype_bases", "1_non_subtype_count", "1_ratio"]
    results_df = results_df.select(cols)
    results_df.write_csv(readable_output_path, separator='\t')

    # run same but 2 vs 3
    snp_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_2_vs_3_fixed.tsv"
    output_path = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples_2_vs_3.tsv"
    # analysis
    results_df = compute_snp_ratios_for_subtype_across_samples(snp_file, bam_dir, target_id, tabulate_alignments_executable)
    results_df.write_csv(output_path, separator='\t')
    # make readable version of results
    results_df = make_results_readable(results_df, snp_file)
    readable_output_path = output_path.replace(".tsv", ".readable.tsv")
    cols = ["sample", "position_in_reference", "total_counts_at_position", "2_subtype_bases", "2_subtype_count", "2_non_subtype_bases", "2_non_subtype_count", "2_ratio"]
    results_df = results_df.select(cols)
    results_df.write_csv(readable_output_path, separator='\t')

if __name__ == "__main__":
    main()