from analysis.Cx_annulirostris.collate_results.process import collate_count_and_linkage_results_by_sample
from analysis.Ae_notoscriptus.collate_results.process import get_Ae_notoscriptus_cluster1_to_3_ratios, combine_chunks_of_linkage_results
import os

def main():
    dataset = "COSMOS3"
    counts_results_1_vs_2or3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples.1_vs_2or3.readable.tsv"
    counts_results_2_vs_3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples_2_vs_3.readable.tsv"
    linkage_results_1_vs_2or3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/linkage_check/{dataset}.Ae_notoscriptus_linkage_check_results.1_vs_2or3.tsv"
    linkage_results_2_vs_3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/linkage_check/{dataset}.Ae_notoscriptus_linkage_check_results_2_vs_3.tsv"
    output_file_1_vs_2or3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/{dataset}.Ae_notoscriptus_collated_results.1_vs_2or3.tsv"
    output_file_2_vs_3 = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/{dataset}.Ae_notoscriptus_collated_results.2_vs_3.tsv"
    output_all_clusters = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/{dataset}.Ae_notoscriptus_collated_results.subtype_ratios_clusters1-3.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_file_1_vs_2or3), exist_ok=True)
    os.makedirs(os.path.dirname(output_file_2_vs_3), exist_ok=True)
    os.makedirs(os.path.dirname(output_all_clusters), exist_ok=True)
    # collate results
    # combine linkage chunks
    chunk_file_template = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/linkage_check/{dataset}.Ae_notoscriptus_linkage_check_results_2_vs_3.chunk="
    chunk_files = []
    for chunk in range(0, 10):
        chunk_files.append(f"{chunk_file_template}{chunk}.tsv")
    combined_chunks = combine_chunks_of_linkage_results(chunk_files)
    combined_chunks.write_csv(linkage_results_2_vs_3, separator='\t')
    chunk_files = [ i.replace("_2_vs_3", ".1_vs_2or3") for i in chunk_files ]
    combined_chunks = combine_chunks_of_linkage_results(chunk_files)
    combined_chunks.write_csv(linkage_results_1_vs_2or3, separator='\t')
    # I think do for 1 vs 2 or 3, then 2 vs 3 and compute 2/3 as proportion of against 1
    collate_count_and_linkage_results_by_sample(counts_results_1_vs_2or3, linkage_results_1_vs_2or3, output_file_1_vs_2or3, count_threshold_for_averages=20)
    collate_count_and_linkage_results_by_sample(counts_results_2_vs_3, linkage_results_2_vs_3, output_file_2_vs_3, count_threshold_for_averages=20)
    # now get ratio for each cluster
    get_Ae_notoscriptus_cluster1_to_3_ratios(output_file_1_vs_2or3, output_file_2_vs_3, output_all_clusters)

if __name__ == "__main__":
    main()