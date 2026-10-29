from analysis.common.paths import RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.check_reads_for_linkage.process import check_if_reads_support_subtypes_by_linkage
from analysis.cryptic_mosquito_population_structure.Ae_notoscriptus.scripts.check_reads_for_linkage.process import get_list_of_notoscriptus_samples
import os, sys

def main():
    chunk = int(sys.argv[1]) if len(sys.argv) > 1 else None
    dataset = "COSMOS3"
    snp_ratio_file = f"{RESULTS_DIR}/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples.1_vs_2or3.readable.tsv"
    snp_ratio_file_2_vs_3 = f"{RESULTS_DIR}/Ae_notoscriptus/subtype_snp_ratios/{dataset}.Ae_notoscriptus_subtype_defining_snp_ratios_across_samples_2_vs_3.readable.tsv"
    bam_dir = f"{RESULTS_DIR}/{dataset}/align_kma_assigned/clustered_0.995_pid_trimmed/Ae_notoscriptus"
    output_path = f"{RESULTS_DIR}/Ae_notoscriptus/linkage_check/{dataset}.Ae_notoscriptus_linkage_check_results.1_vs_2or3.chunk={chunk}.tsv"
    output_path_2_vs_3 = f"{RESULTS_DIR}/Ae_notoscriptus/linkage_check/{dataset}.Ae_notoscriptus_linkage_check_results_2_vs_3.chunk={chunk}.tsv"
    focal_subtype = "1"
    contig = "KM676218.1_Aedes_notoscriptus_isolate_MC04_mitochondrion_complete_genome"
    abundance_results_file = f"{RESULTS_DIR}/COSMOS1-4/species_id_and_abundance/clustered_0.995_pid_trimmed/COSMOS1-4.broad_dataset.basic_filtering.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # analysis
    sample_list = get_list_of_notoscriptus_samples(abundance_results_file)
    results = check_if_reads_support_subtypes_by_linkage(snp_ratio_file, bam_dir, abundance_results_file, focal_subtype, contig, 
                                                         mapq_threshold=20, sample_list=sample_list, chunk=chunk, total_chunks=10)
    results.write_csv(output_path, separator='\t')
    results_2_vs_3 = check_if_reads_support_subtypes_by_linkage(snp_ratio_file_2_vs_3, bam_dir, abundance_results_file, focal_subtype="2", contig=contig, 
                                                                mapq_threshold=20, sample_list=sample_list, chunk=chunk, total_chunks=10)
    results_2_vs_3.write_csv(output_path_2_vs_3, separator='\t')

if __name__ == "__main__":
    main()