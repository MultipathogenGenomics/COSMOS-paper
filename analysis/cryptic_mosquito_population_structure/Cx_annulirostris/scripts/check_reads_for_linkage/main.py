from analysis.common.paths import RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.check_reads_for_linkage.process import check_if_reads_support_subtypes_by_linkage
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.count_subtype_kmers.process import get_list_of_annulirostris_samples
import os

def main():
    dataset = "COSMOS4"
    snp_ratio_file = f"{RESULTS_DIR}/Cx_annulirostris/subtype_snp_ratios/{dataset}.Cx_annulirostris_subtype_defining_snp_ratios_across_samples.readable.tsv"
    bam_dir = f"{RESULTS_DIR}/{dataset}/align_kma_assigned/clustered_0.995_pid_trimmed/Cx_annulirostris"
    output_path = f"{RESULTS_DIR}/Cx_annulirostris/linkage_check/{dataset}.Cx_annulirostris_linkage_check_results.tsv"
    focal_subtype = "S-AUS"
    contig = "mosquito_NC_054323.1-1430-2966-Culex-annulirostris-isolate-W01-mitochondrion"
    abundance_results_file = f"{RESULTS_DIR}/COSMOS1-4/species_id_and_abundance/clustered_0.995_pid_trimmed/COSMOS1-4.broad_dataset.basic_filtering.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # analysis
    sample_list = get_list_of_annulirostris_samples(abundance_results_file)
    results = check_if_reads_support_subtypes_by_linkage(snp_ratio_file, bam_dir, abundance_results_file, focal_subtype, contig, mapq_threshold=20, sample_list=sample_list)
    results.write_csv(output_path, separator='\t')

if __name__ == "__main__":
    main()