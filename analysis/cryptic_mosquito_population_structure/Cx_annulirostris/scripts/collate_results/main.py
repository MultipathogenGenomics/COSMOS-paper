from analysis.common.paths import RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.collate_results.process import collate_count_and_linkage_results_by_sample
import os

def main():
    dataset = "COSMOS4"
    counts_results = f"{RESULTS_DIR}/Cx_annulirostris/subtype_snp_ratios/{dataset}.Cx_annulirostris_subtype_defining_snp_ratios_across_samples.readable.tsv"
    linkage_results = f"{RESULTS_DIR}/Cx_annulirostris/linkage_check/{dataset}.Cx_annulirostris_linkage_check_results.tsv"
    output_file = f"{RESULTS_DIR}/Cx_annulirostris/collate_results/{dataset}.Cx_annulirostris_collated_results.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # collate results
    collate_count_and_linkage_results_by_sample(counts_results, linkage_results, output_file, count_threshold_for_averages=20)

if __name__ == "__main__":
    main()