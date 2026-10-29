from analysis.common.paths import RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.count_subtype_kmers.process import count_subtype_kmers_across_samples
import polars as pl

def main():
    dataset = "COSMOS3"
    kmer_file = f"{RESULTS_DIR}/Cx_annulirostris/Cx_annulirostris.subtype_private_kmers.min_fraction_0.5.tsv"
    counts_dir = f"{RESULTS_DIR}/{dataset}/kmer_counts/Cx_annulirostris"
    abundance_results_file = f"{RESULTS_DIR}/COSMOS1-4/species_id_and_abundance/clustered_0.995_pid_trimmed/COSMOS1-4.broad_dataset.basic_filtering.tsv"
    subtypes_to_test = ["ann AUS", "ann S-AUS", "ann SI", "ann PNG-1", "ann PNG-2"]
    # subtypes_to_test = ["ann S-AUS", "ann AUS"]
    # run analysis
    results_df = count_subtype_kmers_across_samples(kmer_file, counts_dir, abundance_results_file, subtypes_to_test)
    print(results_df['count'].describe())
    # print samples where count > 0
    print(results_df.select(["sample", "kmer", "subtype", "count"]))
    results_df.write_csv(f"{RESULTS_DIR}/Cx_annulirostris/subtype_kmer_counts/{dataset}_Cx_annulirostris_subtype_kmer_counts.min_fraction_0.5.csv")

if __name__ == "__main__":
    main()