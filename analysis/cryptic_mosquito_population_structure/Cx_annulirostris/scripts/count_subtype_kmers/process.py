
import polars as pl
import os

def parse_tsv(file_path: str) -> pl.DataFrame:
    return pl.read_csv(file_path, separator='\t')

def count_kmers(kmer_file: str, counts_file: str, test_subtype: str) -> pl.DataFrame:
    # parse kmer file
    kmer_df = parse_tsv(kmer_file)
    # filter to subtype
    kmer_df = kmer_df.filter(pl.col("subtype") == test_subtype)
    # parse counts file
    counts_df = parse_tsv(counts_file)
    # merge on kmer
    merged_df = kmer_df.join(counts_df, on="kmer", how="left")
    # fill missing counts with 0
    merged_df = merged_df.fill_null(0)
    return merged_df

def get_list_of_annulirostris_samples(abundance_results_file: str) -> list[str]:
    abundance_df = parse_tsv(abundance_results_file)
    samples = abundance_df.filter(pl.col("kma:read_count") > 50).filter(pl.col("target_id").str.contains("annulirostris"))["sample"].unique().to_list()
    return samples

def count_subtype_kmers_across_samples(kmer_file: str, counts_dir: str, abundance_results_file: str, subtypes_to_test: list[str]) -> pl.DataFrame:
    all_results = []
    samples_to_include = get_list_of_annulirostris_samples(abundance_results_file)
    samples = os.listdir(counts_dir)
    for subtype in subtypes_to_test:
        for sample in samples:
            if sample not in samples_to_include:
                continue
            counts_file = os.path.join(counts_dir, sample, f"{sample}.counts_as_table.tsv")
            result_df = count_kmers(kmer_file, counts_file, subtype)
            result_df = result_df.with_columns(pl.lit(sample).alias("sample"))
            # make sure count is int
            result_df = result_df.with_columns(pl.col("count").cast(pl.Int64))
            all_results.append(result_df)
    return pl.concat(all_results)
    