import polars as pl

def parse_tsv(file_path: str) -> pl.DataFrame:
    return pl.read_csv(file_path, separator='\t')

def get_list_of_notoscriptus_samples(abundance_results_file: str) -> list[str]:
    abundance_df = parse_tsv(abundance_results_file)
    samples = abundance_df.filter(pl.col("kma:read_count") > 50).filter(pl.col("target_id").str.contains("notoscriptus"))["sample"].unique().to_list()
    return samples