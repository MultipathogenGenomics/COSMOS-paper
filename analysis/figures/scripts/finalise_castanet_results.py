from analysis.common.paths import DATA_DIR, RESULTS_DIR
import json

import polars as pl
import os
from datetime import date
import json

# make data directory
output_dir = f"{RESULTS_DIR}/final_castanet_data"
os.makedirs(output_dir, exist_ok=True)

# add a read me saying this data is post_filt = False for all bar COSMOS2, which we use post_filt is True due to index hopping
readme_file = os.path.join(output_dir, "README.txt")
with open(readme_file, 'w') as f:
    f.write(
        "Final Castanet data\n"
        f"Date generated: {date.today().isoformat()}\n\n"
        "These data use post_filt = False for all COSMOS batches except COSMOS2.\n"
        "COSMOS2 uses post_filt = True due to index hopping.\n"
    )

# parse files
nc2_post_filt_true = f"{DATA_DIR}/castanetv9.2-output/COSMOS1-4_castanetv9.2_nc2.tsv"
reads_dedup_post_filt_true = f"{DATA_DIR}/castanetv9.2-output/COSMOS1-4_castanetv9.2_reads_dedup.tsv"
nc2_post_filt_false = f"{DATA_DIR}/castanetv9.3_output/COSMOS1-4_castanetv9.3_nc2.post_filt_false.tsv"
reads_dedup_post_filt_false = f"{DATA_DIR}/castanetv9.3_output/COSMOS1-4_castanetv9.3_reads_dedup.post_filt_false.tsv"
config_file = "config/COSMOS1-4.json"

with open(config_file, 'r') as f:
    config = json.load(f)

nc2_post_filt_true_df = pl.read_csv(nc2_post_filt_true, separator="\t")
nc2_post_filt_false_df = pl.read_csv(nc2_post_filt_false, separator="\t")
reads_dedup_post_filt_true_df = pl.read_csv(reads_dedup_post_filt_true, separator="\t")
reads_dedup_post_filt_false_df = pl.read_csv(reads_dedup_post_filt_false, separator="\t")

original_shape = nc2_post_filt_true_df.shape[0] 

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

# filter the post filt false to only have COSMOS1, 3 and 4, and filter the post filt true to only have COSMOS2
nc2_post_filt_false_df = nc2_post_filt_false_df.filter(pl.col("sample_set") != "COSMOS2")
reads_dedup_post_filt_false_df = reads_dedup_post_filt_false_df.filter(pl.col("sample_set") != "COSMOS2")
nc2_post_filt_true_df = nc2_post_filt_true_df.filter(pl.col("sample_set") == "COSMOS2")
reads_dedup_post_filt_true_df = reads_dedup_post_filt_true_df.filter(pl.col("sample_set") == "COSMOS2")

# concat the nc2 and reads_dedup dataframes
final_nc2_df = pl.concat([nc2_post_filt_false_df, nc2_post_filt_true_df], how="diagonal_relaxed")
final_reads_dedup_df = pl.concat([reads_dedup_post_filt_false_df, reads_dedup_post_filt_true_df], how="diagonal_relaxed")

# remove the sample_set column from the final dataframes
final_nc2_df = final_nc2_df.drop("sample_set")
final_reads_dedup_df = final_reads_dedup_df.drop("sample_set")

# assert that the shapes are the same
assert final_nc2_df.shape[0] == final_reads_dedup_df.shape[0], "The number of rows in the final nc2 and reads_dedup dataframes are not the same."

# save to output directory
final_nc2_df.write_csv(os.path.join(output_dir, "COSMOS1-4.final_castanet_nc2.tsv"), separator="\t")
final_reads_dedup_df.write_csv(os.path.join(output_dir, "COSMOS1-4.final_castanet_reads_dedup.tsv"), separator="\t")

# now the same for the virus coverage files
virus_coverage_post_filt_false = f"{RESULTS_DIR}/virus/data/virus_coverage_per_sample.post_filt_false.tsv"
virus_coverage_post_filt_true = f"{RESULTS_DIR}/virus/data/virus_coverage_per_sample.post_filt_true.tsv"
virus_coverage_post_filt_false_df = pl.read_csv(virus_coverage_post_filt_false, separator="\t")
virus_coverage_post_filt_true_df = pl.read_csv(virus_coverage_post_filt_true, separator="\t")

# add sample_set column to virus coverage dataframes
virus_coverage_post_filt_false_df = virus_coverage_post_filt_false_df.with_columns(
    pl.col("sample").replace_strict(sample_to_sampleset).alias("sample_set")
)
virus_coverage_post_filt_true_df = virus_coverage_post_filt_true_df.with_columns(
    pl.col("sample").replace_strict(sample_to_sampleset).alias("sample_set")
)

# filter the post filt false to only have COSMOS1, 3 and 4, and filter the post filt true to only have COSMOS2
virus_coverage_post_filt_false_df = virus_coverage_post_filt_false_df.filter(pl.col("sample_set") != "COSMOS2")
virus_coverage_post_filt_true_df = virus_coverage_post_filt_true_df.filter(pl.col("sample_set") == "COSMOS2")

# concat the virus coverage dataframes
final_virus_coverage_df = pl.concat([virus_coverage_post_filt_false_df, virus_coverage_post_filt_true_df], how="diagonal_relaxed")

# remove the sample_set column from the final virus coverage dataframe
final_virus_coverage_df = final_virus_coverage_df.drop("sample_set")

# save to output directory
final_virus_coverage_df.write_csv(os.path.join(output_dir, "COSMOS1-4.final_castanet_virus_coverage_per_sample.tsv"), separator="\t")