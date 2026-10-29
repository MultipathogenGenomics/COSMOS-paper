# kma and kallisto index and quantification functions for in silico mixtures

import os
import subprocess
import polars as pl
from analysis.common.mosquito_ids import return_species_from_target_id

def build_kallisto_index(fasta_file: str, index_prefix: str) -> None:
    command = f"kallisto index -i {index_prefix} {fasta_file}"
    subprocess.run(command, shell=True, check=True)

def build_kma_index(fasta_file: str, index_prefix: str) -> None:
    command = f"kma index -i {fasta_file} -o {index_prefix}"
    subprocess.run(command, shell=True, check=True)

def run_kallisto_quantification(index_prefix: str, read1: str, read2: str, output_dir: str) -> None:
    command = f"kallisto quant -i {index_prefix} -o {output_dir} -t 4 {read1} {read2}"
    subprocess.run(command, shell=True, check=True)

def run_kma_quantification(index_prefix: str, read1: str, read2: str, output_prefix: str, keep_files: bool = False) -> None:
    command = f"kma -ipe {read1} {read2} -o {output_prefix} -t_db {index_prefix} -t 4 -ef -ex_mode -vcf 2 -1t1 > /dev/null" # check this matches
    subprocess.run(command, shell=True, check=True)
    # delete the frag file to save space
    if keep_files:
        if os.path.exists(f"{output_prefix}.frag.gz"):
            os.remove(f"{output_prefix}.frag.gz")
        # and vcf file
        if os.path.exists(f"{output_prefix}.vcf.gz"):
            os.remove(f"{output_prefix}.vcf.gz")

def parse_kallisto_output(abundance_file: str) -> pl.DataFrame:
    df = pl.read_csv(abundance_file, separator="\t")
    df = df.select([
        pl.col("target_id").cast(pl.String).alias("reference"),
        pl.col("est_counts").cast(pl.Float64).alias("kallisto_estimated_counts"),
        pl.col("tpm").cast(pl.Float64).alias("kallisto_tpm"),
    ])
    # convert reference to species name
    df = df.with_columns(pl.col("reference").map_elements(lambda x: return_species_from_target_id(x).replace(" ", "_")).alias("species"))
    return df

def parse_kma_output(kma_output_file: str) -> pl.DataFrame:
    df = pl.read_csv(kma_output_file, separator="\t", comment_prefix="##")
    df = df.select([
        pl.col("# refSequence").cast(pl.String).alias("reference"),
        pl.col("readCount").cast(pl.Float64).alias("kma_read_count"),
        pl.col("fragmentCount").cast(pl.Float64).alias("kma_fragment_count"),
    ])
    # convert reference to species name
    df = df.with_columns(pl.col("reference").map_elements(lambda x: return_species_from_target_id(x).replace(" ", "_")).alias("species"))
    # get the query and template coverage from the res file
    res_df = pl.read_csv(kma_output_file.replace(".mapstat", ".res"), separator="\t")
    assert len(res_df) == len(df), "Number of rows in res file does not match number of rows in mapstat file"
    res_df = res_df.select([
        pl.col("#Template").cast(pl.String).alias("reference"),
        pl.col("Template_Identity").str.strip_chars().cast(pl.Float64).alias("kma_template_identity"),
        pl.col("Query_Coverage").str.strip_chars().cast(pl.Float64).alias("kma_query_coverage"),
        pl.col("Template_Coverage").str.strip_chars().cast(pl.Float64).alias("kma_template_coverage"),
        pl.col("Query_Identity").str.strip_chars().cast(pl.Float64).alias("kma_query_identity"),
    ])
    df = df.join(res_df, on="reference", how="inner")
    return df

def run_kallisto_and_kma_quantification(index_file: str, fastq_dir: str, output_dir: str) -> None:
    mixture_metadata = pl.read_csv(f"{fastq_dir}/mixture_metadata.tsv", separator='\t')
    # make indexes if they don't exist
    if not os.path.exists(f"{output_dir}/kallisto_index.idx"):
        build_kallisto_index(index_file, f"{output_dir}/kallisto_index.idx")
    if not os.path.exists(f"{output_dir}/kma_index.idx"):
        build_kma_index(index_file, f"{output_dir}/kma_index.idx")
    results_dfs = []
    for row in mixture_metadata.iter_rows(named=True):
        mixture_id = row['mixture_id']
        read1 = f"{fastq_dir}/mixture_{mixture_id}_R1.fastq.gz"
        read2 = f"{fastq_dir}/mixture_{mixture_id}_R2.fastq.gz"
        kallisto_output_dir = f"{output_dir}/kallisto_mixture_{mixture_id}"
        kma_output_prefix = f"{output_dir}/kma_mixture_{mixture_id}"
        run_kallisto_quantification(f"{output_dir}/kallisto_index.idx", read1, read2, kallisto_output_dir)
        run_kma_quantification(f"{output_dir}/kma_index.idx", read1, read2, kma_output_prefix)
        kallisto_df = parse_kallisto_output(f"{kallisto_output_dir}/abundance.tsv")
        kma_df = parse_kma_output(f"{kma_output_prefix}.mapstat")
        kallisto_df = kallisto_df.with_columns(pl.lit(mixture_id).alias("mixture_id"))
        kma_df = kma_df.with_columns(pl.lit(mixture_id).alias("mixture_id"))
        merged_df = kallisto_df.join(kma_df, on=["reference", "species", "mixture_id"], how="inner")
        results_dfs.append(merged_df)
    final_df = pl.concat(results_dfs)
    final_df.write_csv(f"{output_dir}/kallisto_kma_quantification_results.tsv", separator='\t')