# mixing sequence reads in silico functions to test kma and kallisto abundance estimation

import subprocess
from itertools import combinations
from analysis.common.utils import count_reads_in_fastq
from analysis.common.io import parse_json
import numpy as np
import polars as pl
import json

def make_sets_n_choose_k(items, k):
    return list(combinations(items, k))

def get_number_of_reads_to_sample(r1_files: list, proportions: list) -> tuple:
    total_reads_by_sample = np.array([count_reads_in_fastq(f) for f in r1_files])
    props = np.array(proportions)
    total_reads = np.floor(total_reads_by_sample / props).min()
    total_reads_to_sample = total_reads * props
    return list(total_reads_to_sample.astype(int))

def subsample_fastq(input_fastq: str, output_fastq: str, n_reads: int) -> None:
    command = f"seqtk sample -s100 {input_fastq} {n_reads} >> {output_fastq}"
    subprocess.run(command, shell=True, check=True)

def gzip_fastq(input_fastq: str, keep_original: bool = False) -> None:
    command = f"gzip -c {input_fastq} > {input_fastq}.gz"
    subprocess.run(command, shell=True, check=True)
    if not keep_original:
        subprocess.run(f"rm {input_fastq}", shell=True, check=True)

def mix_fastq_files(sample_config: str, n_in_mixture: int, proportions: list, output_dir: str) -> None:
    config = parse_json(sample_config, remove_nest=True)
    samples = config['samples']
    sets = make_sets_n_choose_k(list(samples.keys()), n_in_mixture)
    metadata = []
    for i, s in enumerate(sets):
        reads_to_sample_from_each = get_number_of_reads_to_sample([samples[sample]['R1'] for sample in s], proportions)
        output_prefix = f"{output_dir}/mixture_{i}"
        # make a blank fastq file first
        subprocess.run(f"> {output_prefix}_R1.fastq", shell=True, check=True)
        subprocess.run(f"> {output_prefix}_R2.fastq", shell=True, check=True)
        for sample, n_reads in zip(s, reads_to_sample_from_each):
            # now subsample the reads and append to the blank fastq files
            subsample_fastq(samples[sample]['R1'], f"{output_prefix}_R1.fastq", n_reads)
            subsample_fastq(samples[sample]['R2'], f"{output_prefix}_R2.fastq", n_reads)
        # gzip the output fastq files
        gzip_fastq(f"{output_prefix}_R1.fastq")
        gzip_fastq(f"{output_prefix}_R2.fastq")
        row = {"mixture_id": i}
        for j, (sample, proportion, n_reads) in enumerate(zip(s, proportions, reads_to_sample_from_each), start=1):
            row[f"sample_{j}"] = sample
            row[f"proportion_{j}"] = proportion
            row[f"reads_{j}"] = n_reads
        metadata.append(row)
    metadata_df = pl.DataFrame(metadata)
    metadata_df.write_csv(f"{output_dir}/mixture_metadata.tsv", separator='\t')
        
def make_snakemake_input_json(metadata_file: str, fastq_dir: str, output_file: str) -> None:
    metadata_df = pl.read_csv(metadata_file, separator='\t')
    output_dict = {"sample_sets": {"in_silico_mixtures": {"samples": {}}}}
    for sample in metadata_df.iter_rows(named=True):
        sample_id = sample['mixture_id']
        sample_fastq_files = {
            "R1": f"{fastq_dir}/mixture_{sample_id}_R1.fastq.gz",
            "R2": f"{fastq_dir}/mixture_{sample_id}_R2.fastq.gz",
        }
        output_dict["sample_sets"]["in_silico_mixtures"]["samples"][f"mixture_{sample_id}"] = sample_fastq_files
    with open(output_file, 'w') as f:
        json.dump(output_dict, f, indent=4)
