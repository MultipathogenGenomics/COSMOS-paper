# compute coverage of mt genomes in alignments

import subprocess
import os

def compute_coverage_of_mt_genomes_in_alignments(alignments_dir: str, output_dir: str) -> None:
    for alignment_file in os.listdir(alignments_dir):
        if alignment_file.endswith(".bam"):
            species = alignment_file.replace(".bam", "")
            input_bam = f"{alignments_dir}/{alignment_file}"
            output_coverage = f"{output_dir}/{species}_coverage.txt"
            command = f"bedtools genomecov -ibam {input_bam} -d | gzip -c > {output_coverage}.gz"            
            subprocess.run(command, shell=True, check=True)