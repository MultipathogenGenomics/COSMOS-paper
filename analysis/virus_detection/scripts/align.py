
import os
import subprocess
from analysis.common.io import parse_json
from analysis.virus_detection.scripts.gather_seqs import find_sample_entry_in_multikey_dict

def ensure_fasta_index(fasta: str) -> None:
    if not os.path.exists(f"{fasta}.fai"):
        subprocess.run(f"samtools faidx {fasta}", shell=True, check=True)

def align_reads_with_minimap(r1: str, r2: str, genome: str, output_prefix: str) -> None:
    ensure_fasta_index(genome)
    command = f"minimap2 -t 16 -ax sr {genome} {r1} {r2} | samtools view -@ 16 -b -F 4 - | samtools sort -@ 16 -o {output_prefix}.bam"
    subprocess.run(command, shell=True, check=True)

def make_consensus_with_viralconsensus(input_bam: str, reference_fasta: str, output_fasta: str) -> None:
    ensure_fasta_index(reference_fasta)
    command = f"viral_consensus -i {input_bam} -r {reference_fasta} -o {output_fasta} --min_depth 2"
    subprocess.run(command, shell=True, check=True)

def make_consensus_from_input_ref(sample: str, sample_config: dict, input_ref: str, output_dir: str, output_base: str) -> None:
    sample_config = parse_json(sample_config, remove_nest=True)
    sample_entry = find_sample_entry_in_multikey_dict(sample_config, sample)
    r1 = sample_entry['R1']
    r2 = sample_entry['R2']
    output_prefix = f"{output_dir}/{output_base}"
    align_reads_with_minimap(r1, r2, input_ref, output_prefix)
    make_consensus_with_viralconsensus(f"{output_prefix}.bam", input_ref, f"{output_prefix}.min_depth_2.fasta")

def compute_coverage_at_each_position_using_samtools(input_bam: str, output_tsv: str) -> None:
    command = f"samtools depth -aa {input_bam} > {output_tsv}"
    subprocess.run(command, shell=True, check=True)