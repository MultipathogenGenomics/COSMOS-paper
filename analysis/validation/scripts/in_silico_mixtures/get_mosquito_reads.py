from analysis.common.software import REMOVECLIPPING_SCRIPT
import os
import subprocess
import json
import gzip
import sys
from analysis.common.io import parse_json

def ensure_fasta_index(fasta: str) -> None:
    if not os.path.exists(f"{fasta}.fai"):
        subprocess.run(f"samtools faidx {fasta}", shell=True, check=True)

def align_reads_with_minimap(r1: str, r2: str, genome: str, output_prefix: str) -> None:
    ensure_fasta_index(genome)
    command = f"minimap2 -t 16 -ax sr {genome} {r1} {r2} | samtools view -@ 16 -b -F 4 - | samtools sort -@ 16 -o {output_prefix}.bam"
    subprocess.run(command, shell=True, check=True)

def remove_softclipping_from_bam(input_bam: str, output_bam: str, script=REMOVECLIPPING_SCRIPT) -> None:
    subprocess.run([sys.executable, script, "-f", input_bam, output_bam], check=True)

def filter_bam_to_reads_entirely_within_region(input_bam: str, output_bam: str, contig_string: str) -> None:
    # make a bed file with the contig string
    os.makedirs("tmp", exist_ok=True)
    bedfile = "tmp/region.bed"
    with open(bedfile, 'w') as f:
        for contig in contig_string.split(" "):
            contig_name, region = contig.split(":")
            start, end = region.split("-")
            f.write(f"{contig_name}\t{int(start) - 1}\t{int(end)}\n")
    # filter the bam file to only include reads entirely within the region
    command = f"bedtools intersect -abam {input_bam} -b {bedfile} -f 1.0 -wa -u > {output_bam}"
    subprocess.run(command, shell=True, check=True)
    subprocess.run(f"rm {bedfile}", shell=True, check=True)

def extract_fastq_from_bam(input_bam: str, output_prefix: str) -> None:
        samtools_command = f"""
            samtools index {input_bam}
            samtools view -@ 16 -b -F 4 {input_bam} \
            | samtools sort -@ 16 -n - | bedtools bamtofastq -i - -fq {output_prefix}_R1.fastq -fq2 {output_prefix}_R2.fastq
            """
        subprocess.run(samtools_command, shell=True, check=True)

def filter_fastq_by_read_length(input_fastq: str, output_fastq: str, min_length: int) -> None:
    command = f"fastq-filter -l {min_length} -o {output_fastq} {input_fastq}"
    subprocess.run(command, shell=True, check=True)

def _normalise_read_id(header: str) -> str:
    return header.split()[0].removeprefix("@").removesuffix("/1").removesuffix("/2")

def filter_paired_fastq_by_read_length(r1_fastq: str, r2_fastq: str, output_r1: str, output_r2: str, min_length: int) -> None:
    total = 0
    kept = 0
    with open(r1_fastq) as r1_in, open(r2_fastq) as r2_in, open(output_r1, "w") as r1_out, open(output_r2, "w") as r2_out:
        while True:
            r1_record = [r1_in.readline() for _ in range(4)]
            r2_record = [r2_in.readline() for _ in range(4)]
            if not any(r1_record) and not any(r2_record):
                break
            if not all(r1_record) or not all(r2_record):
                raise ValueError(f"FASTQ files have incomplete or unequal records: {r1_fastq}, {r2_fastq}")
            if _normalise_read_id(r1_record[0]) != _normalise_read_id(r2_record[0]):
                raise ValueError(f"FASTQ files are not paired at record {total + 1}: {r1_record[0].strip()} vs {r2_record[0].strip()}")
            total += 1
            if len(r1_record[1].strip()) >= min_length and len(r2_record[1].strip()) >= min_length:
                r1_out.writelines(r1_record)
                r2_out.writelines(r2_record)
                kept += 1
    print(f"{r1_fastq}, {r2_fastq}: kept {kept}/{total} read pairs with both mates >= {min_length} bp")

def gzip_fastq(input_fastq: str, output_fastq: str) -> None:
    command = f"gzip -c {input_fastq} > {output_fastq}"
    subprocess.run(command, shell=True, check=True)

def get_num_reads_in_fastq(fastq_file: str) -> int:
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as f:
        return sum(1 for _ in f) // 4

def print_num_reads_in_fastq(fastq_file: str) -> None:
    num_reads = get_num_reads_in_fastq(fastq_file)
    print(f"{fastq_file}: {num_reads} reads")

def get_mosquito_reads_from_ss_pools(samples_config_file: str, reference_genome: str, output_dir: str, contig_string: str) -> None:
    config = parse_json(samples_config_file, remove_nest=True)
    samples = config['samples']
    for sample_name, sample_info in samples.items():
        r1 = sample_info['R1']
        r2 = sample_info['R2']
        output_prefix = f"{output_dir}/{sample_name}"
        align_reads_with_minimap(r1, r2, reference_genome, output_prefix)
        remove_softclipping_from_bam(f"{output_prefix}.bam", f"{output_prefix}_no_softclip.bam")
        filter_bam_to_reads_entirely_within_region(f"{output_prefix}_no_softclip.bam", f"{output_prefix}_no_softclip_region.bam", contig_string)
        extract_fastq_from_bam(f"{output_prefix}_no_softclip_region.bam", output_prefix)
        print_num_reads_in_fastq(f"{output_prefix}_R1.fastq")
        print_num_reads_in_fastq(f"{output_prefix}_R2.fastq")
        filter_paired_fastq_by_read_length(
            f"{output_prefix}_R1.fastq",
            f"{output_prefix}_R2.fastq",
            f"{output_prefix}_R1.filtered.fastq",
            f"{output_prefix}_R2.filtered.fastq",
            min_length=100,
        )
        print_num_reads_in_fastq(f"{output_prefix}_R1.filtered.fastq")
        gzip_fastq(f"{output_prefix}_R1.filtered.fastq", f"{output_prefix}_R1.fastq.gz")
        print_num_reads_in_fastq(f"{output_prefix}_R2.filtered.fastq")
        gzip_fastq(f"{output_prefix}_R2.filtered.fastq", f"{output_prefix}_R2.fastq.gz")
        # remove intermediate bam files
        os.remove(f"{output_prefix}.bam")
        os.remove(f"{output_prefix}_no_softclip.bam")
        os.remove(f"{output_prefix}_no_softclip_region.bam")
        os.remove(f"{output_prefix}_no_softclip_region.bam.bai")
        # remove intermediate fastq files
        os.remove(f"{output_prefix}_R1.fastq")
        os.remove(f"{output_prefix}_R1.filtered.fastq")
        os.remove(f"{output_prefix}_R2.fastq")
        os.remove(f"{output_prefix}_R2.filtered.fastq")

def make_json_of_mosquito_fastq_files(original_config_file: str, new_fastq_dir: str, output_json: str) -> None:
    config = parse_json(original_config_file, remove_nest=False)
    first_nest = list(config.keys())[0] + "-mosquito-cox1-extracted"
    original_first_nest = list(config.keys())[0]
    config[first_nest] = config.pop(original_first_nest)
    samples = config[first_nest]['samples']
    samples_to_remove = []
    for sample in samples:
        new_R1 = f"{new_fastq_dir}/{sample}_R1.fastq.gz"
        new_R2 = f"{new_fastq_dir}/{sample}_R2.fastq.gz"
        # check existence of new fastq files
        if not os.path.exists(new_R1) or not os.path.exists(new_R2):
            raise FileNotFoundError(f"New fastq files for sample {sample} not found in {new_fastq_dir}")
        samples[sample]['R1'] = new_R1
        samples[sample]['R2'] = new_R2
        # if n reads in fastq < 2000, remove from config
        # num_reads = get_num_reads_in_fastq(new_R1)
    #     if num_reads < 2000:
    #         print(f"Sample {sample} has only {num_reads} reads in R1 fastq, removing from config")
    #         samples_to_remove.append(sample)
    # for sample in samples_to_remove:
    #     del samples[sample]
    with open(output_json, 'w') as f:
        json.dump(config, f, indent=4)
