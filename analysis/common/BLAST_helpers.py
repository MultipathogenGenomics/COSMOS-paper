import polars as pl
import os
import tempfile
import subprocess
import numpy as np
import io
from Bio import SeqIO

def parse_tsv(tsv_path: str) -> pl.DataFrame:
    return pl.read_csv(tsv_path, separator='\t')

def join_kmer_count_tables(list_of_sample_kmer_files: list) -> pl.DataFrame:
    results = None
    for kmer_file in list_of_sample_kmer_files:
        sample_name = os.path.basename(kmer_file).replace('.counts_as_table.tsv', '')
        df = parse_tsv(kmer_file)
        df = df.rename({"count": sample_name})
        if results is None:
            results = df
        else:
            results = results.join(df, on="kmer", how="outer", coalesce=True)
    # Fill NaN with 0
    results = results.fill_null(0)
    return results

def filter_to_sample_defining_kmers(df: pl.DataFrame, min_count=2, max_count=50) -> pl.DataFrame:
    # filter to where only one sample has non-zero count
    sample_cols = [col for col in df.columns if col != 'kmer']
    non_zero_count = sum((pl.col(col) > 0).cast(pl.Int32) for col in sample_cols) # find where only one sample has non-zero count
    max_more_than_min_count = pl.max_horizontal([pl.col(col) for col in sample_cols]) >= min_count
    max_less_than_max_count = pl.max_horizontal([pl.col(col) for col in sample_cols]) <= max_count
    combined_filter = (non_zero_count == 1) & (max_more_than_min_count) & (max_less_than_max_count)
    df = df.filter(combined_filter)
    return df

def make_fasta_from_kmers(df: pl.DataFrame, output_fasta: str):
    with open(output_fasta, 'w') as f:
        for col in df.columns:
            if col == 'kmer':
                continue
            sample_df = df.filter(pl.col(col) > 0)
            n = 1
            for row in sample_df.iter_rows():
                kmer = row[0]
                f.write(f">{col}_{n}\n{kmer}\n")
                n += 1
    f.close()

def process_unique_kmers(list_of_sample_kmer_files: list, output_table: str, min_count=50, max_count=100):
    combined_df = join_kmer_count_tables(list_of_sample_kmer_files)
    unique_kmers_df = filter_to_sample_defining_kmers(combined_df, min_count=min_count, max_count=max_count)
    unique_kmers_df.write_csv(output_table, separator='\t')
    output_fasta = output_table.replace('.tsv', '.fasta')
    make_fasta_from_kmers(unique_kmers_df, output_fasta)

def make_blast_db(fasta_path: str, db_name: str):
    os.system(f"makeblastdb -in {fasta_path} -dbtype nucl -out {db_name}")

def blast_kmer(kmer_str: str, blast_db: str, output_blast: str):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_fasta:
        tmp_fasta.write(">query\n")
        tmp_fasta.write(kmer_str.strip() + "\n")
        tmp_fasta_path = tmp_fasta.name
    subprocess.run(
            [
                "blastn",
                "-task", "blastn-short",
                "-query",
                tmp_fasta_path,
                "-db",
                blast_db,
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                "-out",
                output_blast,
            ],
            check=True,
        )
    os.remove(tmp_fasta_path)

def extract_best_blast_hit(kmer_str: str, output_blast: str, target_of_interest: str) -> str:
    if not os.path.exists(output_blast) or os.stat(output_blast).st_size == 0:
        return None, None, None
    blast_df = pl.read_csv(output_blast, separator='\t', has_header=False)
    blast_df = blast_df.rename({
        "column_1": "qseqid",
        "column_2": "sseqid",
        "column_3": "pident",
        "column_4": "length",
        "column_5": "mismatch",
        "column_6": "gapopen",
        "column_7": "qstart",
        "column_8": "qend",
        "column_9": "sstart",
        "column_10": "send",
        "column_11": "evalue",
        "column_12": "bitscore"
    })
    # filter for hits to target of interest
    blast_df = blast_df.filter(pl.col("sseqid").str.contains(target_of_interest))
    if blast_df.is_empty():
        return None, None, None
    best_hit = blast_df.sort("bitscore", descending=True)[0] # take best hit
    start_in_target = min(best_hit["sstart"].item(), best_hit["send"].item())
    end_in_target = max(best_hit["sstart"].item(), best_hit["send"].item())
    # if hit doesn't cover whole kmer, extend to cover whole kmer
    kmer_length = len(kmer_str.strip())
    if best_hit["length"].item() < kmer_length:
        if best_hit["sstart"].item() < best_hit["send"].item():
            best_start = max(1, start_in_target - (kmer_length - best_hit["length"].item()))
            best_end = end_in_target
        else:
            best_start = start_in_target
            best_end = end_in_target + (kmer_length - best_hit["length"].item())
    else:
        best_start = start_in_target
        best_end = end_in_target
    return best_hit["sseqid"].item(), best_start, best_end
        
def check_kmer_not_polymorphic_in_sample(sample_bam: str, target_id: str, start: int, end: int, tabulate_alignments_executable: str, mapq_threshold: int=10, polymorphism_threshold: float=0.95) -> bool:
    # use tabulate alignments to get pileup across kmer
    command = f"{tabulate_alignments_executable} -reads {sample_bam} -range {target_id}:{start}-{end} -mq-bins 0 {mapq_threshold}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    header = ["file","chromosome","position","mq>=0:A","mq>=0:C","mq>=0:G","mq>=0:T","mq>=0:-","mq>=0:+",f"mq>={mapq_threshold}:A",f"mq>={mapq_threshold}:C",f"mq>={mapq_threshold}:G",f"mq>={mapq_threshold}:T",f"mq>={mapq_threshold}:-",f"mq>={mapq_threshold}:+"] 
    df = pl.read_csv(io.StringIO(result.stdout), separator='\t', has_header=True)
    # if one nucleotide, where mapq greater than threshold, does not carry the polymorpshism threshold proportion then kmer fails
    desired_cols = [pl.col(col) for col in df.columns if col.startswith(f"mq>={mapq_threshold}:") and col not in ["+","-"] ]
    df = df.with_columns(
        (pl.max_horizontal(desired_cols) / pl.sum_horizontal(desired_cols)).alias("max_allele_proportion")
    )
    mean_allele_proportion = np.nanmean(df["max_allele_proportion"].to_numpy())
    if (df["max_allele_proportion"] < polymorphism_threshold).any() or np.isnan(mean_allele_proportion):
        return False, mean_allele_proportion
    else:
        return True, mean_allele_proportion
    
def check_kmers(
        unique_kmers_file: str, 
        targets_fasta: str, 
        sample_bams: list, 
        target_of_interest: list, 
        tabulate_alignments_executable: str, 
        mapq_threshold: int=10, 
        polymorphism_threshold: float=0.95, 
        path_for_blast: str="targets_db",
        output_path: str="output",
        skip_blast: bool=False) -> pl.DataFrame:
    unique_kmers_df = pl.read_csv(unique_kmers_file, separator='\t')
    make_blast_db(targets_fasta, f"{path_for_blast}/blast_db")
    results = []
    samples = unique_kmers_df.columns[1:] # all columns except kmer column are sample columns
    for row in unique_kmers_df.iter_rows(named=True):
        kmer = row['kmer']
        sample = [col for col in samples if int(row[col]) > 0][0]
        sample_index = samples.index(sample)
        bam = sample_bams[sample_index] # find bam for sample
        target = target_of_interest[sample_index]
        blast_output = f"{path_for_blast}/{sample}_{kmer}.blast"
        if not skip_blast:
            blast_kmer(kmer, f"{path_for_blast}/blast_db", blast_output)
            best_hit, start, end = extract_best_blast_hit(kmer, blast_output, target)
            os.remove(blast_output)
        else:
            best_hit = None
        if best_hit is None:
            best_hit = "no_hit"
            start = None
            end = None
            is_polymorphic = None
            average_allele_proportion = None
        else:
            is_polymorphic, average_allele_proportion = check_kmer_not_polymorphic_in_sample(bam, best_hit, start, end, tabulate_alignments_executable, mapq_threshold, polymorphism_threshold)
        results.append({
            "kmer": kmer,
            "sample": sample,
            "target": best_hit,
            "start": start,
            "end": end,
            "is_polymorphic": is_polymorphic,
            "average_allele_proportion": average_allele_proportion
        })
    results_df = pl.DataFrame(results)
    results_df.write_csv(output_path, separator='\t')
    return results_df

def minimap_reads_to_target(reads_folder: str, target_fasta: str, target: str, output_bam: str, minimap_executable: str = "minimap2", tmp_folder: str = "/tmp"):
    # extract target from fasta using samtools faidx
    target_fasta_index = f"{target_fasta}.fai"
    if not os.path.exists(target_fasta_index):
        subprocess.run(f"samtools faidx {target_fasta}", shell=True, check=True)
    target_fasta_tmp = os.path.join(tmp_folder, f"{target}_temp.fasta")
    subprocess.run(f"samtools faidx {target_fasta} {target} > {target_fasta_tmp}", shell=True, check=True)
    # minimap reads to target
    subprocess.run(f"{minimap_executable} -ax sr {target_fasta_tmp} {reads_folder}/*.fastq.gz | samtools view -bS - | samtools sort -o {output_bam}", shell=True, check=True)
    os.remove(target_fasta_tmp)

def count_kmers_in_mixture(kmers_df: pl.DataFrame, mixture_kmer_counts_file: str) -> pl.DataFrame:
    mixture_counts_df = pl.read_csv(mixture_kmer_counts_file, separator='\t')
    merged_df = kmers_df.join(mixture_counts_df, on="kmer", how="left")
    merged_df = merged_df.rename({"count": "mixture_count"})
    merged_df = merged_df.fill_null(0)
    return merged_df

def bootstrap_median_ratio_between_mixture_and_sample(samples: list, df: pl.DataFrame, bootstrap_size=1000, n_bootstrap=1000):
    results = []
    for i in range(n_bootstrap):
        bootstrap_result = {}
        bootstrap_result["n_bootstrap"] = i
        median_ratios = []
        for sample in samples:
            sample_df = df.filter(pl.col(sample) > 0)
            sample_df = sample_df.with_columns((pl.col("mixture_count") / pl.col(f"{sample}_normalised")).alias("count_ratio"))
            bootstrap_sample = sample_df.sample(bootstrap_size, with_replacement=True)
            median_ratio = bootstrap_sample['count_ratio'].median()
            median_ratios.append(median_ratio)
        median_ratios = np.array(median_ratios)
        for sample in samples:
            bootstrap_result[f"{sample}_median_ratio"] = median_ratios[samples.index(sample)]
            bootstrap_result[f"{sample}_abundance"] = median_ratios[samples.index(sample)] / np.sum(median_ratios)
        results.append(bootstrap_result)
    results_df = pl.DataFrame(results)
    return results_df