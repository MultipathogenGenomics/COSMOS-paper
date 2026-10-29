import polars as pl
import subprocess
import io
import os
import numpy as np

def run_tabulate_alignments(sample_bam: str, position: int, tabulate_alignments_executable: str, target_id: str, mapq_threshold: int) -> pl.DataFrame:
    # use tabulate alignments to get pileup across kmer
    command = f"{tabulate_alignments_executable} -reads {sample_bam} -range {target_id}:{position}-{position} -mq-bins 0 {mapq_threshold}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    header = ["file","chromosome","position","mq>=0:A","mq>=0:C","mq>=0:G","mq>=0:T","mq>=0:-","mq>=0:+",f"mq>={mapq_threshold}:A",f"mq>={mapq_threshold}:C",f"mq>={mapq_threshold}:G",f"mq>={mapq_threshold}:T",f"mq>={mapq_threshold}:-",f"mq>={mapq_threshold}:+"] 
    df = pl.read_csv(io.StringIO(result.stdout), separator='\t', has_header=True)
    return df

def return_nucleotide_counts_from_tabulate_alignments(sample_bam: str, position: int, tabulate_alignments_executable: str, target_id: str, mapq_threshold: int) -> dict:
    df = run_tabulate_alignments(sample_bam, position, tabulate_alignments_executable, target_id, mapq_threshold)
    counts = {
        "A": df[f"mq>={mapq_threshold}:A"].sum(),
        "C": df[f"mq>={mapq_threshold}:C"].sum(),
        "G": df[f"mq>={mapq_threshold}:G"].sum(),
        "T": df[f"mq>={mapq_threshold}:T"].sum(),
        "-": df[f"mq>={mapq_threshold}:-"].sum(),
        "+": df[f"mq>={mapq_threshold}:+"].sum()
    }
    return counts

def compute_snp_ratio_at_site(sample_bam: str, position: int, nucleotideA: str, nucleotideB: str, tabulate_alignments_executable: str, target_id: str, mapq_threshold: int) -> float:
    read_count_df = run_tabulate_alignments(sample_bam, position, tabulate_alignments_executable, target_id, mapq_threshold)
    # find A, C, T, G counts with mapq >= threshold
    A_count = read_count_df[f"mq>={mapq_threshold}:{nucleotideA}"].sum()
    B_count = read_count_df[f"mq>={mapq_threshold}:{nucleotideB}"].sum()
    total_count = A_count + B_count
    return A_count / total_count if total_count > 0 else 0.0

def compute_snp_ratios_for_subtype_across_samples(snp_file: str, bam_dir: str, target_id: str, tabulate_alignments_executable: str) -> pl.DataFrame:
    snp_df = pl.read_csv(snp_file, separator='\t')
    results = []
    samples = os.listdir(bam_dir)
    for sample in samples:
        sample_bam = os.path.join(bam_dir, sample, f"{sample}_kma_mapped_to_genome.sorted.bam")
        for row in snp_df.iter_rows(named=True):
            position = row['position_in_reference'] 
            counts = return_nucleotide_counts_from_tabulate_alignments(sample_bam, position, tabulate_alignments_executable, target_id, mapq_threshold=20)
            counts["sample"] = sample
            counts["position_in_reference"] = position
            results.append(counts)
    counts_df = pl.DataFrame(results)
    results_df = counts_df.join(snp_df, on="position_in_reference", how="left")
    return results_df

def split_bases_with_commas(bases_list: list) -> list:
    split_bases = []
    for bases in bases_list:
        split_bases.extend(bases.split(","))
    return np.unique(np.asarray(split_bases)).tolist()

def make_results_readable(df: pl.DataFrame, snp_file: str) -> pl.DataFrame:
    # make a dict of info from snp file
    snp_df = pl.read_csv(snp_file, separator='\t')
    snp_info_dict = {}
    subtypes = snp_df['subtype'].unique().to_list()
    positions = snp_df['position_in_reference'].unique().to_list()
    for position in positions:
        snp_info_dict[position] = {}
        for subtype in subtypes:
            snp_info_dict[position][subtype] = {}
            subset = snp_df.filter((pl.col("position_in_reference") == position) & (pl.col("subtype") == subtype))
            bases = subset['subtype_base'].unique().to_list()
            non_subtype_bases = subset['non_subtype_bases'].unique().to_list()
            snp_info_dict[position][subtype]["bases"] = split_bases_with_commas(bases)
            snp_info_dict[position][subtype]["non_subtype_bases"] = split_bases_with_commas(non_subtype_bases)
    # now collect results
    all_bases = ["A", "C", "G", "T", "-", "+"]
    # remove duplicated rows based on all_bases, sample, position_in_reference as don't need here
    df = df.unique(subset=all_bases + ["sample", "position_in_reference"], keep="first")
    df = df.with_columns(pl.sum_horizontal(all_bases).alias("total_count"))
    samples = df['sample'].unique().to_list()
    rows = []
    for sample in samples:
        for position in snp_info_dict.keys():
            row = {"sample": sample, "position_in_reference": position}
            total_counts = df.filter((pl.col("sample") == sample) & (pl.col("position_in_reference") == position))["total_count"].unique()
            assert len(total_counts) == 1
            row["total_counts_at_position"] = total_counts[0]
            for subtype in snp_info_dict[position].keys():
                if "ann " in str(subtype):
                    subtype_label = str(subtype).replace("ann ", "")
                else:
                    subtype_label = str(subtype)
                count_for_subtype = 0
                count_for_non_subtype = 0
                for base in snp_info_dict[position][subtype]["bases"]:
                    base = base.upper()
                    count_for_subtype += df.filter((pl.col("sample") == sample) & (pl.col("position_in_reference") == position))[base].sum()
                for base in snp_info_dict[position][subtype]["non_subtype_bases"]:
                    base = base.upper()
                    count_for_non_subtype += df.filter((pl.col("sample") == sample) & (pl.col("position_in_reference") == position))[base].sum()
                row[f"{subtype_label}_subtype_bases"] = ",".join(snp_info_dict[position][subtype]["bases"])
                row[f"{subtype_label}_subtype_count"] = count_for_subtype
                row[f"{subtype_label}_non_subtype_bases"] = ",".join(snp_info_dict[position][subtype]["non_subtype_bases"])
                row[f"{subtype_label}_non_subtype_count"] = count_for_non_subtype
                row[f"{subtype_label}_ratio"] = round(count_for_subtype / row["total_counts_at_position"] if row["total_counts_at_position"] > 0 else 0.0, 3)
            # add any missing subtypes with 0 counts
            for subtype in subtypes:
                if "ann " in str(subtype):
                    subtype_label = str(subtype).replace("ann ", "")
                else:
                    subtype_label = str(subtype)
                if f"{subtype_label}_subtype_count" not in row:
                    row[f"{subtype_label}_subtype_bases"] = ""
                    row[f"{subtype_label}_subtype_count"] = 0
                if f"{subtype_label}_non_subtype_count" not in row:
                    row[f"{subtype_label}_non_subtype_bases"] = ""
                    row[f"{subtype_label}_non_subtype_count"] = 0
                if f"{subtype_label}_ratio" not in row:
                    row[f"{subtype_label}_ratio"] = 0.0
            rows.append(row)
    readable_df = pl.DataFrame(rows)
    return readable_df
