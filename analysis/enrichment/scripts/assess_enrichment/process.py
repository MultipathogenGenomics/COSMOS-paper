# sort coverage outputs into plotting format

import os
import subprocess

import numpy as np
import polars as pl

from analysis.common.io import parse_json

def read_fasta(fasta_file: str) -> dict[str, str]:
    records = {}
    current_name = None
    current_sequence = []
    with open(fasta_file) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_sequence)
                current_name = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)
    if current_name is not None:
        records[current_name] = "".join(current_sequence)
    return records

def build_alignment_position_map(aligned_fasta: str) -> pl.DataFrame:
    rows = []
    for ref_genome, aligned_sequence in read_fasta(aligned_fasta).items():
        position = 0
        for aligned_position, base in enumerate(aligned_sequence, start=1):
            if base == "-":
                continue
            position += 1
            rows.append({
                "ref_genome": ref_genome,
                "position": position,
                "aligned_position": aligned_position,
            })
    return pl.DataFrame(rows)

def add_alignment_positions_to_coverage_data(coverage_df: pl.DataFrame, aligned_fasta: str) -> pl.DataFrame:
    alignment_position_map = build_alignment_position_map(aligned_fasta)
    coverage_with_alignment = coverage_df.join(
        alignment_position_map,
        on=["ref_genome", "position"],
        how="left",
    )
    missing_alignment_positions = coverage_with_alignment.filter(pl.col("aligned_position").is_null())
    if missing_alignment_positions.height > 0:
        missing_refs = missing_alignment_positions.select("ref_genome").unique().to_series().to_list()
        raise ValueError(f"Could not map some coverage positions to aligned positions for references: {missing_refs}")
    return coverage_with_alignment

def concat_coverage_data(sample_config: str, genomes_config: str, coverage_dir: str, genus_defaults_dict: dict) -> pl.DataFrame:
    sample_config_data = parse_json(sample_config, remove_nest=True)
    genomes_config_data = parse_json(genomes_config, remove_nest=True)
    coverage_dfs = []
    for sample in sample_config_data["samples"]:
        if sample in genomes_config_data["references"].keys():
            ref_genome = sample
            ref_is_species = True
        else:
            genus = sample.split("_")[0]
            if genus in genus_defaults_dict:
                ref_genome_file = genus_defaults_dict[genus]
                ref_genome_file_basename = os.path.basename(ref_genome_file)
                ref_genome = ref_genome_file_basename.split("_mt_genome.fasta")[0]
                ref_is_species = False
            else:
                raise ValueError(f"No reference genome found for sample {sample} and no default for genus {genus}")
        cox1_start = genomes_config_data["references"][ref_genome]["cox1_start"]
        cox1_end = genomes_config_data["references"][ref_genome]["cox1_end"]
        coverage_file = f"{coverage_dir}/{sample}_coverage.txt.gz"
        coverage_data = pl.read_csv(coverage_file, separator="\t")
        coverage_data.columns = ["contig", "position", "coverage"]
        coverage_data = coverage_data.with_columns(
            pl.lit(sample).alias("sample"),
            pl.lit(ref_genome).alias("ref_genome"),
            pl.lit(ref_is_species).alias("ref_is_species"),
        )
        coverage_data = coverage_data.with_columns(
            pl.when((pl.col("position") < cox1_end) & (pl.col("position") >= cox1_start))
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("is_cox1")
        )
        coverage_dfs.append(coverage_data)
    results_df = pl.concat(coverage_dfs)
    return results_df

def compute_median_mean_in_cox1(df: pl.DataFrame, buffer=250) -> pl.DataFrame:
    results = []
    for sample in df['sample'].unique():
        sample_df = df.filter(pl.col("sample") == sample)
        cox1_start = sample_df.filter(pl.col("is_cox1"))['position'].min()
        cox1_end = sample_df.filter(pl.col("is_cox1"))['position'].max()
        median_in_cox1 = sample_df.filter(pl.col("is_cox1")).select(pl.col("coverage").median()).item()
        median_out_cox1 = sample_df.filter((pl.col("position") < cox1_start - buffer) | (pl.col("position") > cox1_end + buffer)).select(pl.col("coverage").median()).item()
        mean_in_cox1 = sample_df.filter(pl.col("is_cox1")).select(pl.col("coverage").mean()).item()
        mean_out_cox1 = sample_df.filter((pl.col("position") < cox1_start - buffer) | (pl.col("position") > cox1_end + buffer)).select(pl.col("coverage").mean()).item()
        results.append({
            "sample": sample,
            "median_in_cox1": median_in_cox1,
            "median_out_cox1": median_out_cox1,
            "mean_in_cox1": mean_in_cox1,
            "mean_out_cox1": mean_out_cox1
        })
    results_df = pl.DataFrame(results)
    results_df = results_df.with_columns(
        (pl.col("mean_in_cox1") / pl.col("mean_out_cox1")).alias("mean_fold_enrichment"),
        (pl.col("median_in_cox1") / pl.col("median_out_cox1")).alias("median_fold_enrichment"),
    )
    return results_df


def count_bam_alignments(bam_file: str, region: str | None = None) -> int:
    command = ["samtools", "view", "-c", "-F", "4", bam_file]
    if region is not None:
        command.append(region)
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def compute_enrichment_summary_table(
    coverage_df: pl.DataFrame,
    alignment_map_file: str,
    genomes_config: str,
    output_file: str,
    buffer: int = 250,
) -> pl.DataFrame:
    genomes_config_data = parse_json(genomes_config, remove_nest=True)
    alignment_map_df = pl.read_csv(alignment_map_file, separator="\t")

    rows = []
    for alignment_row in alignment_map_df.to_dicts():
        sample = alignment_row["sample"]
        reference_name = alignment_row["reference_name"]
        sample_coverage_df = coverage_df.filter(pl.col("sample") == sample)
        if sample_coverage_df.is_empty():
            continue

        cox1_start = int(genomes_config_data["references"][reference_name]["cox1_start"])
        cox1_end = int(genomes_config_data["references"][reference_name]["cox1_end"])
        contig = sample_coverage_df.select("contig").to_series().item(0)
        cox1_region = f"{contig}:{cox1_start}-{cox1_end}"

        in_cox1_df = sample_coverage_df.filter(pl.col("is_cox1"))
        out_cox1_df = sample_coverage_df.filter(
            (pl.col("position") < cox1_start - buffer)
            | (pl.col("position") > cox1_end + buffer)
        )

        n_reads_aligned = count_bam_alignments(alignment_row["bam"])
        n_reads_in_cox1 = count_bam_alignments(alignment_row["bam"], cox1_region)
        n_reads_out_of_cox1 = n_reads_aligned - n_reads_in_cox1
        mean_coverage_in_cox1 = in_cox1_df.select(pl.col("coverage").mean()).item()
        mean_coverage_out_of_cox1 = out_cox1_df.select(pl.col("coverage").mean()).item()
        median_coverage_in_cox1 = in_cox1_df.select(pl.col("coverage").median()).item()
        median_coverage_out_of_cox1 = out_cox1_df.select(pl.col("coverage").median()).item()
        rows.append(
            {
                "species_in_pool": sample,
                "reference_aligned_to": reference_name,
                "ref_is_species": alignment_row["ref_is_species"],
                "n_reads_aligned": n_reads_aligned,
                "n_reads_in_cox1": n_reads_in_cox1,
                "n_reads_out_of_cox1": n_reads_out_of_cox1,
                "proportion_reads_in_cox1": n_reads_in_cox1 / n_reads_aligned if n_reads_aligned > 0 else None,
                "proportion_reads_out_of_cox1": n_reads_out_of_cox1 / n_reads_aligned if n_reads_aligned > 0 else None,
                "mean_coverage_in_cox1": mean_coverage_in_cox1,
                "mean_coverage_out_of_cox1": mean_coverage_out_of_cox1,
                "mean_fold_enrichment": mean_coverage_in_cox1 / mean_coverage_out_of_cox1 if mean_coverage_out_of_cox1 > 0 else None,
                "median_coverage_in_cox1": median_coverage_in_cox1,
                "median_coverage_out_of_cox1": median_coverage_out_of_cox1,
                "median_fold_enrichment": median_coverage_in_cox1 / median_coverage_out_of_cox1 if median_coverage_out_of_cox1 > 0 else None,
                "sd_coverage_in_cox1": in_cox1_df.select(pl.col("coverage").std()).item(),
                "sd_coverage_out_of_cox1": out_cox1_df.select(pl.col("coverage").std()).item(),
            }
        )

    summary_df = pl.DataFrame(rows).sort("species_in_pool")
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    summary_df.write_csv(output_file, separator="\t")
    return summary_df


def compute_log_log_fit(summary_df: pl.DataFrame, x: str, y: str) -> dict:
    fit_df = summary_df.select([x, y]).drop_nulls()
    fit_df = fit_df.filter((pl.col(x) > 0) & (pl.col(y) > 0))
    if fit_df.height < 3:
        return {
            "n": fit_df.height,
            "slope": None,
            "intercept": None,
            "coefficient": None,
            "r_squared": None,
            "formula_log10": None,
            "formula_original_scale": None,
            "fold_enrichment_at_mean_x": None,
            "fold_enrichment_at_median_x": None,
        }

    x_values = fit_df[x].to_numpy()
    y_values = fit_df[y].to_numpy()
    log_x = np.log10(x_values)
    log_y = np.log10(y_values)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    fitted_log_y = intercept + slope * log_x
    residual_sum_squares = np.sum((log_y - fitted_log_y) ** 2)
    total_sum_squares = np.sum((log_y - log_y.mean()) ** 2)
    r_squared = 1 - (residual_sum_squares / total_sum_squares) if total_sum_squares > 0 else None
    coefficient = 10 ** intercept

    mean_x = float(np.mean(x_values))
    median_x = float(np.median(x_values))
    mean_y_fit = coefficient * (mean_x ** slope)
    median_y_fit = coefficient * (median_x ** slope)
    return {
        "n": fit_df.height,
        "slope": float(slope),
        "intercept": float(intercept),
        "coefficient": float(coefficient),
        "r_squared": None if r_squared is None else float(r_squared),
        "formula_log10": f"log10({y}) = {intercept:.6g} + {slope:.6g} * log10({x})",
        "formula_original_scale": f"{y} = {coefficient:.6g} * {x}^{slope:.6g}",
        "fold_enrichment_at_mean_x": float(mean_y_fit / mean_x),
        "fold_enrichment_at_median_x": float(median_y_fit / median_x),
    }


def write_overall_enrichment_summary_table(summary_df: pl.DataFrame, output_file: str) -> pl.DataFrame:
    numeric_cols = [
        "n_reads_aligned",
        "n_reads_in_cox1",
        "n_reads_out_of_cox1",
        "proportion_reads_in_cox1",
        "proportion_reads_out_of_cox1",
        "mean_coverage_in_cox1",
        "mean_coverage_out_of_cox1",
        "mean_fold_enrichment",
        "median_coverage_in_cox1",
        "median_coverage_out_of_cox1",
        "median_fold_enrichment",
        "sd_coverage_in_cox1",
        "sd_coverage_out_of_cox1",
    ]
    rows = []
    for col in numeric_cols:
        q1 = summary_df.select(pl.col(col).quantile(0.25)).item()
        q3 = summary_df.select(pl.col(col).quantile(0.75)).item()
        rows.append(
            {
                "summary_type": "metric_distribution",
                "metric": col,
                "n_alignments": summary_df.select(pl.col(col).drop_nulls().len()).item(),
                "mean": summary_df.select(pl.col(col).mean()).item(),
                "sd": summary_df.select(pl.col(col).std()).item(),
                "median": summary_df.select(pl.col(col).median()).item(),
                "q1": q1,
                "q3": q3,
                "iqr": q3 - q1 if q1 is not None and q3 is not None else None,
                "value": "",
            }
        )

    fit = compute_log_log_fit(
        summary_df,
        x="mean_coverage_out_of_cox1",
        y="mean_coverage_in_cox1",
    )
    for metric in [
        "formula_log10",
        "formula_original_scale",
        "n",
        "slope",
        "intercept",
        "coefficient",
        "r_squared",
        "fold_enrichment_at_mean_x",
        "fold_enrichment_at_median_x",
    ]:
        rows.append(
            {
                "summary_type": "mean_coverage_log_log_fit",
                "metric": metric,
                "n_alignments": fit["n"],
                "mean": None,
                "sd": None,
                "median": None,
                "q1": None,
                "q3": None,
                "iqr": None,
                "value": "" if fit[metric] is None else str(fit[metric]),
            }
        )

    overall_df = pl.DataFrame(rows)
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    overall_df.write_csv(output_file, separator="\t")
    return overall_df

