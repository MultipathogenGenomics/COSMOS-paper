from analysis.common.paths import RESULTS_DIR
# process and combine virus data for plotting

import os
from pathlib import Path

import polars as pl
from analysis.virus_detection.scripts.align import compute_coverage_at_each_position_using_samtools
from analysis.common.io import parse_json
from analysis.virus_detection.scripts.gather_seqs import find_sample_entry_in_multikey_dict
import subprocess

def filter_pcr_positive(data: pl.DataFrame) -> pl.DataFrame:
    return data.filter(pl.col("ct_value").cast(pl.Float64, strict=False).fill_null(0) > 0)


def filter_dedup_read_positive(data: pl.DataFrame) -> pl.DataFrame:
    return data.filter(pl.col("dedup_reads").cast(pl.Float64, strict=False).fill_null(0) > 0)


def _get_first_value_or_none(data: pl.DataFrame, sample: str, sample_column: str, value_column: str):
    matched_rows = data.filter(pl.col(sample_column) == sample)
    if matched_rows.is_empty():
        return None
    return matched_rows.select(value_column).to_series().item()

def combine_ct_dedup_and_nc2_data(cosmos_metadata_file: str, nc2_data_file: str, dedup_reads_data_file: str, virus_acronyms: dict, nc2_groups: list, nc2_labels: list) -> pl.DataFrame:
    cosmos_metadata = pl.read_csv(cosmos_metadata_file, separator="\t")
    nc2_data = pl.read_csv(nc2_data_file, separator="\t")
    dedup_reads_data = pl.read_csv(dedup_reads_data_file, separator="\t")
    # change SINDV ct values to NaN as PCR was not performed for this virus
    if "ct_SINDV" in cosmos_metadata.columns:
        cosmos_metadata = cosmos_metadata.with_columns(pl.lit(None).alias("ct_SINDV"))
    # add an empty col for LSV Cts as PCR was not performed for this virus
    if "ct_LSV" not in cosmos_metadata.columns:
        cosmos_metadata = cosmos_metadata.with_columns(pl.lit(None).alias("ct_LSV"))
    # make a format of sample, virus, ct, dedup_reads, nc2
    samples = cosmos_metadata.select("sample").unique().to_series().to_list()
    results = []
    for sample in samples:
        for virus, acronym in virus_acronyms.items():
            ct_col = f"ct_{acronym}"
            ct_value = _get_first_value_or_none(cosmos_metadata, sample, "sample", ct_col)
            dedup_reads_value = _get_first_value_or_none(dedup_reads_data, sample, "sampleid", virus)
            nc2_value = _get_first_value_or_none(nc2_data, sample, "sampleid", virus)
            results.append({"sample": sample, "virus": virus, "ct_value": ct_value, "dedup_reads": dedup_reads_value, "nc2": nc2_value})
    # add a column of genome coverage label based on nc2 value is in group
    df = pl.DataFrame(results)
    df = add_genome_coverage_category(df, nc2_groups, nc2_labels)
    return df


def add_genome_coverage_category(data: pl.DataFrame, nc2_groups: list, nc2_labels: list[str]) -> pl.DataFrame:
    if len(nc2_groups) != len(nc2_labels):
        raise ValueError("nc2_groups and nc2_labels must be the same length")
    return data.with_columns(
        pl.col("nc2").cast(pl.Float64, strict=False).fill_null(0).map_elements(
            lambda nc2: next(
                (
                    nc2_labels[index]
                    for index, (lower, upper) in enumerate(nc2_groups)
                    if lower <= nc2 < upper or (index == len(nc2_groups) - 1 and lower <= nc2 <= upper)
                ),
                nc2_labels[0],
            ),
            return_dtype=pl.Utf8,
        ).alias("genome_coverage_category")
    )


def compute_genome_coverage_category_counts(data: pl.DataFrame, nc2_groups: list, nc2_labels: list[str], virus_acronyms: dict) -> pl.DataFrame:
    category_data = add_genome_coverage_category(data, nc2_groups, nc2_labels)
    counts = category_data.group_by(["virus", "genome_coverage_category"]).agg(pl.len().alias("count"))
    pivot = counts.pivot(index="virus", columns="genome_coverage_category", values="count").fill_null(0)
    for label in nc2_labels:
        if label not in pivot.columns:
            pivot = pivot.with_columns(pl.lit(0).alias(label))
    return (
        pivot.select(["virus"] + nc2_labels)
        .with_columns(pl.col("virus").replace(virus_acronyms).alias("virus_acronym"))
        .select(["virus", "virus_acronym"] + nc2_labels)
        .sort("virus_acronym")
    )


def _diagnostic_summary_row(data: pl.DataFrame, virus: str | None = None) -> dict:
    if virus is not None:
        data = data.filter(pl.col("virus") == virus)
    pcr_positive = pl.col("pcr_positive")
    capture_positive = pl.col("capture_positive")
    summary = data.select(
        ((pcr_positive) & (capture_positive)).sum().alias("true_positive"),
        ((pcr_positive) & (~capture_positive)).sum().alias("false_negative"),
        ((~pcr_positive) & (capture_positive)).sum().alias("false_positive"),
        ((~pcr_positive) & (~capture_positive)).sum().alias("true_negative"),
    ).row(0, named=True)
    tp = summary["true_positive"]
    fn = summary["false_negative"]
    fp = summary["false_positive"]
    tn = summary["true_negative"]
    sensitivity_denominator = tp + fn
    specificity_denominator = tn + fp
    return {
        "virus": virus if virus is not None else "overall",
        "true_positive": tp,
        "false_negative": fn,
        "false_positive": fp,
        "true_negative": tn,
        "sensitivity": tp / sensitivity_denominator if sensitivity_denominator > 0 else None,
        "specificity": tn / specificity_denominator if specificity_denominator > 0 else None,
        "pcr_positive_n": sensitivity_denominator,
        "pcr_negative_n": specificity_denominator,
    }


def compute_pcr_capture_sensitivity_specificity(data: pl.DataFrame, virus_acronyms: dict, capture_positive_nc2_threshold: float = 0.05) -> pl.DataFrame:
    diagnostic_data = data.with_columns(
        (pl.col("ct_value").cast(pl.Float64, strict=False).fill_null(0) > 0).alias("pcr_positive"),
        (pl.col("nc2").fill_null(0) >= capture_positive_nc2_threshold).alias("capture_positive"),
    )
    viruses = sorted(diagnostic_data.select("virus").unique().to_series().to_list())
    rows = [_diagnostic_summary_row(diagnostic_data, virus) for virus in viruses]
    rows.append(_diagnostic_summary_row(diagnostic_data))
    return (
        pl.DataFrame(rows)
        .with_columns(pl.col("virus").replace({**virus_acronyms, "overall": "overall"}).alias("virus_acronym"))
        .select(["virus", "virus_acronym", "true_positive", "false_negative", "false_positive", "true_negative", "pcr_positive_n", "pcr_negative_n", "sensitivity", "specificity"])
        .rename(
            {
                "true_positive": "COSMOS_pos_PCR_pos",
                "false_negative": "COSMOS_neg_PCR_pos",
                "false_positive": "COSMOS_pos_PCR_neg",
                "true_negative": "COSMOS_neg_PCR_neg",
            }
        )
    )


def compute_virus_positive_summary_table(data: pl.DataFrame, nc2_groups: list, nc2_labels: list[str], virus_acronyms: dict, capture_positive_nc2_threshold: float = 0.05) -> pl.DataFrame:
    pcr_positive_data = filter_pcr_positive(data)
    category_counts = compute_genome_coverage_category_counts(pcr_positive_data, nc2_groups, nc2_labels, virus_acronyms)
    diagnostic_counts = compute_pcr_capture_sensitivity_specificity(data, virus_acronyms, capture_positive_nc2_threshold)
    return (
        diagnostic_counts.filter(pl.col("virus") != "overall")
        .select(["virus", "virus_acronym", "pcr_positive_n", "COSMOS_pos_PCR_pos", "COSMOS_neg_PCR_pos", "sensitivity"])
        .filter(pl.col("pcr_positive_n") > 0)
        .rename(
            {
                "pcr_positive_n": "PCR positive samples",
                "COSMOS_pos_PCR_pos": f"capture positive samples (NC2 >= {capture_positive_nc2_threshold:g})",
                "COSMOS_neg_PCR_pos": f"capture missed samples (NC2 < {capture_positive_nc2_threshold:g})",
                "sensitivity": "capture sensitivity",
            }
        )
        .join(category_counts, on=["virus", "virus_acronym"], how="left")
        .sort("virus_acronym")
    )


def compute_ct_summary_by_genome_coverage_category(data: pl.DataFrame, nc2_groups: list, nc2_labels: list[str]) -> pl.DataFrame:
    pcr_positive_data = filter_pcr_positive(data)
    category_data = add_genome_coverage_category(pcr_positive_data, nc2_groups, nc2_labels)
    category_order = pl.DataFrame({"genome_coverage_category": nc2_labels, "category_order": list(range(len(nc2_labels)))})
    return (
        category_data.group_by("genome_coverage_category")
        .agg(
            pl.len().alias("n"),
            pl.col("ct_value").mean().alias("mean_ct"),
            pl.col("ct_value").std().alias("sd_ct"),
            pl.col("ct_value").median().alias("median_ct"),
            pl.col("ct_value").min().alias("min_ct"),
            pl.col("ct_value").max().alias("max_ct"),
        )
        .join(category_order, on="genome_coverage_category", how="right")
        .drop("category_order")
    )


def compute_virus_codetection_matrix(data: pl.DataFrame, virus_acronyms: dict, capture_positive_nc2_threshold: float = 0.05) -> pl.DataFrame:
    codetection_data = data.with_columns(
        (pl.col("ct_value").cast(pl.Float64, strict=False).fill_null(0) > 0).alias("pcr_positive"),
        (pl.col("nc2").cast(pl.Float64, strict=False).fill_null(0) >= capture_positive_nc2_threshold).alias("capture_positive"),
        pl.col("virus").replace(virus_acronyms).alias("virus_acronym"),
    ).with_columns(
        pl.when(pl.col("pcr_positive") & pl.col("capture_positive"))
        .then(pl.lit("PCR and COSMOS"))
        .when(pl.col("pcr_positive") & ~pl.col("capture_positive"))
        .then(pl.lit("PCR only"))
        .when(~pl.col("pcr_positive") & pl.col("capture_positive"))
        .then(pl.lit("COSMOS only"))
        .otherwise(pl.lit("neither"))
        .alias("detection_category")
    )
    detected_samples = (
        codetection_data.group_by("sample")
        .agg(
            (pl.col("pcr_positive") | pl.col("capture_positive")).sum().alias("n_detected_viruses"),
            (pl.col("pcr_positive") & pl.col("capture_positive")).sum().alias("n_pcr_and_cosmos"),
            pl.col("capture_positive").sum().alias("n_cosmos_detected"),
            pl.col("pcr_positive").sum().alias("n_pcr_detected"),
        )
        .filter(pl.col("n_detected_viruses") > 0)
    )
    detected_viruses = (
        codetection_data.group_by("virus_acronym")
        .agg((pl.col("pcr_positive") | pl.col("capture_positive")).sum().alias("n_detected_samples"))
        .filter(pl.col("n_detected_samples") > 0)
        .sort(["n_detected_samples", "virus_acronym"], descending=[True, False])
        .get_column("virus_acronym")
        .to_list()
    )
    matrix = (
        codetection_data.join(detected_samples, on="sample", how="inner")
        .filter(pl.col("virus_acronym").is_in(detected_viruses))
        .pivot(index="sample", columns="virus_acronym", values="detection_category")
        .fill_null("neither")
        .join(detected_samples, on="sample", how="left")
        .sort(
            ["n_detected_viruses", "n_pcr_and_cosmos", "n_cosmos_detected", "n_pcr_detected", "sample"],
            descending=[True, True, True, True, False],
        )
    )
    return matrix.select(["sample"] + detected_viruses)

def generate_matrix(data: pl.DataFrame, virus_acronyms: dict, val_column: str = "ct_value", nan_value=0) -> pl.DataFrame:
    # make pivot table of virus and ct_value
    ct_matrix = (
        data.with_columns(pl.col("virus").replace(virus_acronyms).alias("virus_acronym"))
        .pivot(index="sample", columns="virus_acronym", values=val_column)
        .fill_null(nan_value)
    )
    return ct_matrix

def write_virus_codetection_matrix(data: pl.DataFrame, output_file: str, virus_acronyms: dict, capture_positive_nc2_threshold: float = 0.05) -> None:
    compute_virus_codetection_matrix(data, virus_acronyms, capture_positive_nc2_threshold).write_csv(output_file, separator="\t")

def write_virus_count_tables(
    data: pl.DataFrame,
    output_dir: str,
    virus_acronyms: dict,
    nc2_groups: list,
    nc2_labels: list[str],
    capture_positive_nc2_threshold: float = 0.05,
) -> None:
    dedup_read_positive_data = filter_dedup_read_positive(data)
    genome_coverage_counts = compute_genome_coverage_category_counts(dedup_read_positive_data, nc2_groups, nc2_labels, virus_acronyms)
    diagnostic_counts = compute_pcr_capture_sensitivity_specificity(data, virus_acronyms, capture_positive_nc2_threshold=capture_positive_nc2_threshold)
    virus_positive_summary = compute_virus_positive_summary_table(data, nc2_groups, nc2_labels, virus_acronyms, capture_positive_nc2_threshold)
    ct_summary = compute_ct_summary_by_genome_coverage_category(data, nc2_groups, nc2_labels)
    genome_coverage_counts.write_csv(f"{output_dir}/virus_genome_coverage_category_counts.tsv", separator="\t")
    diagnostic_counts.write_csv(f"{output_dir}/pcr_vs_capture_sensitivity_specificity.tsv", separator="\t")
    virus_positive_summary.write_csv(f"{output_dir}/virus_positive_summary_by_virus.tsv", separator="\t")
    ct_summary.write_csv(f"{output_dir}/ct_summary_by_genome_coverage_category.tsv", separator="\t")

def get_coverage_across_bam(input_bam: str, output_tsv: str) -> pl.DataFrame:
    if not os.path.exists(output_tsv):
        compute_coverage_at_each_position_using_samtools(input_bam, output_tsv)
    else:
        print(f"Coverage file {output_tsv} already exists. Skipping coverage computation.")
    df = pl.read_csv(output_tsv, separator="\t", has_header=False)
    df = df.rename({"column_1": "chrom", "column_2": "position", "column_3": "depth"})
    return df

def get_virus_coverage_from_castanet_outputs(virus_name: str, castanet_nc2_file: str, sample_config: str, virus_acronyms: dict, detection_threshold: float = 0.5) -> pl.DataFrame:
    # find samples with the virus detected in the castanet nc2 file
    # read the castanet nc2 file
    nc2_df = pl.read_csv(castanet_nc2_file, separator="\t")
    # convert cols to floats
    for col in nc2_df.columns:
        if col != "sampleid":
            nc2_df = nc2_df.with_columns(pl.col(col).cast(pl.Float64))
    # get sampleid where virus is detected (i.e. nc2 > 0)
    detected_samples = set(nc2_df.filter(pl.col(virus_name) > detection_threshold)["sampleid"].to_list())
    # parse sample config to get file path
    sample_config = parse_json(sample_config, remove_nest=True)
    depth_results = []
    for sample in detected_samples:
        sample_entry = find_sample_entry_in_multikey_dict(sample_config, sample)
        fastq_file = sample_entry['R1']
        fastq_dir = os.path.dirname(fastq_file)
        head_dir = os.path.dirname(fastq_dir)
        castanet_sample_folder = os.path.join(
            head_dir, "castanet", "v9.2", "trimmed", f"{sample}"
        )
        castanet_depth_folder = os.path.join(castanet_sample_folder, "Depth_output")

        if virus_name == "umbrevirus":
            # Castanet's standard depth file contains only one LSV segment.
            # Extract all three contigs from the sample BAM and concatenate
            # them in this stable order instead.
            lsv_segments = {
                "segment_1": "UmbreVirus_262256411645004674050341512477535811351",
                "segment_2": "UmbreVirus_217043466430343676518769285178128316378",
                "segment_3": "UmbreVirus_200112682654376413090417675282467369124",
            }
            castanet_output_bam = find_castanet_bam_with_lsv_contigs(
                castanet_sample_folder,
                lsv_segments,
            )
            lsv_tmp_dir = os.path.join(
                "tmp", "virus_lsv_coverage", sample
            )
            depth_df = get_depth_across_lsv_and_concat(
                castanet_output_bam,
                lsv_segments,
                sample,
                lsv_tmp_dir,
            )
            depth_results.append(
                depth_df.select(["sample", "virus", "position", "depth"])
            )
            continue

        depth_by_pos_file = os.path.join(castanet_depth_folder, f"{virus_name}-{sample}_depth_by_pos.csv")
        if os.path.exists(depth_by_pos_file):
            depth_df = pl.read_csv(depth_by_pos_file, separator=",", has_header=False) # parse csv
            depth_df = depth_df.transpose(include_header=False) # transpose
            depth_df.columns = ["depth", "dedup_depth"] # rename cols
            depth_df = depth_df.with_columns(pl.lit(sample).alias("sample"), pl.lit(virus_name).alias("virus")) # add sample and virus cols
            depth_df = depth_df.with_columns(pl.Series("position", range(1, depth_df.height + 1))) # add position col
            depth_results.append(depth_df)
        else:
            print(f"Depth file not found for sample {sample} and virus {virus_name}: {depth_by_pos_file}")
    if depth_results:
        combined_depth_df = pl.concat(depth_results)
        return combined_depth_df.select(["sample", "virus", "position", "depth"])
    else:
        return pl.DataFrame({"sample": [], "virus": [], "position": [], "depth": []})

def prepare_virus_coverage_to_plot(virus_coverage_df: pl.DataFrame, virus_acronyms: dict) -> pl.DataFrame:
    # replace virus names with acronyms
    virus_coverage_df = virus_coverage_df.with_columns(pl.col("virus").replace(virus_acronyms).alias("virus_acronym"))
    # convert position to int and depth to float
    virus_coverage_df = virus_coverage_df.with_columns(
        pl.col("position").cast(pl.Int64),
        pl.col("depth").cast(pl.Float64),
    )
    # average across samples for each virus and position
    virus_coverage_avg_df = virus_coverage_df.group_by(["virus_acronym", "position"]).agg(
        pl.col("depth").mean().alias("mean_depth"),
        pl.col("depth").std(ddof=1).alias("depth_sd"),
        pl.col("sample").n_unique().alias("n_samples"),
    )
    # calculate the normal-approximation 95% confidence interval half-width
    virus_coverage_avg_df = virus_coverage_avg_df.with_columns(
        pl.when(pl.col("n_samples") >= 2)
        .then(1.96 * pl.col("depth_sd") / pl.col("n_samples").sqrt())
        .otherwise(None)
        .alias("depth_ci95")
    )
    virus_coverage_avg_df = virus_coverage_avg_df.with_columns(
        (pl.col("mean_depth") - pl.col("depth_ci95"))
        .clip(lower_bound=0)
        .alias("depth_ci95_lower"),
        (pl.col("mean_depth") + pl.col("depth_ci95")).alias("depth_ci95_upper"),
    )
    # median, q25 and q75
    virus_coverage_avg_df = virus_coverage_df.group_by(["virus_acronym", "position"]).agg(
        pl.col("depth").median().alias("median_depth"),
        pl.col("depth").quantile(0.25).alias("depth_q25"),
        pl.col("depth").quantile(0.75).alias("depth_q75"),
    ).join(virus_coverage_avg_df, on=["virus_acronym", "position"], how="left")
    return virus_coverage_avg_df.sort(["virus_acronym", "position"])

def build_virus_coverage_plotting_data(virus_list: list[str], castanet_nc2_file: str, sample_config: str, virus_acronyms: dict, detection_threshold: float = 0.5, metadata_file: str = None) -> pl.DataFrame:
    coverage_results = []
    for virus in virus_list:
        coverage_df = get_virus_coverage_from_castanet_outputs(virus, castanet_nc2_file, sample_config, virus_acronyms, detection_threshold)
        coverage_df = coverage_df.with_columns(
            pl.col("sample").cast(pl.String),
            pl.col("virus").cast(pl.String),
            pl.col("position").cast(pl.Int64),
            pl.col("depth").cast(pl.Float64, strict=False),
        )
        coverage_results.append(coverage_df)
    # coverage_results = pl.concat(coverage_results)
    # coverage_results = coverage_results.with_columns(pl.col("virus").replace(virus_acronyms).alias("virus_acronym"))
    # return coverage_results
    if coverage_results:
        combined_coverage_df = pl.concat(coverage_results)
        if metadata_file is not None:
            metadata_df = pl.read_csv(metadata_file, separator="\t")
            combined_coverage_df = combined_coverage_df.join(metadata_df.select("sample"), on="sample", how="inner") # filter to samples in metadata
        combined_coverage_df.select(
            "sample",
            "virus",
            "position",
            pl.col("depth").alias("coverage"),
        ).write_csv(
            f"{RESULTS_DIR}/virus/data/virus_coverage_per_sample.tsv",
            separator="\t",
        )
        prepared_coverage_df = prepare_virus_coverage_to_plot(combined_coverage_df, virus_acronyms)
        return prepared_coverage_df
    else:
        return pl.DataFrame({"virus_acronym": [], "position": [], "mean_depth": [], "depth_ci95_lower": [], "depth_ci95_upper": [], "depth_sd": [], "median_depth": [], "depth_q25": [], "depth_q75": [], "n_samples": []})

def find_lsv_contigs_in_castanet_bam(input_bam: str, lsv_segments: dict[str, str]) -> dict[str, str]:
    result = subprocess.run(
        ["samtools", "idxstats", input_bam],
        check=True,
        capture_output=True,
        text=True,
    )
    bam_contigs = [
        line.split("\t")[0]
        for line in result.stdout.splitlines()
        if line
    ]
    matched_contigs = {}
    for segment, accession in lsv_segments.items():
        matches = [
            contig
            for contig in bam_contigs
            if accession in contig
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one BAM contig containing {accession}, "
                f"but found {matches}"
            )
        matched_contigs[segment] = matches[0]
    return matched_contigs

def get_depth_across_lsv_and_concat(castanet_output_bam: str, lsv_segments: dict[str, str], sample: str, tmp_dir: str) -> pl.DataFrame:
    os.makedirs(tmp_dir, exist_ok=True)
    lsv_contigs = find_lsv_contigs_in_castanet_bam(castanet_output_bam, lsv_segments)
    depth_results = []
    position_offset = 0
    for segment, contig in lsv_contigs.items():
        # filter bam
        filtered_bam = os.path.join(tmp_dir, f"{sample}_{segment}.bam")
        subprocess.run(
            ["samtools", "view", "-b", "-o", filtered_bam, castanet_output_bam, contig],
            check=True,
        )
        depth_file = os.path.join(tmp_dir, f"{sample}_{segment}_depth.tsv")
        depth_df = get_coverage_across_bam(filtered_bam, depth_file)
        depth_df = depth_df.filter(pl.col("chrom") == contig)
        depth_df = depth_df.with_columns(
            pl.col("position").cast(pl.Int64).alias("segment_position"),
            (pl.col("position").cast(pl.Int64) + position_offset).alias("position"),
            pl.col("depth").cast(pl.Float64),
            pl.lit(segment).alias("segment"),
            pl.lit(sample).alias("sample"),
            pl.lit("umbrevirus").alias("virus"),
        )
        depth_results.append(depth_df)
        position_offset += depth_df.height
    if depth_results:
        return pl.concat(depth_results).select(
            ["sample", "virus", "segment", "segment_position", "position", "depth"]
        )
    return pl.DataFrame(
        schema={
            "sample": pl.String,
            "virus": pl.String,
            "segment": pl.String,
            "segment_position": pl.Int64,
            "position": pl.Int64,
            "depth": pl.Float64,
        }
    )


def find_castanet_bam_with_lsv_contigs(
    castanet_sample_folder: str,
    lsv_segments: dict[str, str],
) -> str:
    matching_bams = []
    for bam in Path(castanet_sample_folder).glob("*.bam"):
        try:
            find_lsv_contigs_in_castanet_bam(str(bam), lsv_segments)
        except (ValueError, subprocess.CalledProcessError):
            continue
        matching_bams.append(str(bam))

    if len(matching_bams) != 1:
        raise ValueError(
            f"Expected one BAM containing all LSV contigs in "
            f"{castanet_sample_folder}, but found {matching_bams}"
        )
    return matching_bams[0]
