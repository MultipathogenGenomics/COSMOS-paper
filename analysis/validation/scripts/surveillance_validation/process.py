from datetime import datetime

import numpy as np
import polars as pl

from analysis.common.abundance import read_abundance_for_sample_sets


def normalise_location_expr(col: str) -> pl.Expr:
    return pl.col(col).cast(pl.Utf8).str.strip_chars().str.to_lowercase().str.replace_all(r"\s+", " ")


def parse_nsw_health_surveillance_counts(
    nsw_counts_file: str,
    target_species: dict[str, str],
) -> pl.DataFrame:
    species_to_keep = list(target_species.values()) + ["weekly_totals"]
    counts_df = pl.read_csv(nsw_counts_file)
    counts_df = counts_df.with_columns(
        pl.col("count_avg").cast(pl.Float64, strict=False).alias("count_avg"),
        pl.col("trap_date").str.strptime(pl.Date, "%d/%m/%Y", strict=False).alias("trap_date_parsed"),
        normalise_location_expr("location").alias("location_key"),
    )
    counts_df = counts_df.filter(
        pl.col("species").is_in(species_to_keep)
        & pl.col("count_avg").is_not_null()
        & pl.col("trap_date_parsed").is_not_null()
    )
    counts_df = counts_df.group_by(
        [
            "location",
            "location_key",
            "epi_week",
            "trap_date",
            "trap_date_parsed",
            "week_ending_saturday",
            "site_type",
            "lat",
            "lon",
            "species",
        ]
    ).agg(pl.col("count_avg").mean().alias("count_avg"))
    counts_wide = counts_df.pivot(
        index=[
            "location",
            "location_key",
            "epi_week",
            "trap_date",
            "trap_date_parsed",
            "week_ending_saturday",
            "site_type",
            "lat",
            "lon",
        ],
        columns="species",
        values="count_avg",
    )

    for species in target_species.values():
        if species not in counts_wide.columns:
            counts_wide = counts_wide.with_columns(pl.lit(0.0).alias(species))
    if "weekly_totals" not in counts_wide.columns:
        counts_wide = counts_wide.with_columns(pl.lit(None, dtype=pl.Float64).alias("weekly_totals"))

    counts_wide = counts_wide.rename(
        {
            "weekly_totals": "weekly_total",
            "aedes_vigilax": "aedes_vigilax_count",
            "culex_annulirostris": "culex_annulirostris_count",
        }
    )
    for species_slug in target_species.values():
        count_col = f"{species_slug}_count"
        if count_col not in counts_wide.columns and species_slug in counts_wide.columns:
            counts_wide = counts_wide.rename({species_slug: count_col})

    proportion_exprs = []
    for species_slug in target_species.values():
        count_col = f"{species_slug}_count"
        proportion_col = f"{species_slug}_proportion"
        proportion_exprs.append(
            pl.when(pl.col("weekly_total") > 0)
            .then(pl.col(count_col) / pl.col("weekly_total"))
            .otherwise(None)
            .alias(proportion_col)
        )
    return counts_wide.with_columns(proportion_exprs)


def parse_cosmos_trap_metadata(metadata_file: str) -> pl.DataFrame:
    return (
        pl.read_csv(metadata_file, separator="\t")
        .with_columns(
            pl.col("date-set").str.strptime(pl.Date, "%Y-%m-%d", strict=False).alias("trap_date_parsed"),
            normalise_location_expr("location").alias("location_key"),
        )
    )


def read_sample_read_assignment_summaries(estimate_abundance_dir: str, sample_sets: list[str]) -> pl.DataFrame:
    summary_dfs = []
    for sample_set in sample_sets:
        summary_file = f"{estimate_abundance_dir}/{sample_set}/sample_read_assignment_summary.tsv"
        summary_df = pl.read_csv(summary_file, separator="\t")
        summary_dfs.append(summary_df)
    return pl.concat(summary_dfs, how="diagonal")


def read_filter_pass_read_counts(estimate_abundance_dir: str, sample_sets: list[str]) -> pl.DataFrame:
    read_count_dfs = []
    for sample_set in sample_sets:
        read_count_file = f"{estimate_abundance_dir}/{sample_set}/filter_pass_read_counts.tsv"
        read_count_df = pl.read_csv(read_count_file, separator="\t")
        read_count_df = read_count_df.with_columns(pl.lit(sample_set).alias("batch"))
        read_count_dfs.append(read_count_df)
    read_count_df = pl.concat(read_count_dfs, how="diagonal")
    read_count_cols = [col for col in read_count_df.columns if col not in ["sample", "batch"]]
    read_count_df = read_count_df.with_columns([pl.col(col).fill_null(0) for col in read_count_cols])
    return read_count_df.select(["sample", "batch"] + read_count_cols)


def prepare_surveillance_validation_data(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    abundance_filename: str,
    cosmos_metadata_file: str,
    nsw_counts_file: str,
    target_species: dict[str, str],
    assignment_reference: str | None = None,
    disagreement_threshold: float = 0.25,
) -> dict:
    abundance_df = read_abundance_for_sample_sets(
        estimate_abundance_dir,
        sample_sets,
        abundance_filename,
    )
    metadata_df = parse_cosmos_trap_metadata(cosmos_metadata_file)
    nsw_counts_df = parse_nsw_health_surveillance_counts(nsw_counts_file, target_species)
    read_assignment_summary_df = read_sample_read_assignment_summaries(estimate_abundance_dir, sample_sets)
    filter_pass_read_counts_df = read_filter_pass_read_counts(estimate_abundance_dir, sample_sets)

    abundance_cols = ["sample", "batch"] + [species for species in target_species if species in abundance_df.columns]
    missing_species = [species for species in target_species if species not in abundance_df.columns]
    if missing_species:
        raise ValueError(f"Missing target species in abundance data: {missing_species}")

    cosmos_df = metadata_df.join(abundance_df.select(abundance_cols), on="sample", how="left")
    matched_df = cosmos_df.join(
        nsw_counts_df,
        on=["location_key", "trap_date_parsed"],
        how="inner",
        suffix="_nsw",
    )

    validation_rows = []
    for species_name, species_slug in target_species.items():
        validation_rows.append(
            matched_df.select(
                [
                    "sample",
                    "batch",
                    "location",
                    "location-site",
                    "location2",
                    "date-set",
                    "trap_date",
                    "epi_week",
                    pl.col("mosquitoes").alias("n_mosquitoes"),
                    "weekly_total",
                    pl.lit(species_name).alias("species"),
                    pl.col(species_name).alias("estimated_abundance"),
                    pl.col(f"{species_slug}_count").alias("surveillance_count"),
                    pl.col(f"{species_slug}_proportion").alias("surveillance_proportion"),
                ]
            )
        )
    validation_df = pl.concat(validation_rows)
    validation_df = validation_df.filter(
        pl.col("estimated_abundance").is_not_null()
        & pl.col("surveillance_proportion").is_not_null()
    )
    correlation_df = compute_surveillance_validation_correlations(validation_df)
    target_read_count_rows = []
    for species_name in target_species:
        target_read_count_rows.append(
            filter_pass_read_counts_df.select(
                [
                    "sample",
                    pl.lit(species_name).alias("species"),
                    pl.col(species_name).alias("target_filter_pass_reads"),
                ]
            )
        )
    target_read_count_df = pl.concat(target_read_count_rows)
    disagreement_df = (
        validation_df.with_columns(
            (pl.col("estimated_abundance") - pl.col("surveillance_proportion")).abs().alias("absolute_difference"),
            (pl.col("estimated_abundance") - pl.col("surveillance_proportion")).alias("signed_difference"),
        )
        .join(
            read_assignment_summary_df.select(
                [
                    "sample",
                    "kallisto:n_processed",
                    "kallisto:n_assigned",
                    "kallisto:proportion_assigned",
                    "kma:n_processed",
                    "kma:n_assigned",
                    "kma:proportion_assigned",
                ]
            ),
            on="sample",
            how="left",
        )
        .join(target_read_count_df, on=["sample", "species"], how="left")
        .filter(pl.col("absolute_difference") >= disagreement_threshold)
        .sort("absolute_difference", descending=True)
    )

    unmatched_cosmos_df = cosmos_df.join(
        nsw_counts_df.select(["location_key", "trap_date_parsed"]).unique(),
        on=["location_key", "trap_date_parsed"],
        how="anti",
    )
    unmatched_nsw_df = nsw_counts_df.join(
        cosmos_df.select(["location_key", "trap_date_parsed"]).unique(),
        on=["location_key", "trap_date_parsed"],
        how="anti",
    )

    return {
        "abundance_df": abundance_df,
        "read_assignment_summary_df": read_assignment_summary_df,
        "filter_pass_read_counts_df": filter_pass_read_counts_df,
        "cosmos_metadata_df": metadata_df,
        "nsw_counts_df": nsw_counts_df,
        "matched_df": matched_df,
        "validation_df": validation_df,
        "correlation_df": correlation_df,
        "disagreement_df": disagreement_df,
        "unmatched_cosmos_df": unmatched_cosmos_df,
        "unmatched_nsw_df": unmatched_nsw_df,
        "target_species": target_species,
        "disagreement_threshold": disagreement_threshold,
        "sample_sets": sample_sets,
        "assignment_reference": assignment_reference,
        "estimate_abundance_dir": estimate_abundance_dir,
        "abundance_filename": abundance_filename,
        "cosmos_metadata_file": cosmos_metadata_file,
        "nsw_counts_file": nsw_counts_file,
    }


def write_surveillance_validation_tables(validation_data: dict, output_dir: str) -> None:
    validation_data["validation_df"].write_csv(f"{output_dir}/surveillance_validation_data.tsv", separator="\t")
    validation_data["correlation_df"].write_csv(f"{output_dir}/surveillance_validation_correlations.tsv", separator="\t")
    validation_data["disagreement_df"].write_csv(f"{output_dir}/surveillance_validation_disagreements.tsv", separator="\t")
    validation_data["matched_df"].write_csv(f"{output_dir}/matched_cosmos_nsw_traps.tsv", separator="\t")
    validation_data["unmatched_cosmos_df"].write_csv(f"{output_dir}/unmatched_cosmos_traps.tsv", separator="\t")
    validation_data["unmatched_nsw_df"].write_csv(f"{output_dir}/unmatched_nsw_collections.tsv", separator="\t")


def rank_values(values: np.ndarray) -> np.ndarray:
    sorter = np.argsort(values)
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[sorter]
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[sorter[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return ranks


def compute_correlation_values(x_values: np.ndarray, y_values: np.ndarray) -> dict:
    x_mean = float(np.mean(x_values)) if len(x_values) > 0 else None
    y_mean = float(np.mean(y_values)) if len(y_values) > 0 else None
    x_variance = float(np.var(x_values, ddof=1)) if len(x_values) > 1 else None
    y_variance = float(np.var(y_values, ddof=1)) if len(y_values) > 1 else None
    x_y_covariance = float(np.cov(x_values, y_values, ddof=1)[0, 1]) if len(x_values) > 1 else None
    concordance_correlation_coefficient = None
    if (
        x_mean is not None
        and y_mean is not None
        and x_variance is not None
        and y_variance is not None
        and x_y_covariance is not None
    ):
        ccc_denominator = x_variance + y_variance + (x_mean - y_mean) ** 2
        if ccc_denominator != 0:
            concordance_correlation_coefficient = (2 * x_y_covariance) / ccc_denominator
    if len(x_values) < 3 or np.std(x_values) == 0 or np.std(y_values) == 0:
        return {
            "x_mean": x_mean,
            "y_mean": y_mean,
            "x_variance": x_variance,
            "y_variance": y_variance,
            "concordance_correlation_coefficient": concordance_correlation_coefficient,
            "pearson_r": None,
            "pearson_r_squared": None,
            "pearson_p_value": None,
            "spearman_rho": None,
            "spearman_p_value": None,
        }
    try:
        from scipy import stats

        pearson_result = stats.pearsonr(x_values, y_values)
        spearman_result = stats.spearmanr(x_values, y_values)
        pearson_r = float(pearson_result.statistic)
        pearson_p_value = float(pearson_result.pvalue)
        spearman_rho = float(spearman_result.statistic)
        spearman_p_value = float(spearman_result.pvalue)
    except ImportError:
        pearson_r = float(np.corrcoef(x_values, y_values)[0, 1])
        pearson_p_value = None
        spearman_rho = float(np.corrcoef(rank_values(x_values), rank_values(y_values))[0, 1])
        spearman_p_value = None
    return {
        "x_mean": x_mean,
        "y_mean": y_mean,
        "x_variance": x_variance,
        "y_variance": y_variance,
        "concordance_correlation_coefficient": concordance_correlation_coefficient,
        "pearson_r": pearson_r,
        "pearson_r_squared": pearson_r ** 2,
        "pearson_p_value": pearson_p_value,
        "spearman_rho": spearman_rho,
        "spearman_p_value": spearman_p_value,
    }


def compute_surveillance_validation_correlations(validation_df: pl.DataFrame) -> pl.DataFrame:
    correlation_rows = []
    for species in validation_df.select("species").unique().sort("species").to_series().to_list():
        species_df = validation_df.filter(pl.col("species") == species).select(
            ["estimated_abundance", "surveillance_proportion"]
        ).drop_nulls()
        x_values = species_df["surveillance_proportion"].to_numpy()
        y_values = species_df["estimated_abundance"].to_numpy()
        correlation_values = compute_correlation_values(x_values, y_values)
        correlation_rows.append(
            {
                "species": species,
                "x": "NSW Health surveillance proportion",
                "y": "COSMOS estimated abundance",
                "n": len(x_values),
                **correlation_values,
            }
        )
    return pl.DataFrame(correlation_rows)


def write_surveillance_validation_description(
    validation_data: dict,
    output_file: str,
    plot_output_file: str,
) -> None:
    validation_df = validation_data["validation_df"]
    matched_df = validation_data["matched_df"]
    with open(output_file, "w") as handle:
        handle.write(f"created_at\t{datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"plot_output_file\t{plot_output_file}\n")
        if validation_data["assignment_reference"] is not None:
            handle.write(f"assignment_reference\t{validation_data['assignment_reference']}\n")
        handle.write(f"estimate_abundance_dir\t{validation_data['estimate_abundance_dir']}\n")
        handle.write(f"sample_sets\t{'; '.join(validation_data['sample_sets'])}\n")
        handle.write(f"abundance_filename\t{validation_data['abundance_filename']}\n")
        handle.write(f"cosmos_metadata_file\t{validation_data['cosmos_metadata_file']}\n")
        handle.write(f"nsw_counts_file\t{validation_data['nsw_counts_file']}\n")
        handle.write("join_keys\tCOSMOS location + date-set to NSW Health location + trap_date\n")
        handle.write("surveillance_total\tweekly_totals row from NSW Health species column\n")
        handle.write(f"disagreement_threshold\t{validation_data['disagreement_threshold']}\n")
        handle.write(f"n_cosmos_metadata_traps\t{validation_data['cosmos_metadata_df'].height}\n")
        handle.write(f"n_abundance_samples\t{validation_data['abundance_df'].height}\n")
        handle.write(f"n_nsw_collection_rows\t{validation_data['nsw_counts_df'].height}\n")
        handle.write(f"n_matched_cosmos_traps\t{matched_df.height}\n")
        handle.write(f"n_validation_rows\t{validation_df.height}\n")
        handle.write(f"correlation_file\tsurveillance_validation_correlations.tsv\n")
        handle.write(f"n_disagreement_rows\t{validation_data['disagreement_df'].height}\n")
        handle.write(f"n_unmatched_cosmos_traps\t{validation_data['unmatched_cosmos_df'].height}\n")
        handle.write(f"n_unmatched_nsw_collections\t{validation_data['unmatched_nsw_df'].height}\n")
        for species in validation_data["target_species"]:
            n_species = validation_df.filter(pl.col("species") == species).height
            handle.write(f"n_validation_rows_{species.replace(' ', '_')}\t{n_species}\n")
