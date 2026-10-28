from pathlib import Path

import numpy as np
import polars as pl


def read_abundance_for_sample_sets(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    abundance_filename: str = "abundance_within_filter_pass_reads.tsv",
) -> pl.DataFrame:
    abundance_dfs = []
    for sample_set in sample_sets:
        abundance_file = Path(estimate_abundance_dir) / sample_set / abundance_filename
        abundance_df = pl.read_csv(abundance_file, separator="\t")
        abundance_df = abundance_df.with_columns(pl.lit(sample_set).alias("batch"))
        abundance_dfs.append(abundance_df)

    combined_df = pl.concat(abundance_dfs, how="diagonal")
    abundance_cols = [col for col in combined_df.columns if col not in ["sample", "batch"]]
    combined_df = combined_df.with_columns([pl.col(col).fill_null(0) for col in abundance_cols])
    return combined_df.select(["sample", "batch"] + abundance_cols)


def normalise_abundance_species_name(species_name: str) -> str:
    normalised_name = species_name.strip().replace("_", " ")
    if normalised_name.startswith("Ochlerotatus "):
        normalised_name = normalised_name.replace("Ochlerotatus ", "Aedes ", 1)
    return normalised_name


def get_abundance_value(row: dict, abundance_column_lookup: dict[str, str], species_name: str) -> float:
    abundance_column = abundance_column_lookup.get(normalise_abundance_species_name(species_name).lower())
    if abundance_column is None:
        return 0.0
    value = row.get(abundance_column)
    return 0.0 if value is None else float(value)


def compute_mixture_abundance_error_records(
    abundance_df: pl.DataFrame,
    mixture_metadata_df: pl.DataFrame,
    tool: str,
    relative_error_margin: float = 0.10,
) -> pl.DataFrame:
    abundance_column_lookup = {
        normalise_abundance_species_name(column).lower(): column
        for column in abundance_df.columns
        if column not in ["sample", "batch"]
    }
    n_components = len([column for column in mixture_metadata_df.columns if column.startswith("proportion_")])
    mixture_df = (
        abundance_df
        .with_columns(pl.col("sample").str.extract(r"mixture_([0-9]+)", 1).cast(pl.Int64).alias("mixture_id"))
        .join(mixture_metadata_df, on="mixture_id", how="inner")
    )

    rows = []
    for row in mixture_df.to_dicts():
        mixture_records = []
        for index in range(1, n_components + 1):
            expected_abundance = float(row[f"proportion_{index}"])
            species = normalise_abundance_species_name(str(row[f"sample_{index}"]))
            estimated_abundance = get_abundance_value(row, abundance_column_lookup, species)
            error = estimated_abundance - expected_abundance
            relative_error = error / expected_abundance if expected_abundance > 0 else None
            mixture_records.append(
                {
                    "tool": tool,
                    "mixture_id": row["mixture_id"],
                    "species": species,
                    "expected_abundance": expected_abundance,
                    "estimated_abundance": estimated_abundance,
                    "error_estimated_minus_expected": error,
                    "error_expected_minus_estimated": expected_abundance - estimated_abundance,
                    "absolute_error": abs(error),
                    "squared_error": error**2,
                    "relative_error": relative_error,
                    "absolute_relative_error": abs(relative_error) if relative_error is not None else None,
                    "fold_estimated_over_expected": estimated_abundance / expected_abundance if expected_abundance > 0 else None,
                    "within_relative_error_margin": abs(relative_error) <= relative_error_margin if relative_error is not None else None,
                    "relative_error_margin": relative_error_margin,
                }
            )
        mixture_rmse = float(np.sqrt(np.mean([record["squared_error"] for record in mixture_records])))
        for record in mixture_records:
            record["mixture_rmse"] = mixture_rmse
            rows.append(record)
    return pl.DataFrame(rows)


def summarise_mixture_abundance_errors(error_df: pl.DataFrame) -> pl.DataFrame:
    group_cols = ["tool", "species"] if "tool" in error_df.columns else ["species"]

    def summarise(df: pl.DataFrame, extra_cols: dict) -> dict:
        summary = df.select(
            pl.len().alias("n_estimates"),
            pl.col("mixture_id").n_unique().alias("n_mixtures"),
            pl.col("expected_abundance").mean().alias("mean_expected_abundance"),
            pl.col("estimated_abundance").mean().alias("mean_estimated_abundance"),
            pl.col("estimated_abundance").median().alias("median_estimated_abundance"),
            pl.col("error_estimated_minus_expected").mean().alias("mean_error_estimated_minus_expected"),
            pl.col("error_estimated_minus_expected").median().alias("median_error_estimated_minus_expected"),
            pl.col("absolute_error").mean().alias("mean_absolute_error"),
            pl.col("absolute_error").median().alias("median_absolute_error"),
            pl.col("absolute_relative_error").mean().alias("mean_absolute_relative_error"),
            pl.col("absolute_relative_error").median().alias("median_absolute_relative_error"),
            pl.col("mixture_rmse").mean().alias("mean_mixture_rmse"),
            pl.col("within_relative_error_margin").sum().alias("n_within_relative_error_margin"),
            (~pl.col("within_relative_error_margin")).sum().alias("n_outside_relative_error_margin"),
            pl.col("relative_error_margin").first().alias("relative_error_margin"),
        ).row(0, named=True)
        summary["rmse"] = float(np.sqrt(df["squared_error"].mean()))
        summary["proportion_within_relative_error_margin"] = (
            summary["n_within_relative_error_margin"] / summary["n_estimates"]
            if summary["n_estimates"] > 0
            else None
        )
        return {**extra_cols, **summary}

    rows = []
    if "tool" in error_df.columns:
        for tool in error_df["tool"].unique().to_list():
            tool_df = error_df.filter(pl.col("tool") == tool)
            rows.append(summarise(tool_df, {"tool": tool, "summary_level": "overall", "species": "all"}))
    else:
        rows.append(summarise(error_df, {"summary_level": "overall", "species": "all"}))

    for group_values, group_df in error_df.group_by(group_cols):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        extra_cols = dict(zip(group_cols, group_values))
        extra_cols["summary_level"] = "species"
        rows.append(summarise(group_df, extra_cols))

    preferred_cols = [
        "tool",
        "summary_level",
        "species",
        "n_estimates",
        "n_mixtures",
        "mean_expected_abundance",
        "mean_estimated_abundance",
        "median_estimated_abundance",
        "rmse",
        "mean_mixture_rmse",
        "mean_error_estimated_minus_expected",
        "median_error_estimated_minus_expected",
        "mean_absolute_error",
        "median_absolute_error",
        "mean_absolute_relative_error",
        "median_absolute_relative_error",
        "relative_error_margin",
        "n_within_relative_error_margin",
        "n_outside_relative_error_margin",
        "proportion_within_relative_error_margin",
    ]
    output_df = pl.DataFrame(rows)
    return output_df.select([col for col in preferred_cols if col in output_df.columns])
