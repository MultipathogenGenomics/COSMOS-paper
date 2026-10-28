"""Shared Cx. annulirostris map membership and site averages (no plotting imports)."""
from pathlib import Path

import polars as pl

ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR = 0.0
AGG_SITE_FOR_AVERAGE_BUFFER = 0.01
ANNU_MIN_COUNT_THRESHOLD = 20
SITE_AVERAGES_FILE = "results/publication_figures/tables/cx_annulirostris_site_averages.tsv"


def annu_abundance_filter():
    abundance = pl.col("Culex annulirostris proportion")
    threshold = ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR
    return abundance.is_finite() & ((abundance >= 0) if threshold == 0 else (abundance > threshold))


def add_annu_site_clusters(df, buffer=AGG_SITE_FOR_AVERAGE_BUFFER):
    """Bin coordinates on the figure's rounded grid, not by pairwise distance."""
    return df.with_columns(pl.concat_str(
        [(pl.col("longitude") / buffer).round().cast(pl.Int64),
         (pl.col("latitude") / buffer).round().cast(pl.Int64)],
        separator="_",
    ).alias("agg_site"))


def average_annu_sites(points):
    """Average each plotted sample equally; retain collated ratio precision.

    IDs are deterministic for a given set of groups, ordered west to east then
    south to north. agg_site records the actual grid key for traceability.
    """
    return (add_annu_site_clusters(points).group_by("agg_site").agg(
        pl.col("latitude").mean().alias("site_latitude"),
        pl.col("longitude").mean().alias("site_longitude"),
        pl.col("avg_subtype_ratio").mean().alias("mean_saus_ratio"),
        pl.col("avg_subtype_ratio").std().alias("std_saus_ratio"),
        pl.col("Culex annulirostris proportion").mean().alias("mean_annu_abundance"),
        pl.col("Culex annulirostris proportion").std().alias("std_annu_abundance"),
        pl.len().alias("n_samples"),
        pl.col("sample").sort().str.join("; ").alias("samples"),
        pl.col("location").drop_nulls().unique().sort().str.join("; ").alias("locations"),
        pl.col("location-site").drop_nulls().unique().sort().str.join("; ").alias("site_names"),
        pl.col("sample").filter(pl.col("MVEV") > 0.05).sort().str.join("; ").alias("mvev_samples"),
    ).sort("site_longitude", "site_latitude", "agg_site")
      .with_row_index("site_id", offset=1))


def write_annu_site_averages(sites, output_file=SITE_AVERAGES_FILE):
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    sites.write_csv(output_file, separator="\t")


def load_annu_plot_points(
    metadata_file="data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv",
    abundance_file="results/publication_figures/data/combined_abundance_data.tsv",
    virus_file="data/final_castanet_data/COSMOS1-4.final_castanet_nc2.tsv",
):
    """Reproduce the map's point selection from text inputs, without rendering."""
    ratios = pl.concat([
        pl.read_csv(f"data/cryptic_subpops/COSMOS{i}.Cx_annulirostris_collated_results.tsv", separator="\t")
        .select("sample", "min_count", "avg_subtype_ratio")
        for i in range(1, 5)
    ])
    metadata = pl.read_csv(metadata_file, separator="\t").select(
        "sample", "latitude", "longitude", "location", "location-site")
    abundance = pl.read_csv(abundance_file).select(
        "sample", pl.col("Culex annulirostris").alias("Culex annulirostris proportion"))
    virus = pl.read_csv(virus_file, separator="\t").select(
        pl.col("sampleid").alias("sample"), pl.col("murrayvalleyencephalitisvirus").alias("MVEV"))
    return (ratios.filter(pl.col("min_count") > ANNU_MIN_COUNT_THRESHOLD)
            .join(metadata, on="sample", how="inner", validate="1:1")
            .join(abundance, on="sample", how="left", validate="1:1")
            .join(virus, on="sample", how="left", validate="1:1")
            .with_columns(pl.col("latitude", "longitude", "avg_subtype_ratio").fill_nan(None))
            .drop_nulls(["latitude", "longitude", "avg_subtype_ratio"])
            .filter(annu_abundance_filter()))
