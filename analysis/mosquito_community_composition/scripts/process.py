# Prepare data for abundance geospatial plots.

from datetime import datetime

import numpy as np
import polars as pl
from scipy import stats
from sklearn.decomposition import PCA
from umap import UMAP

from ...common.abundance import read_abundance_for_sample_sets


COASTAL_THRESHOLD_KM = 20
URBAN_CODES = [1, 2, 3, 4, 5]
FRINGE_CODES = [6, 7]
RURAL_CODES = [8, 9, 10, 11]
URBAN_RURAL_HEATMAP_ORDER = ["Rural", "Met_Fringe", "Urban"]


def normalise_species_name(species_name: str) -> str:
    normalised_name = species_name.strip()
    if normalised_name.startswith("Ochlerotatus "):
        normalised_name = normalised_name.replace("Ochlerotatus ", "Aedes ", 1)
    return normalised_name


def load_abundance_and_metadata(abundance_file: str, metadata_file: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    abundance_df = pl.read_csv(abundance_file, separator="\t")
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    return abundance_df, metadata_df


def load_estimated_abundance_and_metadata(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    metadata_file: str,
    abundance_filename: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    abundance_df = read_abundance_for_sample_sets(
        estimate_abundance_dir=estimate_abundance_dir,
        sample_sets=sample_sets,
        abundance_filename=abundance_filename,
    )
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    return abundance_df, metadata_df


def filter_to_most_common(df: pl.DataFrame, top_n: int = 10) -> pl.DataFrame:
    species_cols = [col for col in df.columns if col not in ["sample", "batch"]]
    species_sums = df.select([pl.col(col).sum().alias(col) for col in species_cols]).row(0, named=True)
    top_species = sorted(species_sums, key=species_sums.get, reverse=True)[:top_n]
    keep_cols = ["sample"] + (["batch"] if "batch" in df.columns else [])
    return df.select(keep_cols + top_species)



def centre_log_ratio_transform(data: np.ndarray, pseudo_count: float = 1e-9) -> np.ndarray:
    data = data.copy()
    data += pseudo_count
    log_data = np.log(data)
    geometric_means = np.exp(np.mean(log_data, axis=1)).reshape(-1, 1)
    clr_data = log_data - np.log(geometric_means)
    return clr_data


def add_geospatial_categories(metadata_df: pl.DataFrame) -> pl.DataFrame:
    return metadata_df.with_columns(
        pl.when(pl.col("distance_to_coast_km") < COASTAL_THRESHOLD_KM)
        .then(pl.lit(1))
        .when(pl.col("distance_to_coast_km") >= COASTAL_THRESHOLD_KM)
        .then(pl.lit(0))
        .otherwise(pl.lit(None))
        .alias("Coastal"),
        pl.when(pl.col("council_olg_code").is_in(URBAN_CODES))
        .then(pl.lit("Urban"))
        .when(pl.col("council_olg_code").is_in(FRINGE_CODES))
        .then(pl.lit("Met_Fringe"))
        .when(pl.col("council_olg_code").is_in(RURAL_CODES))
        .then(pl.lit("Rural"))
        .otherwise(pl.lit(None))
        .alias("urban_rural"),
    )


def get_abundance_matrix(df: pl.DataFrame, metadata_cols: list[str]) -> np.ndarray:
    return df.select(pl.exclude(["sample", "batch"] + metadata_cols)).to_numpy()


def compute_umap_embedding(abundance_df: pl.DataFrame, metadata_df: pl.DataFrame, top_n: int = 15) -> pl.DataFrame:
    from umap import UMAP

    plot_df = filter_to_most_common(abundance_df, top_n=top_n)
    plot_df = plot_df.join(metadata_df.select(["sample", "Coastal", "urban_rural", "season"]), on="sample", how="inner")
    print("Running UMAP on abundance data for samples with metadata, N =", plot_df.shape[0])

    umap_data = plot_df.select(pl.exclude(["sample", "batch", "Coastal", "urban_rural", "season"])).to_numpy()
    umap_data_transformed = centre_log_ratio_transform(umap_data)
    embedding = UMAP(random_state=42).fit_transform(umap_data_transformed)

    embedding_df = pl.DataFrame(
        {
            "UMAP1": embedding[:, 0],
            "UMAP2": embedding[:, 1],
            "sample": plot_df["sample"],
        }
    )
    return embedding_df.join(metadata_df.select(["sample", "Coastal", "urban_rural"]), on="sample", how="left")


def get_top_species(abundance_df: pl.DataFrame, top_n: int) -> list[str]:
    species_cols = [col for col in abundance_df.columns if col not in ["sample", "batch"]]
    species_sums = abundance_df.select([pl.col(col).sum().alias(col) for col in species_cols]).row(0, named=True)
    return sorted(species_sums, key=species_sums.get, reverse=True)[:top_n]


def compute_species_trap_counts(abundance_df: pl.DataFrame, metadata_df: pl.DataFrame) -> pl.DataFrame:
    species_cols = [
        col
        for col in abundance_df.columns
        if col not in ["sample", "batch", "n_non_confidently_assigned_reads"]
    ]
    abundance_df = (
        metadata_df.select("sample")
        .unique()
        .join(abundance_df.select(["sample"] + species_cols), on="sample", how="left")
    )
    n_total_traps = abundance_df.height
    rows = []
    for species_col in species_cols:
        species_abundance = abundance_df.select(
            pl.col(species_col).cast(pl.Float64, strict=False).fill_null(0).alias("abundance")
        )
        present_abundance = species_abundance.filter(pl.col("abundance") > 0)
        n_traps_present = present_abundance.height
        rows.append(
            {
                "species": normalise_species_name(species_col),
                "abundance_column": species_col,
                "n_traps_present": n_traps_present,
                "n_total_traps": n_total_traps,
                "proportion_traps_present": n_traps_present / n_total_traps if n_total_traps > 0 else None,
                "mean_abundance_in_present_traps": present_abundance.select(pl.col("abundance").mean()).item() if n_traps_present > 0 else None,
                "median_abundance_in_present_traps": present_abundance.select(pl.col("abundance").median()).item() if n_traps_present > 0 else None,
            }
        )
    return pl.DataFrame(rows).sort(["n_traps_present", "species"], descending=[True, False])


def write_species_trap_counts(species_trap_counts: pl.DataFrame, output_file: str) -> None:
    species_trap_counts.write_csv(output_file, separator="\t")


def write_simpsons_diversity_index(sdi_df: pl.DataFrame, output_file: str) -> None:
    sdi_df.write_csv(output_file, separator="\t")


def compute_simpsons_diversity_index_by_council_type(sdi_df: pl.DataFrame) -> pl.DataFrame:
    return (
        sdi_df.group_by("council_type")
        .agg(
            pl.len().alias("n_traps"),
            pl.col("Simpsons_Diversity_Index").mean().alias("mean_simpsons_diversity_index"),
            pl.col("Simpsons_Diversity_Index").median().alias("median_simpsons_diversity_index"),
            pl.col("Simpsons_Diversity_Index").std().alias("sd_simpsons_diversity_index"),
            pl.col("Simpsons_Diversity_Index").quantile(0.25).alias("q1_simpsons_diversity_index"),
            pl.col("Simpsons_Diversity_Index").quantile(0.75).alias("q3_simpsons_diversity_index"),
        )
        .with_columns(
            (pl.col("q3_simpsons_diversity_index") - pl.col("q1_simpsons_diversity_index")).alias("iqr_simpsons_diversity_index")
        )
        .sort("mean_simpsons_diversity_index", descending=True)
    )


def compute_pca_embedding_and_loadings(
    abundance_df: pl.DataFrame,
    metadata_df: pl.DataFrame,
    top_n: int = 15,
    n_components: int = 5,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    plot_df = filter_to_most_common(abundance_df, top_n=top_n)
    plot_df = plot_df.join(metadata_df.select(["sample", "Coastal", "urban_rural", "season"]), on="sample", how="inner")
    species = [col for col in plot_df.columns if col not in ["sample", "batch", "Coastal", "urban_rural", "season"]]
    pca_data = plot_df.select(species).to_numpy()
    pca_data_transformed = centre_log_ratio_transform(pca_data)
    n_components = min(n_components, len(species), plot_df.height)
    pca = PCA(n_components=n_components, random_state=42)
    pca_embedding = pca.fit_transform(pca_data_transformed)
    embedding_df = pl.DataFrame(
        {
            "PCA1": pca_embedding[:, 0],
            "PCA2": pca_embedding[:, 1],
            "sample": plot_df["sample"],
        }
    )
    loadings_df = pl.DataFrame({"species": species})
    loadings_df = loadings_df.with_columns(
        [
            pl.Series(f"PC{component_index + 1}", pca.components_[component_index, :])
            for component_index in range(n_components)
        ]
    )
    variance_df = pl.DataFrame(
        {
            "principal_component": [f"PC{component_index + 1}" for component_index in range(n_components)],
            "component_number": list(range(1, n_components + 1)),
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_explained_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        }
    )
    return (
        embedding_df.join(
            metadata_df.select(["sample", "Coastal", "urban_rural", "distance_to_coast_km"]),
            on="sample",
            how="left",
        ),
        loadings_df,
        variance_df,
    )


def prepare_embedding_comparison_data(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    metadata_file: str,
    abundance_filename: str,
    virus_nc2_file: str,
    virus_acronyms: dict[str, str],
    s_aus_file_template: str,
    embedding_top_n: int = 15,
    s_aus_min_count: int = 20,
) -> pl.DataFrame:
    abundance_df, metadata_df = load_estimated_abundance_and_metadata(
        estimate_abundance_dir,
        sample_sets,
        metadata_file,
        abundance_filename,
    )
    metadata_df = add_geospatial_categories(metadata_df)
    plot_df = filter_to_most_common(abundance_df, top_n=embedding_top_n).join(
        metadata_df.select("sample"),
        on="sample",
        how="inner",
    )
    species = [column for column in plot_df.columns if column not in ["sample", "batch"]]
    transformed_abundance = centre_log_ratio_transform(plot_df.select(species).to_numpy())

    pca_embedding = PCA(n_components=2, random_state=42).fit_transform(transformed_abundance)

    umap_embedding = UMAP(random_state=42).fit_transform(transformed_abundance)
    embedding_df = pl.DataFrame(
        {
            "sample": plot_df["sample"],
            "PCA1": pca_embedding[:, 0],
            "PCA2": pca_embedding[:, 1],
            "UMAP1": umap_embedding[:, 0],
            "UMAP2": umap_embedding[:, 1],
        }
    )

    nc2_df = pl.read_csv(virus_nc2_file, separator="\t")
    missing_virus_columns = [column for column in virus_acronyms.keys() if column not in nc2_df.columns]
    if missing_virus_columns:
        raise ValueError(f"Virus nc2 columns not found: {missing_virus_columns}")
    # if nc2 > 0.05 OR ct > 0 then is a virus positive
    virus_status_df = nc2_df.select(["sampleid"] + list(virus_acronyms.keys())).join(
        metadata_df.select(["sample"] + [f"ct_{short_virus}" for short_virus in virus_acronyms.values() if f"ct_{short_virus}" in metadata_df.columns]),
        left_on="sampleid",
        right_on="sample",
    )
    nc2_columns = [col for col in virus_acronyms.keys() if col in virus_status_df.columns]
    ct_columns = [f"ct_{short_virus}" for short_virus in virus_acronyms.values() if f"ct_{short_virus}" in virus_status_df.columns]
    virus_status_df = virus_status_df.with_columns(
        (
            pl.any_horizontal(
                pl.col(col).cast(pl.Float64, strict=False).fill_null(0) > 0.05
                for col in nc2_columns
            )
            | pl.any_horizontal(
                pl.col(col).cast(pl.Float64, strict=False).fill_null(0) > 0
                for col in ct_columns
            )
        )
        .cast(pl.Int8)
        .alias("virus_positive")
    ).rename({"sampleid": "sample"})

    # virus_status_df = (
    #     nc2_df.select(["sampleid"] + virus_columns)
    #     .with_columns(
    #         pl.max_horizontal(
    #             [pl.col(column).cast(pl.Float64, strict=False).fill_null(0) for column in virus_columns]
    #         ).alias("maximum_virus_nc2")
    #     )
    #     .select(
    #         pl.col("sampleid").alias("sample"),
    #         (pl.col("maximum_virus_nc2") > 0.05).alias("virus_positive"),
    #     )
    # )

    s_aus_dfs = [
        pl.read_csv(s_aus_file_template.format(dataset=dataset), separator="\t").select(
            "sample",
            pl.col("min_count").cast(pl.Float64, strict=False),
            pl.col("avg_subtype_ratio").cast(pl.Float64, strict=False).alias("s_aus_proportion"),
        )
        for dataset in sample_sets
    ]
    s_aus_df = (
        pl.concat(s_aus_dfs)
        .filter(
            (pl.col("min_count") > s_aus_min_count)
            & pl.col("s_aus_proportion").is_not_null()
        )
        .join(
            metadata_df.select(["sample", "latitude", "longitude"]),
            on="sample",
            how="inner",
        )
        .filter(pl.col("latitude").is_not_null() & pl.col("longitude").is_not_null())
        .group_by("sample")
        .agg(pl.col("s_aus_proportion").first())
    )

    diversity_df = compute_sdi_by_sample(abundance_df, metadata_df).select(
        "sample",
        pl.col("Simpsons_Diversity_Index").alias("diversity_1_minus_D"),
    )

    # add season 
    season_df = metadata_df.select(["sample", "season"])

    return (
        embedding_df.join(
            metadata_df.select(
                "sample",
                "Coastal",
                pl.when(pl.col("council_olg_code").is_in([1, 2, 3]))
                .then(pl.lit("Met"))
                .when(pl.col("council_olg_code").is_in([4, 5]))
                .then(pl.lit("Urban"))
                .when(pl.col("council_olg_code").is_in([6, 7]))
                .then(pl.lit("Met-fringe"))
                .when(pl.col("council_olg_code").is_in([8, 9, 10, 11]))
                .then(pl.lit("Rural"))
                .otherwise(pl.lit(None))
                .alias("region_category"),
            ),
            on="sample",
            how="left",
        )
        .join(virus_status_df, on="sample", how="left")
        .join(s_aus_df, on="sample", how="left")
        .join(diversity_df, on="sample", how="left")
        .join(season_df, on="sample", how="left")
        .with_columns(pl.col("virus_positive").fill_null(False))
    )


def compute_pc_distance_to_coast_associations(embedding_df: pl.DataFrame) -> pl.DataFrame:
    pc_cols = [col for col in embedding_df.columns if col.startswith("PCA")]
    rows = []
    for pc_col in pc_cols:
        model_df = embedding_df.select(
            [
                pl.col(pc_col).cast(pl.Float64, strict=False),
                pl.col("distance_to_coast_km").cast(pl.Float64, strict=False),
            ]
        ).drop_nulls()
        if model_df.height < 3:
            rows.append(
                {
                    "principal_component": pc_col.replace("PCA", "PC"),
                    "n": model_df.height,
                    "model_status": "not_fit_too_few_samples",
                }
            )
            continue
        result = stats.linregress(
            model_df["distance_to_coast_km"].to_numpy(),
            model_df[pc_col].to_numpy(),
        )
        rows.append(
            {
                "principal_component": pc_col.replace("PCA", "PC"),
                "n": model_df.height,
                "model": f"{pc_col} ~ distance_to_coast_km",
                "model_status": "fit",
                "slope_pc_units_per_km": result.slope,
                "intercept": result.intercept,
                "r_value": result.rvalue,
                "r_squared": result.rvalue**2,
                "p_value": result.pvalue,
                "slope_standard_error": result.stderr,
                "intercept_standard_error": result.intercept_stderr,
            }
        )
    return pl.DataFrame(rows)


def write_pc_distance_to_coast_associations(association_df: pl.DataFrame, output_file: str) -> None:
    association_df.write_csv(output_file, separator="\t")


def compute_simpsons_diversity_index(proportions: list[float]) -> float:
    assert np.isclose(sum(proportions), 1), f"Proportions should sum to 1, but sum to {sum(proportions)}"
    D = [p**2 for p in proportions]
    return 1 - sum(D)


def compute_sdi_by_sample(abundance_df: pl.DataFrame, metadata_df: pl.DataFrame) -> pl.DataFrame:
    df = abundance_df.join(
        metadata_df.select(["sample", "Coastal", "urban_rural", "season", "council_olg_code", "council_type"]),
        on="sample",
        how="inner",
    )
    print("Computing SDI by council type, N =", df.shape[0])

    metadata_cols = ["batch", "Coastal", "urban_rural", "season", "council_olg_code", "council_type"]
    sdi_values = []
    for row in df.iter_rows(named=True):
        species_proportions = [row[col] for col in df.columns if col not in ["sample"] + metadata_cols]
        sdi_values.append(
            {
                "sample": row["sample"],
                "Coastal": row["Coastal"],
                "urban_rural": row["urban_rural"],
                "season": row["season"],
                "Simpsons_Diversity_Index": compute_simpsons_diversity_index(species_proportions),
                "council_olg_code": row["council_olg_code"],
                "council_type": row["council_type"],
            }
        )
    return pl.DataFrame(sdi_values)


def prepare_heatmap_data(
    abundance_df: pl.DataFrame,
    metadata_df: pl.DataFrame,
    features_to_order_by: list[str],
    custom_mosquito_order: list[str] | None = None,
    min_abundance: float = 0.01,
) -> dict:
    cols_to_heatmap = [col for col in abundance_df.columns if col not in ["sample", "batch", "n_non_confidently_assigned_reads"]]
    joined = abundance_df.join(metadata_df, on="sample", how="inner")
    print("Preparing abundance heatmap, N =", joined.shape[0])
    sort_expressions = [
        (
            pl.col(feature).replace_strict(
                {label: index for index, label in enumerate(URBAN_RURAL_HEATMAP_ORDER)},
                default=len(URBAN_RURAL_HEATMAP_ORDER),
            )
            if feature == "urban_rural"
            else pl.col(feature)
        )
        for feature in features_to_order_by
    ]
    joined = joined.sort(sort_expressions)

    if custom_mosquito_order is not None:
        custom_mosquito_order = [normalise_species_name(species) for species in custom_mosquito_order]
        normalised_column_to_column = {normalise_species_name(col): col for col in cols_to_heatmap}
        cols_to_heatmap = [
            normalised_column_to_column[species]
            for species in custom_mosquito_order
            if species in normalised_column_to_column
        ]
        missing_cols = [species for species in custom_mosquito_order if species not in normalised_column_to_column]
        if missing_cols:
            print(f"Warning: the following columns in the custom mosquito order were not found in the data and will be skipped: {missing_cols}")

    heatmap_data = joined.select(cols_to_heatmap).to_numpy()
    heatmap_data[heatmap_data < min_abundance] = 0
    non_zero_cols = [index for index in range(heatmap_data.shape[1]) if heatmap_data[:, index].sum() > 0]
    heatmap_data = heatmap_data[:, non_zero_cols]
    species_labels = [normalise_species_name(cols_to_heatmap[index]) for index in non_zero_cols]

    heatmap_data = heatmap_data.T
    if custom_mosquito_order is None:
        species_totals = heatmap_data.sum(axis=1)
        species_order = np.argsort(species_totals)[::-1]
        heatmap_data = heatmap_data[species_order, :]
        species_labels = [species_labels[index] for index in species_order]

    ordered_metadata = joined.select(["sample"] + features_to_order_by)
    return {
        "heatmap_data": heatmap_data,
        "species_labels": species_labels,
        "metadata": ordered_metadata,
        "features": features_to_order_by,
    }


def prepare_abundance_geospatial_data(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    metadata_file: str,
    abundance_filename: str = "abundance_within_filter_pass_reads.tsv",
    heatmap_top_n: int = 10,
    embedding_top_n: int = 15,
    embedding_method: str = "UMAP",
    pca_n_components: int = 5,
) -> dict:
    abundance_df, metadata_df = load_estimated_abundance_and_metadata(
        estimate_abundance_dir,
        sample_sets,
        metadata_file,
        abundance_filename,
    )
    metadata_df = add_geospatial_categories(metadata_df)
    heatmap_abundance_df = filter_to_most_common(abundance_df, top_n=heatmap_top_n)
    custom_mosquito_order = [
        "Culex annulirostris",
        "Aedes notoscriptus",
        "Culex australicus",
        "Anopheles annulipes",
        "Culex pipiens",
        "Mansonia uniformis",
        "Aedes vigilax",
        "Culex sitiens",
        "Aedes procax",
        "Coquillettidia linealis",
    ]
    embedding_method = embedding_method.upper()
    if embedding_method == "UMAP":
        embedding_df = compute_umap_embedding(abundance_df, metadata_df, top_n=embedding_top_n)
        pca_loadings_df = None
        pca_variance_df = None
    elif embedding_method == "PCA":
        embedding_df, pca_loadings_df, pca_variance_df = compute_pca_embedding_and_loadings(
            abundance_df,
            metadata_df,
            top_n=embedding_top_n,
            n_components=pca_n_components,
        )
    else:
        raise ValueError("embedding_method must be either 'UMAP' or 'PCA'")

    sdi_df = compute_sdi_by_sample(abundance_df, metadata_df)

    return {
        "heatmap": prepare_heatmap_data(
            heatmap_abundance_df,
            metadata_df,
            features_to_order_by=["Coastal", "urban_rural"],
            custom_mosquito_order=custom_mosquito_order,
        ),
        "embedding": embedding_df,
        "embedding_method": embedding_method,
        "pca_loadings": pca_loadings_df,
        "pca_variance": pca_variance_df,
        "sdi": sdi_df,
        "sdi_by_council_type": compute_simpsons_diversity_index_by_council_type(sdi_df),
        "species_trap_counts": compute_species_trap_counts(abundance_df, metadata_df),
        "pc_distance_to_coast_associations": compute_pc_distance_to_coast_associations(embedding_df)
        if embedding_method == "PCA"
        else None,
        "description": {
            "n_abundance_samples": abundance_df.height,
            "n_metadata_samples": metadata_df.height,
            "n_samples_with_abundance_and_metadata": abundance_df.join(metadata_df.select("sample"), on="sample", how="inner").height,
            "heatmap_top_n": heatmap_top_n,
            "embedding_top_n": embedding_top_n,
            "embedding_method": embedding_method,
            "pca_n_components": pca_n_components,
            "heatmap_top_species": get_top_species(abundance_df, heatmap_top_n),
            "embedding_top_species": get_top_species(abundance_df, embedding_top_n),
            "custom_heatmap_species_order": custom_mosquito_order,
        },
    }


def write_abundance_geospatial_description(
    geospatial_data: dict,
    output_file: str,
    estimate_abundance_dir: str,
    sample_sets: list[str],
    metadata_file: str,
    abundance_filename: str,
    output_plot_file: str,
) -> None:
    description = geospatial_data["description"]
    heatmap_species = geospatial_data["heatmap"]["species_labels"]
    embedding_species = description["embedding_top_species"]
    sdi_species = [
        col
        for col in read_abundance_for_sample_sets(
            estimate_abundance_dir,
            sample_sets,
            abundance_filename,
        ).columns
        if col not in ["sample", "batch"]
    ]
    lines = [
        f"date\t{datetime.now().isoformat(timespec='seconds')}",
        f"output_plot_file\t{output_plot_file}",
        f"estimate_abundance_dir\t{estimate_abundance_dir}",
        f"sample_sets\t{','.join(sample_sets)}",
        f"abundance_filename\t{abundance_filename}",
        f"metadata_file\t{metadata_file}",
        f"n_abundance_samples\t{description['n_abundance_samples']}",
        f"n_metadata_samples\t{description['n_metadata_samples']}",
        f"n_samples_with_abundance_and_metadata\t{description['n_samples_with_abundance_and_metadata']}",
        f"n_heatmap_samples\t{geospatial_data['heatmap']['metadata'].height}",
        f"n_embedding_samples\t{geospatial_data['embedding'].height}",
        f"n_sdi_samples\t{geospatial_data['sdi'].height}",
        f"heatmap_top_n\t{description['heatmap_top_n']}",
        f"embedding_top_n\t{description['embedding_top_n']}",
        f"embedding_method\t{description['embedding_method']}",
        f"pca_n_components\t{description['pca_n_components']}",
        f"n_heatmap_species_plotted\t{len(heatmap_species)}",
        f"heatmap_species_plotted\t{'; '.join(heatmap_species)}",
        f"n_heatmap_top_species_before_custom_order\t{len(description['heatmap_top_species'])}",
        f"heatmap_top_species_before_custom_order\t{'; '.join(description['heatmap_top_species'])}",
        f"custom_heatmap_species_order\t{'; '.join(description['custom_heatmap_species_order'])}",
        f"n_embedding_species\t{len(embedding_species)}",
        f"embedding_species\t{'; '.join(embedding_species)}",
        f"n_sdi_species\t{len(sdi_species)}",
        f"sdi_species\t{'; '.join(sdi_species)}",
        f"coastal_threshold_km\t{COASTAL_THRESHOLD_KM}",
        f"urban_codes\t{URBAN_CODES}",
        f"fringe_codes\t{FRINGE_CODES}",
        f"rural_codes\t{RURAL_CODES}",
    ]
    with open(output_file, "w") as handle:
        handle.write("\n".join(lines) + "\n")
