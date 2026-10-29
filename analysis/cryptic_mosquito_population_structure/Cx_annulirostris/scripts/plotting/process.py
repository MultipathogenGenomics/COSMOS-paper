# Prepare data for cryptic subpopulation maps.

import geopandas as gpd
import pandas as pd
import os

def load_metadata(metadata_file: str) -> pd.DataFrame:
    return pd.read_csv(metadata_file, sep="\t")


def load_nsw_act_boundary(geopkg_file: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(geopkg_file)
    nsw_act = gdf[gdf["STATE_NAME_2021"].isin(["New South Wales", "Australian Capital Territory"])]
    return nsw_act.dissolve(by="STATE_NAME_2021")


def collate_subpopulation_results(file_template: str, datasets: list[str]) -> pd.DataFrame:
    dfs = []
    for dataset in datasets:
        df_file = file_template.format(dataset=dataset)
        dfs.append(pd.read_csv(df_file, sep="\t").assign(dataset=dataset))
    return pd.concat(dfs, ignore_index=True)


def prepare_cx_annulirostris_coastal_association(metadata_file: str, file_template: str, datasets: list[str]) -> pd.DataFrame:
    """Reproduce the coastal analysis independently of map coordinate filters."""
    subtype_df = collate_subpopulation_results(file_template, datasets)
    metadata = load_metadata(metadata_file)[["sample", "distance_to_coast_km"]]
    joined = subtype_df.merge(metadata, on="sample", how="inner", validate="many_to_one")
    for column in ("min_count", "avg_subtype_ratio", "distance_to_coast_km"):
        joined[column] = pd.to_numeric(joined[column], errors="coerce")
    return joined.loc[joined["min_count"] > 20].replace(
        [float("inf"), -float("inf")], float("nan")
    ).dropna(subset=["avg_subtype_ratio", "distance_to_coast_km"])


def derive_ae_notoscriptus_population_results(config: dict, datasets: list[str]) -> pd.DataFrame:
    one_vs_23 = collate_subpopulation_results(config["one_vs_23_file_template"], datasets)
    two_vs_3 = collate_subpopulation_results(config["two_vs_3_file_template"], datasets)
    one_vs_23["proportion_cluster1"] = pd.to_numeric(one_vs_23["avg_subtype_ratio"], errors="coerce")
    two_vs_3["proportion_cluster2_given_23"] = pd.to_numeric(two_vs_3["avg_subtype_ratio"], errors="coerce")
    merged = one_vs_23[
        [
            "sample",
            "dataset",
            "proportion_cluster1",
            "min_count",
            "n_snps_which_meet_threshold",
        ]
    ].merge(
        two_vs_3[
            [
                "sample",
                "dataset",
                "proportion_cluster2_given_23",
                "min_count",
                "n_snps_which_meet_threshold",
            ]
        ],
        on=["sample", "dataset"],
        how="inner",
        suffixes=("_1_vs_23", "_2_vs_3"),
    )
    merged = merged[
        (merged["n_snps_which_meet_threshold_1_vs_23"] >= config["one_vs_23_min_snps"])
        & (merged["n_snps_which_meet_threshold_2_vs_3"] >= config["two_vs_3_min_snps"])
    ].copy()
    merged["proportion_cluster2"] = (1 - merged["proportion_cluster1"]) * merged["proportion_cluster2_given_23"]
    merged["proportion_cluster3"] = (1 - merged["proportion_cluster1"]) * (1 - merged["proportion_cluster2_given_23"])
    return merged


def load_panel_subpopulation_results(config: dict, datasets: list[str]) -> pd.DataFrame:
    if config.get("source") == "ae_notoscriptus_pairwise":
        return derive_ae_notoscriptus_population_results(config, datasets)
    return collate_subpopulation_results(config["file_template"], datasets)


def prepare_subpopulation_points(
    subtype_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    map_boundary: gpd.GeoDataFrame,
    ratio_columns: list[str],
    min_count_threshold: int | None = None,
    complement_ratio_column: str | None = None,
    complement_output_column: str | None = None,
) -> gpd.GeoDataFrame:
    if min_count_threshold is not None:
        subtype_df = subtype_df[subtype_df["min_count"] > min_count_threshold]

    subtype_df = subtype_df.join(metadata_df.set_index("sample"), on="sample", how="left")
    if complement_ratio_column is not None and complement_output_column is not None:
        subtype_df[complement_output_column] = 1 - subtype_df[complement_ratio_column]
    subtype_df = subtype_df.dropna(subset=["latitude", "longitude"] + ratio_columns)
    points = gpd.GeoDataFrame(
        subtype_df.copy(),
        geometry=gpd.points_from_xy(subtype_df["longitude"], subtype_df["latitude"]),
        crs="EPSG:4326",
    ).to_crs(map_boundary.crs)
    return points


def load_estimated_abundance_for_sample_sets(
    estimate_abundance_dir: str,
    sample_sets: list[str],
    abundance_filename: str,
) -> pd.DataFrame:
    dfs = []
    for sample_set in sample_sets:
        abundance_file = f"{estimate_abundance_dir}/{sample_set}/{abundance_filename}"
        dfs.append(pd.read_csv(abundance_file, sep="\t").assign(batch=sample_set))
    return pd.concat(dfs, ignore_index=True)


def compute_detection_summary(
    metadata_df: pd.DataFrame,
    abundance_df: pd.DataFrame,
    panels: list[dict],
) -> pd.DataFrame:
    summary_rows = []
    metadata_samples = metadata_df[["sample"]].drop_duplicates()
    for panel in panels:
        species_col = panel["species_abundance_column"]
        abundance_joined = metadata_samples.merge(abundance_df[["sample", species_col]], on="sample", how="left")
        n_species_detected = (pd.to_numeric(abundance_joined[species_col], errors="coerce").fillna(0) > 0).sum()
        points = panel["points"]
        summary_rows.append(
            {
                "species": panel["label"],
                "n_metadata_traps": metadata_samples["sample"].nunique(),
                "n_species_detected_by_abundance": int(n_species_detected),
                "n_subpopulation_estimates": points["sample"].nunique(),
                "filters": panel.get("summary_filters", ""),
            }
        )
    return pd.DataFrame(summary_rows)


def write_subpopulation_detection_summary(
    subpopulation_data: dict,
    output_file: str,
    estimate_abundance_dir: str,
    sample_sets: list[str],
    abundance_filename: str,
) -> None:
    abundance_df = load_estimated_abundance_for_sample_sets(
        estimate_abundance_dir,
        sample_sets,
        abundance_filename,
    )
    summary = compute_detection_summary(
        subpopulation_data["metadata"],
        abundance_df,
        subpopulation_data["panels"],
    )
    summary.to_csv(output_file, sep="\t", index=False)


def prepare_cryptic_subpopulation_data(
    geopkg_file: str,
    metadata_file: str,
    panel_configs: list[dict],
    datasets: list[str],
) -> dict:
    metadata_df = load_metadata(metadata_file)
    map_boundary = load_nsw_act_boundary(geopkg_file)
    panels = []
    for config in panel_configs:
        subtype_df = load_panel_subpopulation_results(config, datasets)
        points = prepare_subpopulation_points(
            subtype_df,
            metadata_df,
            map_boundary,
            config["ratio_columns"],
            min_count_threshold=config.get("min_count_threshold"),
            complement_ratio_column=config.get("complement_ratio_column"),
            complement_output_column=config.get("complement_output_column"),
        )
        panels.append({**config, "subtype_df": subtype_df, "points": points})
    return {
        "map_boundary": map_boundary,
        "metadata": metadata_df,
        "panels": panels,
    }

def combine_Cx_annulirostris_subpopulation_kmer_count_data(
    results_dir: str,
    samples_to_keep: list[str] | None = None,
    subtype_fraction_threshold: float = 0.8,
) -> pd.DataFrame:
    kmer_count_files = sorted([f for f in os.listdir(results_dir) if f.endswith(".csv")])
    concated_df = None
    for f in kmer_count_files:
        df = pd.read_csv(os.path.join(results_dir, f))
        df = df[df["subtype_fraction"] >= subtype_fraction_threshold]
        df = df.fillna(0)
        if concated_df is None:
            concated_df = df
        else:
            concated_df = pd.concat([concated_df, df], ignore_index=True)
    if samples_to_keep is not None:
        concated_df = concated_df[concated_df["sample"].isin(samples_to_keep)]
    subtype_order = ["S-AUS", "AUS", "PNG-1", "PNG-2", "SI"]
    concated_df["subtype"] = concated_df["subtype"].str.replace("ann ", "", regex=False)
    subtype_kmer_counts = concated_df.groupby("subtype")["kmer"].nunique().to_dict()
    concated_df = concated_df.pivot_table(index="sample", columns="subtype", values="count", aggfunc="sum", fill_value=0)
    if samples_to_keep is not None:
        concated_df = concated_df.reindex(samples_to_keep, fill_value=0)
    concated_df = concated_df.reindex(columns=subtype_order, fill_value=0)
    concated_df.attrs["subtype_kmer_counts"] = {subtype: subtype_kmer_counts.get(subtype, 0) for subtype in subtype_order}
    concated_df = concated_df.assign(_s_aus_total=concated_df["S-AUS"])
    concated_df = concated_df.sort_values("_s_aus_total", ascending=False).drop(columns="_s_aus_total")
    return concated_df
