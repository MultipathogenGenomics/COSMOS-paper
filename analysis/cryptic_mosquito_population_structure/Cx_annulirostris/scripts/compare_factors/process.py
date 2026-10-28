import polars as pl
from analysis.Cx_annulirostris.compare_factors.plots import scatter_snp_ratios_and_spatial_factors, virus_swarmplot_and_ct_scatter, scatter_lat_lon_coloured_by_snp_ratios, scatter_lat_lon_coloured_by_snp_ratio_split_by_season

def compare_subtype_snp_ratio_with_spatial_factors(snp_ratio_files: list, spatial_factors_file: str, spatial_factors_list: list, output_file: str, subtype_col: str="avg_subtype_ratio") -> None:
    # concat snp ratio dfs
    snp_dfs = [pl.read_csv(f, separator='\t') for f in snp_ratio_files]
    snp_df = pl.concat(snp_dfs)
    # combine with spatial factors
    spatial_df = pl.read_csv(spatial_factors_file, separator='\t')
    combined_df = snp_df.join(spatial_df, on="sample", how="left")
    # make plot
    scatter_snp_ratios_and_spatial_factors(combined_df, spatial_factors_list, subtype_col, output_file)

def generate_per_sample_spatial_factor_table(location_to_latlon_file: str, sampleid_to_pt_file:str, pt_to_location_file: str, pt_to_dist_to_coast_file: str) -> pl.DataFrame:
    # parse files - csv
    location_to_latlon_df = pl.read_csv(location_to_latlon_file, encoding="latin-1", ignore_errors=True)
    sampleid_to_pt_df = pl.read_csv(sampleid_to_pt_file, encoding="latin-1", ignore_errors=True)
    pt_to_location_df = pl.read_csv(pt_to_location_file, encoding="latin-1", ignore_errors=True)
    pt_to_dist_to_coast_df = pl.read_csv(pt_to_dist_to_coast_file, encoding="latin-1", ignore_errors=True)
    # cast pt to int for consistent joining
    sampleid_to_pt_df = sampleid_to_pt_df.with_columns(pl.col("pt").cast(pl.String))
    pt_to_location_df = pt_to_location_df.with_columns(pl.col("pt").cast(pl.String))
    pt_to_dist_to_coast_df = pt_to_dist_to_coast_df.with_columns(pl.col("pt").cast(pl.String))
    # combine pt to latlon via location
    pt_to_latlon_df = pt_to_location_df.join(location_to_latlon_df, on=["location-site", "location2"], how="left")
    # combine to sample id via pt
    sampleid_to_latlon_df = sampleid_to_pt_df.join(pt_to_latlon_df, on="pt", how="left")
    # combine with distance to coast
    sampleid_to_latlon_df = sampleid_to_latlon_df.join(pt_to_dist_to_coast_df, on="pt", how="left")
    # filter to sample id, pt, lat, lon, distance to coast, location, location-site
    sampleid_to_latlon_df = sampleid_to_latlon_df.rename({"sampleid": "sample"})
    sampleid_to_latlon_df = sampleid_to_latlon_df.select(["sample", "pt", "lat", "lon", "distance_from_coast", "location", "location-site"])
    return sampleid_to_latlon_df

def generate_per_sample_temporal_factor_table(sampleid_to_pt_file:str, pt_to_season_file: str) -> pl.DataFrame:
    # parse files - csv
    sampleid_to_pt_df = pl.read_csv(sampleid_to_pt_file, encoding="latin-1", ignore_errors=True)
    pt_to_season_df = pl.read_csv(pt_to_season_file, encoding="latin-1", ignore_errors=True)
    # cast pt to int for consistent joining
    sampleid_to_pt_df = sampleid_to_pt_df.with_columns(pl.col("pt").cast(pl.String))
    pt_to_season_df = pt_to_season_df.with_columns(pl.col("pt").cast(pl.String))
    # combine pt to season
    pt_to_season_df = pt_to_season_df.select(["pt", "season"])
    # combine to sample id via pt
    sampleid_to_season_df = sampleid_to_pt_df.join(pt_to_season_df, on="pt", how="left")
    # filter to sample id, pt, season
    sampleid_to_season_df = sampleid_to_season_df.rename({"sampleid": "sample"})
    sampleid_to_season_df = sampleid_to_season_df.select(["sample", "pt", "season"])
    return sampleid_to_season_df

def generate_per_sample_virus_file(sampleid_to_pt_file: str, pt_to_virus_file: str, virus_List: list=["RRV", "BFV", "MVEV", "KUN", "KOK", "STRV", "EHV", "JEV"]) -> pl.DataFrame:
    sampleid_to_pt_df = pl.read_csv(sampleid_to_pt_file, encoding="latin-1", ignore_errors=True)
    pt_to_virus_df = pl.read_csv(pt_to_virus_file, encoding="latin-1", ignore_errors=True)
    pt_to_virus_df = pt_to_virus_df.with_columns(pl.col("pt").cast(pl.String))
    sampleid_to_virus_df = sampleid_to_pt_df.join(pt_to_virus_df, on="pt", how="left")
    sampleid_to_virus_df = sampleid_to_virus_df.rename({"sampleid": "sample"})
    sampleid_to_virus_df = sampleid_to_virus_df.select(["sample", "pt"] + virus_List + [f"ct_{virus}" for virus in virus_List])
    return sampleid_to_virus_df

def compare_subtypes_between_virus_present_absent(snp_ratio_files: list, virus_file: str, virus: str, output_file: str, subtype_col: str="avg_subtype_ratio") -> None:
    virus_df = pl.read_csv(virus_file, separator='\t')
    # concat snp ratio dfs
    snp_dfs = [pl.read_csv(f, separator='\t') for f in snp_ratio_files]
    snp_df = pl.concat(snp_dfs)
    combined_df = snp_df.join(virus_df.select(["sample", virus, f"ct_{virus}"]), on="sample", how="left")
    fig = virus_swarmplot_and_ct_scatter(combined_df, virus_col=virus, snp_ratio_col=subtype_col, ct_col=f"ct_{virus}")
    fig.savefig(output_file, dpi=300)

def visualise_spatial_distribution_of_snp_ratios(snp_ratio_files: list, spatial_factors_file: str, lat_col: str, lon_col: str, snp_ratio_col: str, output_file: str) -> None:
    # concat snp ratio dfs
    snp_dfs = [pl.read_csv(f, separator='\t') for f in snp_ratio_files]
    snp_df = pl.concat(snp_dfs)
    # filter where all cols but sample are 0
    ratio_cols = [col for col in snp_df.columns if col != "sample"]
    snp_df = snp_df.filter(pl.sum_horizontal(ratio_cols) != 0)
    # combine with spatial factors
    spatial_df = pl.read_csv(spatial_factors_file, separator='\t')
    combined_df = snp_df.join(spatial_df, on="sample", how="left")
    # make plot
    scatter_lat_lon_coloured_by_snp_ratios(combined_df, lat_col=lat_col, lon_col=lon_col, snp_ratio_col=snp_ratio_col, output_file=output_file)

def visualise_spatial_distribution_of_snp_ratios_by_season(snp_ratio_files: list, spatial_factors_file: str, temporal_factors_file: str, lat_col: str, lon_col: str, snp_ratio_col: str, season_col: str, output_file: str) -> None:
    # concat snp ratio dfs
    snp_dfs = [pl.read_csv(f, separator='\t') for f in snp_ratio_files]
    snp_df = pl.concat(snp_dfs)
    # filter where all cols but sample are 0
    ratio_cols = [col for col in snp_df.columns if col != "sample"]
    snp_df = snp_df.filter(pl.sum_horizontal(ratio_cols) != 0)
    # combine with spatial factors
    spatial_df = pl.read_csv(spatial_factors_file, separator='\t')
    combined_df = snp_df.join(spatial_df, on="sample", how="left")
    # combine with temporal factors
    temporal_df = pl.read_csv(temporal_factors_file, separator='\t')
    combined_df = combined_df.join(temporal_df, on="sample", how="left")
    # make plot
    fig = scatter_lat_lon_coloured_by_snp_ratio_split_by_season(combined_df, lat_col=lat_col, lon_col=lon_col, snp_ratio_col=snp_ratio_col, season_col=season_col)
    fig.savefig(output_file, dpi=300)