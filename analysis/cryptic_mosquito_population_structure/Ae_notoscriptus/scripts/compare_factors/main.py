from analysis.Cx_annulirostris.compare_factors.process import (compare_subtype_snp_ratio_with_spatial_factors, 
                                                               generate_per_sample_spatial_factor_table, 
                                                               generate_per_sample_virus_file,
                                                               generate_per_sample_temporal_factor_table,
                                                               compare_subtypes_between_virus_present_absent, 
                                                               visualise_spatial_distribution_of_snp_ratios,
                                                               visualise_spatial_distribution_of_snp_ratios_by_season)
import os

def main():
    # inputs
    snp_ratio_files = [
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/COSMOS1.Ae_notoscriptus_collated_results.subtype_ratios_clusters1-3.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/COSMOS2.Ae_notoscriptus_collated_results.subtype_ratios_clusters1-3.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/COSMOS3.Ae_notoscriptus_collated_results.subtype_ratios_clusters1-3.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/collate_results/COSMOS4.Ae_notoscriptus_collated_results.subtype_ratios_clusters1-3.tsv"
    ]
    spatial_factors_list = [
        "lat",
        "lon",
        "distance_from_coast"
    ]
    location_to_latlon_file = "/Users/afor0966/Library/CloudStorage/Dropbox/COSMOS_Manuscript/Data/locs.csv"
    sampleid_to_pt_file = "/Users/afor0966/Library/CloudStorage/Dropbox/COSMOS_Manuscript/Data/sampleid2pt.csv"
    pt_to_location_file = "/Users/afor0966/Library/CloudStorage/Dropbox/COSMOS_Manuscript/Data/COSMOS1+COSMOS2+COSMOS3+COSMOS4_metadata_sequenced_traps.csv"
    pt_to_dist_to_coast_file = "/Users/afor0966/Library/CloudStorage/Dropbox/DataSnapshot2025/Adil/data progression/points_with_distance_to_freshwater.csv"
    # out files
    spatial_factors_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/spatial_factors_per_sample.tsv"
    virus_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/virus_presence_per_sample.tsv"
    temporal_factors_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/temporal_factors_per_sample.tsv"
    # geo factors
    generate_per_sample_spatial_factor_table(
        location_to_latlon_file=location_to_latlon_file,
        sampleid_to_pt_file=sampleid_to_pt_file,
        pt_to_location_file=pt_to_location_file,
        pt_to_dist_to_coast_file=pt_to_dist_to_coast_file
    ).write_csv(spatial_factors_file, separator='\t')
    # virus factors
    generate_per_sample_virus_file(
        sampleid_to_pt_file=sampleid_to_pt_file,
        pt_to_virus_file=pt_to_location_file
    ).write_csv(virus_file, separator='\t')
    # temporal factors
    generate_per_sample_temporal_factor_table(
        sampleid_to_pt_file=sampleid_to_pt_file,
        pt_to_season_file=pt_to_location_file
    ).write_csv(temporal_factors_file, separator='\t')
    # iterate over clusters
    for cluster in ["1", "2", "3"]:
        spatial_png = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/scatter.Ae_notoscriptus_cluster{cluster}_ratio_vs_geo.cluster={cluster}.png"
        spatial_dist_png = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/scatter.Ae_notoscriptus_cluster{cluster}_ratio_vs_geo_spatial_dist.cluster={cluster}.png"
        virus_png = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/compare_factors/swarmplot.Ae_notoscriptus_cluster{cluster}_ratio_by_RRV_presence.cluster={cluster}.png"
        # make sure output dir exists
        os.makedirs(os.path.dirname(spatial_png), exist_ok=True)
        os.makedirs(os.path.dirname(spatial_dist_png), exist_ok=True)
        # plots
        compare_subtype_snp_ratio_with_spatial_factors(
            snp_ratio_files, spatial_factors_file, spatial_factors_list, spatial_png, subtype_col=f"proportion_cluster{cluster}"
        )
        compare_subtypes_between_virus_present_absent(
            snp_ratio_files,
            virus_file,
            virus="RRV",
            output_file=virus_png,
            subtype_col=f"proportion_cluster{cluster}"
        )
        visualise_spatial_distribution_of_snp_ratios(
            snp_ratio_files, spatial_factors_file, lat_col="lat", lon_col="lon", 
            output_file=spatial_dist_png,
            snp_ratio_col=f"proportion_cluster{cluster}"
        )
        visualise_spatial_distribution_of_snp_ratios_by_season(
            snp_ratio_files, spatial_factors_file, temporal_factors_file, lat_col="lat", lon_col="lon", season_col="season",
            output_file=spatial_dist_png.replace("geo_spatial_dist", "geo_spatial_dist_by_season"),
            snp_ratio_col=f"proportion_cluster{cluster}"
        )

if __name__ == "__main__":
    main()