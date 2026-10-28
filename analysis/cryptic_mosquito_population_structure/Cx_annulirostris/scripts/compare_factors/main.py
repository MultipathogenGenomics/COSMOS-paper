from analysis.Cx_annulirostris.compare_factors.process import compare_subtype_snp_ratio_with_spatial_factors, generate_per_sample_spatial_factor_table, generate_per_sample_virus_file, compare_subtypes_between_virus_present_absent, visualise_spatial_distribution_of_snp_ratios
import os

def main():
    # inputs
    snp_ratio_files = [
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/collate_results/COSMOS1.Cx_annulirostris_collated_results.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/collate_results/COSMOS2.Cx_annulirostris_collated_results.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/collate_results/COSMOS3.Cx_annulirostris_collated_results.tsv",
        "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/collate_results/COSMOS4.Cx_annulirostris_collated_results.tsv"
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
    spatial_factors_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/compare_factors/spatial_factors_per_sample.tsv"
    virus_file = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/compare_factors/virus_presence_per_sample.tsv"
    spatial_png = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/compare_factors/scatter.Cx_annulirostris_S-AUS_ratio_vs_geo.png"
    spatial_dist_png = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/compare_factors/scatter.Cx_annulirostris_S-AUS_ratio_vs_geo_spatial_dist.png"
    virus_png = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris/compare_factors/swarmplot.Cx_annulirostris_S-AUS_ratio_by_RRV_presence.png"
    # make sure output dir exists
    os.makedirs(os.path.dirname(spatial_png), exist_ok=True)
    os.makedirs(os.path.dirname(spatial_dist_png), exist_ok=True)
    # analysis
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
    # plots
    compare_subtype_snp_ratio_with_spatial_factors(
        snp_ratio_files, spatial_factors_file, spatial_factors_list, spatial_png
    )
    compare_subtypes_between_virus_present_absent(
        snp_ratio_files,
        virus_file,
        virus="RRV",
        output_file=virus_png
    )
    visualise_spatial_distribution_of_snp_ratios(
        snp_ratio_files, spatial_factors_file, lat_col="lat", lon_col="lon", snp_ratio_col="avg_subtype_ratio", output_file=spatial_dist_png
    )

if __name__ == "__main__":
    main()