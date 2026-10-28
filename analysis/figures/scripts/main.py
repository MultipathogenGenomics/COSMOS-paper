from .process import make_mt_genome_coverage_file, parse_annulirostris_mt_coverage, prepare_heatmap_data
from .enrichment_and_estimation_figure import plot_enrichment_and_estimation, plot_enrichment_and_nc2, quick_plot_of_nc2_of_cox1_across_traps
from .results_figure import plot_results_figure
from .virus_scatters_and_heatmap import plot_virus_scatters_and_heatmap, supplementary_nc2_curve_figure, wuhan6_supp_figure, bfv_supp_figure
from ..common.abundance import read_abundance_for_sample_sets
import os
import polars as pl

def main():
    plot_dir = "results/publication_figures/plots"
    data_dir = "results/publication_figures/data"
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    metadata_file = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    virus_nc2_file = "data/final_castanet_data/COSMOS1-4.final_castanet_nc2.tsv"

    ## Figure 1 - infographic

    # ## Figure 2 (previous figure 1)

    mitochondria_coverage_file = os.path.join(data_dir, "mitochondria_coverage_plotting_data.csv")
    if not os.path.exists(mitochondria_coverage_file):
        make_mt_genome_coverage_file(mitochondria_coverage_file)

    # use TG Cx annu alignments instead
    mitochondria_annu_coverage_file = os.path.join(data_dir, "mitochondria_annu_coverage_plotting_data.csv")
    if not os.path.exists(mitochondria_annu_coverage_file):
        parse_annulirostris_mt_coverage(mitochondria_annu_coverage_file, metadata_file=metadata_file)

    heatmap_data_file = os.path.join(data_dir, "heatmap_data.pkl")
    if not os.path.exists(heatmap_data_file):
        prepare_heatmap_data(heatmap_data_file)

    # # virus_coverage_file = "results/virus/plots/virus_coverage_plotting_data.csv"
    # virus_coverage_file = "/Users/afor0966/Documents/code/COSMOS-analysis/data/final_castanet_data/COSMOS1-4.final_castanet_virus_coverage_per_sample.tsv"
    # output_file = f"{plot_dir}/enrichment_and_estimation.png"
    # plot_enrichment_and_estimation(mitochondria_annu_coverage_file, virus_coverage_file, heatmap_data_file, output_file)
    # virus_per_sample_coverage_file = "data/final_castanet_data/COSMOS1-4.final_castanet_virus_coverage_per_sample.tsv"
    # plot_enrichment_and_nc2(
    #     mitochondria_coverage_file=mitochondria_annu_coverage_file,
    #     virus_coverage_file=virus_per_sample_coverage_file,
    #     metadata_file=metadata_file,
    #     virus_nc2_file=virus_nc2_file,
    #     output_file=f"{plot_dir}/enrichment_and_nc2.final.png",
    #     fit_bin_medians=True,
    # )

    # # Figure 3 - new
    surveillance_validation_file = "results/surveillance-validation/surveillance_validation_data.tsv"
    heatmap_data_file = os.path.join(data_dir, "heatmap_data.pkl")
    virus_nc2_file = "data/final_castanet_data/COSMOS1-4.final_castanet_nc2.tsv"
    # supplementary_nc2_curve_figure(
    #     metadata_file,
    #     virus_nc2_file,
    #     f"{plot_dir}/supplementary_nc2_curve.final.png",
    #     fit_bin_medians=True,
    # )
    virus_dedup_file = "data/final_castanet_data/COSMOS1-4.final_castanet_reads_dedup.tsv"
    # virus_nc2_file = "data/castanetv9.3_output/COSMOS1-4_castanetv9.3_nc2.post_filt_false.tsv"
    # virus_dedup_file = "data/castanetv9.3_output/COSMOS1-4_castanetv9.3_reads_dedup.post_filt_false.tsv"
    output_file = f"{plot_dir}/validation_scatters_and_heatmap.final.png"
    plot_virus_scatters_and_heatmap(metadata_file, virus_nc2_file, virus_dedup_file, surveillance_validation_file, heatmap_data_file, output_file)

    # Figure 4 (previous figure 2)

    # SDI_file = "results/abundance-geospatial/tables/simpsons_diversity_index.tsv"
    # abundance_dir = "results/estimate-abundance/COSMOS-plus-SS-pools"
    # abundance_sample_sets = ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]
    # # virus_nc2_file = "data/castanetv9.2-output/COSMOS1-4_castanetv9.2_nc2.tsv"
    # NSW_map_file = "data/aus_maps/ASGS_Ed3_Non_ABS_Structures_GDA2020_updated_2025.gpkg"

    # combined_abundance_file = os.path.join(data_dir, "combined_abundance_data.tsv")
    # if not os.path.exists(combined_abundance_file):
    #     combined_abundance_df = read_abundance_for_sample_sets(abundance_dir, abundance_sample_sets, "abundance_within_filter_pass_reads.tsv")
    #     combined_abundance_df.write_csv(combined_abundance_file)

    # # supp figure - quick plot of nc2 of cox1 across traps
    # quick_plot_of_nc2_of_cox1_across_traps(
    #     nc2_file=virus_nc2_file,
    #     metadata_file=metadata_file,
    #     mosquito_abundance_results_file=combined_abundance_file,
    #     output_file=f"{plot_dir}/supplementary_cox1_nc2_across_traps.final.png",
    # )

    bfv_supp_figure()
    wuhan6_supp_figure()

    # plot_results_figure(
    #     metadata_file=metadata_file,
    #     SDI_results_file=SDI_file,
    #     mosquito_abundance_file=combined_abundance_file,
    #     virus_nc2_file=virus_nc2_file,
    #     NSW_map_file=NSW_map_file,
    #     out_file=f"{plot_dir}/results_figure.final.png"
    # )

if __name__ == "__main__":
    main()
