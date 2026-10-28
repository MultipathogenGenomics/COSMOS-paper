from .process import (
    prepare_abundance_geospatial_data,
    prepare_embedding_comparison_data,
    write_abundance_geospatial_description,
    write_pc_distance_to_coast_associations,
    write_simpsons_diversity_index,
    write_species_trap_counts,
)
from .plot import (
    plot_abundance_geospatial_overview,
    plot_embedding_comparison,
    plot_pc1_distance_to_coast_association,
    plot_pca_loadings,
    plot_pca_variance_explained,
)
import os


def main():
    assignment_reference = "COSMOS-plus-SS-pools"
    estimate_abundance_dir = f"results/estimate-abundance/{assignment_reference}"
    sample_sets = ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]
    abundance_filename = "abundance_within_filter_pass_reads.tsv"
    embedding_method = "PCA"
    pca_n_components = 5
    metadata_file = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    # virus_nc2_file = "data/castanetv9.2-output/COSMOS1-4_castanetv9.2_nc2.tsv"
    virus_nc2_file = "data/final_castanet_data/COSMOS1-4.final_castanet_nc2.tsv" # updated to final castanet data
    s_aus_file_template = "data/cryptic_subpops/{dataset}.Cx_annulirostris_collated_results.tsv"
    output_dir = "results/abundance-geospatial/plots"
    table_dir = "results/abundance-geospatial/tables"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(table_dir, exist_ok=True)

    geospatial_data = prepare_abundance_geospatial_data(
        estimate_abundance_dir,
        sample_sets,
        metadata_file,
        abundance_filename,
        embedding_method=embedding_method,
        pca_n_components=pca_n_components,
    )
    output_plot_file = f"{output_dir}/abundance_geospatial_overview.png"
    plot_abundance_geospatial_overview(geospatial_data, output_plot_file)
    embedding_comparison_df = prepare_embedding_comparison_data(
        estimate_abundance_dir,
        sample_sets,
        metadata_file,
        abundance_filename,
        virus_nc2_file,
        virus_acronyms = {
            "stratfordvirus": "STRV",
            "rossrivervirus": "RRV",
            "murrayvalleyencephalitisvirus": "MVEV",
            "kunjinvirus": "KUN",
            "kokoberavirus": "KOK",
            "edgehillvirus": "EHV",
            "japaneseencephalitisvirus": "JEV",
            "barmahforestvirus": "BFV",
            "sindbisvirus": "SINDV",
            "umbrevirus": "LSV",
        },
        s_aus_file_template=s_aus_file_template,
    )
    plot_embedding_comparison(
        embedding_comparison_df,
        [
            {
                "column": "region_category",
                "label": "Region category",
                "type": "categorical",
                "order": ["Met", "Urban", "Met-fringe", "Rural"],
                "palette": {
                    "Met": "#0072B2",
                    "Urban": "#E69F00",
                    "Met-fringe": "#CC79A7",
                    "Rural": "#009E73",
                },
            },
            {
                "column": "Coastal",
                "label": "Coastal/Inland",
                "type": "categorical",
                "order": [1, 0],
                "labels": {1: "coastal", 0: "inland"},
                "palette": {1: "#8EBFE8", 0: "#8FD19E"},
            },
            {
                "column": "virus_positive",
                "label": "Virus +ve/-ve",
                "type": "categorical",
                "order": [True, False],
                "labels": {True: "virus detected", False: "no virus detected"},
                "palette": {True: "#D55E00", False: "#D9D9D9"},
            },
            {
                "column": "s_aus_proportion",
                "label": "Cx. annulirostris S-AUS proportion",
                "type": "continuous",
                "cmap": "viridis",
                "vmin": 0,
                "vmax": 1,
            },
            {
                "column": "diversity_1_minus_D",
                "label": "Diversity index (1 − D)",
                "type": "continuous",
                "cmap": "magma",
                "vmin": 0,
                "vmax": 1,
            },
            {
                "column": "season",
                "label": "Season",
                "type": "categorical",
                "order": [2122, 2223, 2324],
                "labels": {2122: "2021/22", 2223: "2022/23", 2324: "2023/24"},
                "palette": {2122: "#8EBFE8", 2223: "#8FD19E", 2324: "#D55E00"}
            }
        ],
        f"{output_dir}/abundance_geospatial_embedding_comparison.png",
    )
    write_abundance_geospatial_description(
        geospatial_data,
        f"{output_dir}/abundance_geospatial_overview.description.txt",
        estimate_abundance_dir,
        sample_sets,
        metadata_file,
        abundance_filename,
        output_plot_file,
    )
    write_species_trap_counts(
        geospatial_data["species_trap_counts"],
        f"{table_dir}/species_trap_counts.tsv",
    )
    write_simpsons_diversity_index(
        geospatial_data["sdi"],
        f"{table_dir}/simpsons_diversity_index.tsv",
    )
    write_simpsons_diversity_index(
        geospatial_data["sdi_by_council_type"],
        f"{table_dir}/simpsons_diversity_index_by_council_type.tsv",
    )
    if geospatial_data["embedding_method"] == "PCA":
        plot_pca_loadings(
            geospatial_data["pca_loadings"],
            f"{output_dir}/abundance_geospatial_pca_loadings.png",
        )
        plot_pca_variance_explained(
            geospatial_data["pca_variance"],
            f"{output_dir}/abundance_geospatial_pca_variance_explained.png",
        )
        write_pc_distance_to_coast_associations(
            geospatial_data["pc_distance_to_coast_associations"],
            f"{table_dir}/abundance_geospatial_pc_distance_to_coast.tsv",
        )
        plot_pc1_distance_to_coast_association(
            geospatial_data["embedding"],
            f"{output_dir}/abundance_geospatial_pc1_distance_to_coast.png",
        )


if __name__ == "__main__":
    main()
