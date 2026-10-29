from analysis.common.paths import RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.process import prepare_cx_annulirostris_coastal_association
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.plot import plot_cx_annulirostris_coastal_association
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.process import prepare_cryptic_subpopulation_data, write_subpopulation_detection_summary, combine_Cx_annulirostris_subpopulation_kmer_count_data
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.plot import plot_ae_notoscriptus_population_latitude_associations, plot_cryptic_subpopulations_overview, plot_sequence_dendrogram, heatmap_Cx_annulirostris_subtype_kmer_detections, plot_iqtree_tree
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.trees import prepare_sequence_tree, make_phylogenetic_trees_for_species
import os
import pandas as pd


def main():
    geopkg_file = "data/aus_maps/ASGS_Ed3_Non_ABS_Structures_GDA2020_updated_2025.gpkg"
    metadata_file = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    output_dir = f"{RESULTS_DIR}/cryptic-subpopulations/plots"
    table_dir = f"{RESULTS_DIR}/cryptic-subpopulations/tables"
    datasets = ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]
    estimate_abundance_dir = f"{RESULTS_DIR}/estimate-abundance/COSMOS-plus-SS-pools"
    abundance_filename = "abundance_within_filter_pass_reads.tsv"
    tree_output_dir = f"{RESULTS_DIR}/cryptic-subpopulations/trees"

    subtype_colours = {
        "S-AUS": "#FFDF50",
        "AUS": "#ED4B00",
        "PNG-1": "#3F3F3F",
        "PNG-2": "#717171A8",
        "SI": "#BCBBBBBE",
        "population 1": "#0072B2",  # blue
        "population 2": "#009E73",  # green/teal
        "population 3": "#EC8CEE",
    }

    panel_configs = [
        {
            "label": "Culex annulirostris",
            "species_abundance_column": "Culex annulirostris",
            "file_template": "data/cryptic_subpops/{dataset}.Cx_annulirostris_collated_results.tsv",
            "ratio_columns": ["avg_subtype_ratio", "non_subtype_ratio"],
            "legend_labels": ["S-AUS", "AUS"],
            "colours": [subtype_colours["S-AUS"], subtype_colours["AUS"]],
            "gradient_colours": [subtype_colours["AUS"], subtype_colours["S-AUS"]],
            "min_count_threshold": 20,
            "complement_ratio_column": "avg_subtype_ratio",
            "complement_output_column": "non_subtype_ratio",
            "summary_filters": "metadata samples only; min_count > 20; valid avg_subtype_ratio; latitude and longitude present",
            "to_add_to_legend": {"PNG-1": subtype_colours["PNG-1"], "PNG-2": subtype_colours["PNG-2"], "SI": subtype_colours["SI"]},
            "draw_type": "gradient"
        },
        {
            "label": "Aedes notoscriptus",
            "species_abundance_column": "Aedes notoscriptus",
            "source": "ae_notoscriptus_pairwise",
            "one_vs_23_file_template": "data/cryptic_subpops/{dataset}.Ae_notoscriptus_collated_results.1_vs_2or3.tsv",
            "two_vs_3_file_template": "data/cryptic_subpops/{dataset}.Ae_notoscriptus_collated_results.2_vs_3.tsv",
            "one_vs_23_min_snps": 12,
            "two_vs_3_min_snps": 5,
            "ratio_columns": ["proportion_cluster1", "proportion_cluster2", "proportion_cluster3"],
            "legend_labels": ["population 1", "population 2", "population 3"],
            "colours": [subtype_colours["population 1"], subtype_colours["population 2"], subtype_colours["population 3"]],
            "summary_filters": "metadata samples only; valid 1-vs-2/3 and 2-vs-3 ratios; n_snps_which_meet_threshold >= 12 for 1-vs-2/3; n_snps_which_meet_threshold >= 5 for 2-vs-3; latitude and longitude present",
            "to_add_to_legend": None,
            "draw_type": "pie"
        },
    ]

    tree_configs = [
        {
            "species_name": "Culex annulirostris",
            "output_prefix": "cx_annulirostris",
            "sequence_fastas": [
                "data/references/lineages/Cx_annulirostris_COX1_fragments_Hemmerter_et_al.fasta",
            ],
            # "pid_threshold": 0.97,
            "min_shared_sites": 100,
            "leaf_annotation": "cx_type",
            "metadata_file": "data/references/lineages/Cx_annulirostris_COX1_lineages_Hemmerter_et_al.tsv",
        },
        {
            "species_name": "Aedes notoscriptus",
            "output_prefix": "ae_notoscriptus",
            "sequence_fastas": [
                "data/references/lineages/Ae_notoscriptus_COX1_fragments_Endersby_et_al.fasta",
                "data/references/lineages/Ae_notoscriptus_COX1_fragments_from_reference_panel.fasta",
            ],
            "pid_threshold": 0.97,
            "min_shared_sites": 100,
            "leaf_annotation": "ae_population",
        },
    ]

    # os.makedirs(output_dir, exist_ok=True)
    # # os.makedirs(table_dir, exist_ok=True)
    # subpopulation_data = prepare_cryptic_subpopulation_data(
    #     geopkg_file,
    #     metadata_file,
    #     panel_configs,
    #     datasets,
    # )
    # plot_cryptic_subpopulations_overview(
    #     subpopulation_data,
    #     f"{output_dir}/cryptic_subpopulations_overview.png",
    # )
    # # write_subpopulation_detection_summary(
    # #     subpopulation_data,
    # #     f"{table_dir}/subpopulation_detection_summary.tsv",
    # #     estimate_abundance_dir,
    # #     datasets,
    # #     abundance_filename,
    # # )
    # ae_panel = next(panel for panel in subpopulation_data["panels"] if panel["label"] == "Aedes notoscriptus")
    # plot_ae_notoscriptus_population_latitude_associations(
    #     ae_panel,
    #     f"{output_dir}/ae_notoscriptus_populations_vs_latitude.png",
    #     f"{output_dir}/ae_notoscriptus_populations_vs_latitude.stats.tsv",
    # )
    cx_config = panel_configs[0]
    coastal_data = prepare_cx_annulirostris_coastal_association(metadata_file, cx_config["file_template"], datasets)
    plot_cx_annulirostris_coastal_association(
        coastal_data,
        f"{output_dir}/cx_annulirostris_s_aus_vs_distance_from_coast.png",
        f"{output_dir}/cx_annulirostris_s_aus_vs_distance_from_coast.stats.tsv",
        metadata_file=metadata_file,
        input_files=[cx_config["file_template"].format(dataset=dataset) for dataset in datasets],
    )
    # supplementary_tables = []
    # for tree_config in tree_configs:
    #     tree_data = prepare_sequence_tree(
    #         tree_config,
    #         tree_output_dir,
    #     )
    #     supplementary_tables.append(tree_data["supplementary_table"])
    #     plot_sequence_dendrogram(
    #         tree_data,
    #         f"{output_dir}/{tree_config['output_prefix']}_sequence_dendrogram.png",
    #     )
    # pd.concat(supplementary_tables, ignore_index=True).to_csv(
    #     f"{table_dir}/sequence_tree_metadata.tsv",
    #     sep="\t",
    #     index=False,
    # )

    # plot a heatmap of Cx annu subpop kmer counts
    os.makedirs(output_dir, exist_ok=True)
    subpopulation_data = prepare_cryptic_subpopulation_data(
        geopkg_file,
        metadata_file,
        panel_configs,
        datasets,
    )
    kmer_counts_dir = "data/cryptic_subpops/subtype_kmer_counts"
    cx_panel = next(panel for panel in subpopulation_data["panels"] if panel["label"] == "Culex annulirostris")
    kmer_count_df = combine_Cx_annulirostris_subpopulation_kmer_count_data(
        kmer_counts_dir,
        samples_to_keep=cx_panel["points"]["sample"].to_list(),
        subtype_fraction_threshold=0.95,
    )
    heatmap_Cx_annulirostris_subtype_kmer_detections(
        kmer_count_df,
        f"{output_dir}/Cx_annulirostris_subtype_kmer_detections.png",
        s_aus_proportions=cx_panel["points"].set_index("sample")["avg_subtype_ratio"],
    )

    # let's draw a phylogenetic tree of each species
    # make_phylogenetic_trees_for_species(tree_configs, tree_output_dir)
    label_metadata = f"{RESULTS_DIR}/cryptic-subpopulations/tables/sequence_tree_metadata.tsv"
    for tree_config in tree_configs:
        tree_file = os.path.join(tree_output_dir, f"{tree_config['output_prefix']}", f"{tree_config['output_prefix']}_phylogenetic_tree.treefile")
        plot_iqtree_tree(
            tree_file,
            f"{output_dir}/{tree_config['output_prefix']}_iqtree.png",
            label_metadata,
            subtype_colours=subtype_colours,
        )

if __name__ == "__main__":
    main()
