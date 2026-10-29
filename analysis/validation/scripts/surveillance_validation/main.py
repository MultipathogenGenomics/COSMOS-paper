from analysis.common.paths import RESULTS_DIR
from analysis.validation.scripts.surveillance_validation.process import (
    prepare_surveillance_validation_data,
    write_surveillance_validation_description,
    write_surveillance_validation_tables,
)
from analysis.validation.scripts.surveillance_validation.plot import (
    plot_pairwise_species_abundance_scatter,
    plot_surveillance_validation,
    write_pairwise_species_abundance_description,
)
import os


def main():
    assignment_reference = "COSMOS-plus-SS-pools"
    estimate_abundance_dir = f"{RESULTS_DIR}/estimate-abundance/{assignment_reference}"
    sample_sets = ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]
    abundance_filename = "abundance_within_filter_pass_reads.tsv"
    cosmos_metadata_file = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    nsw_counts_file = "data/COSMOS_metadata/NSW_health.mosq_counts.csv"
    output_dir = f"{RESULTS_DIR}/surveillance-validation"
    disagreement_threshold = 0.25
    pairwise_scatter_top_n = 10
    target_species = {
        "Aedes vigilax": "aedes_vigilax",
        "Culex annulirostris": "culex_annulirostris",
    }

    os.makedirs(output_dir, exist_ok=True)
    validation_data = prepare_surveillance_validation_data(
        estimate_abundance_dir,
        sample_sets,
        abundance_filename,
        cosmos_metadata_file,
        nsw_counts_file,
        target_species,
        assignment_reference=assignment_reference,
        disagreement_threshold=disagreement_threshold,
    )
    write_surveillance_validation_tables(validation_data, output_dir)

    output_plot_file = f"{output_dir}/surveillance_validation_scatter.png"
    plot_surveillance_validation(validation_data["validation_df"], output_plot_file)
    write_surveillance_validation_description(
        validation_data,
        f"{output_dir}/surveillance_validation_scatter.description.txt",
        output_plot_file,
    )

    pairwise_scatter_output = f"{output_dir}/pairwise_species_abundance_scatter.png"
    plot_pairwise_species_abundance_scatter(
        validation_data["abundance_df"],
        pairwise_scatter_output,
        top_n=pairwise_scatter_top_n,
    )
    write_pairwise_species_abundance_description(
        f"{output_dir}/pairwise_species_abundance_scatter.description.txt",
        pairwise_scatter_output,
        assignment_reference,
        abundance_filename,
        pairwise_scatter_top_n,
        validation_data["abundance_df"],
    )


if __name__ == "__main__":
    main()
