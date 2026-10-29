from analysis.common.paths import RESULTS_DIR
from analysis.enrichment.scripts.assess_enrichment.align import align_reads_to_mt_genomes, run_msa_of_mt_genomes, write_alignment_reference_map
from analysis.enrichment.scripts.assess_enrichment.coverage import compute_coverage_of_mt_genomes_in_alignments
from analysis.enrichment.scripts.assess_enrichment.process import (
    add_alignment_positions_to_coverage_data,
    compute_enrichment_summary_table,
    compute_median_mean_in_cox1,
    concat_coverage_data,
    write_overall_enrichment_summary_table,
)
from analysis.enrichment.scripts.assess_enrichment.plot import plot_enrichment_over_COX1, plot_enrichment_over_COX1_with_BFV_culture
import os

def main():
    # inputs
    genomes_config = "config/mt_genomes.json"
    samples_config = "config/single_species_controls.trimmed.json"
    genus_defaults = {
        "Aedes": "data/mt_genomes/Aedes_vigilax_mt_genome.fasta",
        "Culex": "data/mt_genomes/Culex_annulirostris_mt_genome.fasta",
        "Anopheles": "data/mt_genomes/Anopheles_annulipes_mt_genome.fasta",
        "Mansonia": "data/mt_genomes/Mansonia_uniformis_mt_genome.fasta",
        "Coquillettidia": "data/mt_genomes/Coquillettidia_linealis_mt_genome.fasta",
        "Ochlerotatus": "data/mt_genomes/Aedes_vigilax_mt_genome.fasta"
    }
    BFV_results = "data/BFVculture-dilutions/bfv_expt_COSMOS3.csv" # for scatter plot
    # write the BAM-to-reference map without remaking BAMs
    os.makedirs(f"{RESULTS_DIR}/assess-enrichment/alignments", exist_ok=True)
    write_alignment_reference_map(samples_config, genomes_config, genus_defaults, f"{RESULTS_DIR}/assess-enrichment/alignments")
    # # align reads to mt genomes
    # align_reads_to_mt_genomes(samples_config, genomes_config, genus_defaults, "results/assess-enrichment/alignments")
    # run_msa_of_mt_genomes(genomes_config, "results/assess-enrichment/alignments")
    # # compute coverage of mt genomes in alignments
    # os.makedirs("results/assess-enrichment/coverage", exist_ok=True)
    # compute_coverage_of_mt_genomes_in_alignments("results/assess-enrichment/alignments", "results/assess-enrichment/coverage")
    # collate data
    coverage_df = concat_coverage_data(samples_config, genomes_config, f"{RESULTS_DIR}/assess-enrichment/coverage", genus_defaults)
    coverage_df = add_alignment_positions_to_coverage_data(coverage_df, f"{RESULTS_DIR}/assess-enrichment/alignments/references_aligned.fasta")
    fold_enrichment_df = compute_median_mean_in_cox1(coverage_df)
    enrichment_summary_df = compute_enrichment_summary_table(
        coverage_df,
        f"{RESULTS_DIR}/assess-enrichment/alignments/alignment_reference_map.tsv",
        genomes_config,
        f"{RESULTS_DIR}/assess-enrichment/tables/enrichment_summary.tsv",
    )
    write_overall_enrichment_summary_table(
        enrichment_summary_df,
        f"{RESULTS_DIR}/assess-enrichment/tables/enrichment_summary_overall.tsv",
    )
    # plot coverage across COX1
    os.makedirs(f"{RESULTS_DIR}/assess-enrichment/plots", exist_ok=True)
    plot_enrichment_over_COX1(coverage_df, fold_enrichment_df, f"{RESULTS_DIR}/assess-enrichment/plots")
    plot_enrichment_over_COX1_with_BFV_culture(coverage_df, BFV_results, f"{RESULTS_DIR}/assess-enrichment/plots")

if __name__ == "__main__":
    main()
