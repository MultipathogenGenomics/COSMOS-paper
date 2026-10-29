from analysis.validation.scripts.in_silico_mixtures.get_mosquito_reads import get_mosquito_reads_from_ss_pools, make_json_of_mosquito_fastq_files
from analysis.validation.scripts.in_silico_mixtures.mix import mix_fastq_files, make_snakemake_input_json
from analysis.validation.scripts.in_silico_mixtures.estimate_abundance import run_kallisto_and_kma_quantification
from analysis.validation.scripts.in_silico_mixtures.abundance_error import compute_abundance_error
from analysis.validation.scripts.in_silico_mixtures.plot import plot_in_silico_mixtures_results
from analysis.common.utils import filter_fasta_by_ids
from analysis.common.io import parse_json
import os

def main():
    # inputs
    sample_config = "config/single_species_controls.trimmed.json" # swap back to not trimmed if reverting
    extraction_reference_file = "data/references/COSMOS_original/castanet_mosquito_cox1_seqs.fasta"
    extraction_contig_string = "mosquitocox1_cluster0:472-1402 mosquitocox1_cluster1:472-1402 mosquitocox1_cluster2:472-1402" # matches COX1 trimming in the assignment reference file
    assignment_reference_file = "data/references/COSMOS_original/mozz_2024-11-01_min1kb_AF_trial_update_plus_additional_COX1s.no_LEGS.msa_trimmed_no_gaps.fasta"
    representative_ids = "data/references/COSMOS_original/mozz_2024-11-01_min1kb_AF_trial_update_plus_additional_COX1s.no_LEGS.pid_0.995.representatives.json"
    n_in_mixture = 4
    proportions = [0.25, 0.25, 0.25, 0.25]
    # get mosquito reads from the single species pools
    # os.makedirs("results/in-silico-mixtures/mosquito_fastq", exist_ok=True)
    # get_mosquito_reads_from_ss_pools(sample_config, extraction_reference_file, "results/in-silico-mixtures/mosquito_fastq", extraction_contig_string)
    updated_sample_json = "config/single_species_controls.mosquito_cox1_extracted_fastq.json"
    # make_json_of_mosquito_fastq_files(sample_config, "results/in-silico-mixtures/mosquito_fastq", updated_sample_json)
    # # filter the assignment reference file to only include the species in the sample config
    # os.makedirs("results/in-silico-mixtures/reference", exist_ok=True)
    # ids_to_keep = parse_json(representative_ids, remove_nest=True)['clustered_0.995_pid_trimmed']
    # filter_fasta_by_ids(assignment_reference_file, set(ids_to_keep), f"results/in-silico-mixtures/reference/representatives.fasta")
    # # make the in silico mixtures
    # os.makedirs("results/in-silico-mixtures/mixture_fastq", exist_ok=True)
    # mix_fastq_files(updated_sample_json, n_in_mixture, proportions, "results/in-silico-mixtures/mixture_fastq")
    # # run kallisto and kma quantification
    # os.makedirs("results/in-silico-mixtures/quantification", exist_ok=True)
    # run_kallisto_and_kma_quantification(f"results/in-silico-mixtures/reference/representatives.fasta", "results/in-silico-mixtures/mixture_fastq", "results/in-silico-mixtures/quantification")
    # # compute abundance estimation error
    # os.makedirs("results/in-silico-mixtures/error", exist_ok=True)
    # compute_abundance_error("results/in-silico-mixtures/mixture_fastq/mixture_metadata.tsv", "results/in-silico-mixtures/quantification/kallisto_kma_quantification_results.tsv", "results/in-silico-mixtures/error/abundance_error.tsv")
    # # plot results
    # os.makedirs("results/in-silico-mixtures/plots", exist_ok=True)
    # plot_in_silico_mixtures_results("results/in-silico-mixtures/error/abundance_error.tsv", "results/in-silico-mixtures/quantification/kallisto_kma_quantification_results.tsv", "results/in-silico-mixtures/mixture_fastq/mixture_metadata.tsv", "results/in-silico-mixtures/plots")
    # # make a json with in silico mixtures to run in snakemake
    # make_snakemake_input_json("results/in-silico-mixtures/mixture_fastq/mixture_metadata.tsv", "results/in-silico-mixtures/mixture_fastq", "config/in_silico_mixtures.json")
if __name__ == "__main__":
    main()
