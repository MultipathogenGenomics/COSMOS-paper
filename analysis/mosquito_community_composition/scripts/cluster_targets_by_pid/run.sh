#!/bin/bash

run_name=$1
fasta_file=$2
pid_thresholds=$3  # comma-separated list of PID thresholds, e.g. "0.9,0.95"
allowed_gap_size=$4  # integer, e.g. 3

mkdir -p results/cluster_targets_by_pid/"${run_name}"

# echo the parameters into a log file for record keeping
echo -e "Run name: ${run_name}\nFasta file: ${fasta_file}\nPID thresholds: ${pid_thresholds}\nAllowed gap size: ${allowed_gap_size}\n" > results/cluster_targets_by_pid/"${run_name}"/run_parameters.log

echo -e "\nClustering targets by PID for run: ${run_name} using thresholds ${pid_thresholds}\n"

for pid_threshold in ${pid_thresholds//,/ }; do

    # make a directory for results
    mkdir -p results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"
    mkdir -p results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/trimmed
    mkdir -p results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/not_trimmed

    # get basename for output file naming
    fasta_basename=$(basename "${fasta_file}" .fasta)

    # make msa of fasta using mafft
    if [ -f results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".msa.fasta ]; then
        echo -e "\nMSA already exists, skipping mafft step...\n"
    else
        echo -e "\nRunning mafft to create MSA...\n"
        mafft --auto "${fasta_file}" > results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".msa.fasta
    fi

    # run script
    python -m analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.main \
        --msa results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".msa.fasta \
        --trim_fasta \
        --pid_threshold "${pid_threshold}" \
        --output_prefix results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/trimmed/"${fasta_basename}".pid_"${pid_threshold}" \
        --get_representatives \
        --get_consensuses \
        --plot_tree \
        --plot_distributions \
        --allowed_gap_size $allowed_gap_size
        # --remove_target_for_testing "mosquito_NC_025473.1-1432-2968-Aedes-notoscriptus-isolate-MC04-mitochondrion" \
        # --remove_thresholds 0.80 0.88 0.90 0.92 0.95 0.98 \
        # --remove_output_json_prefix results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/trimmed/Ae_notoscriptus_MC04_representative_targets_by_"${pid_threshold}"_pid

    python -m analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.main \
        --msa results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".msa.fasta \
        --pid_threshold "${pid_threshold}" \
        --output_prefix results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/not_trimmed/"${fasta_basename}".pid_"${pid_threshold}" \
        --get_representatives \
        --get_consensuses \
        --plot_tree \
        --plot_distributions \
        --allowed_gap_size $allowed_gap_size     
        # --remove_target_for_testing "mosquito_NC_025473.1-1432-2968-Aedes-notoscriptus-isolate-MC04-mitochondrion" \
        # --remove_thresholds 0.80 0.88 0.90 0.92 0.95 0.98 \
        # --remove_output_json_prefix results/cluster_targets_by_pid/"${run_name}"/pid_threshold_"${pid_threshold}"/not_trimmed/Ae_notoscriptus_MC04_representative_targets_by_"${pid_threshold}"_pid
done

echo -e "Combining consensus FASTA files across thresholds for non_trimmed using seqkit...\n"

seqkit rmdup -n < <(cat results/cluster_targets_by_pid/"${run_name}"/pid_threshold_*/not_trimmed/*.consensus_sequences.fasta) \
    > results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".all_pid_thresholds.non_trimmed.consensus_sequences.fasta

echo -e "and now for trimmed...\n"

seqkit rmdup -n < <(cat results/cluster_targets_by_pid/"${run_name}"/pid_threshold_*/trimmed/*.consensus_sequences.fasta) \
    > results/cluster_targets_by_pid/"${run_name}"/"${fasta_basename}".all_pid_thresholds.trimmed.consensus_sequences.fasta

echo "finito!"