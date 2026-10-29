from pathlib import Path
from analysis.common.paths import RESULTS_DIR
import argparse
import numpy as np
from analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.process import process_fasta_and_cluster, generate_representative_targets_json, trim_fasta_to_no_gaps_in_msa, remove_targets_by_identity, build_consensus_sequence_for_clusters, trim_overhang
from analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.plots import build_tree_plot_from_fasta, build_pid_distribution_plots_from_fasta

def parse_args():
    parser = argparse.ArgumentParser(description="Cluster target sequences by percent identity (PID).")
    # required args
    parser.add_argument("--msa", type=str, help="Path to the input FASTA file containing target sequences.")
    parser.add_argument("--pid_threshold", type=float, default=0.95, help="Percent identity threshold for clustering (default: 0.95).")
    parser.add_argument("--output_prefix", type=str, default=str(RESULTS_DIR / "cluster_targets_by_pid" / "clusters"), help="Prefix for output files (default: results/cluster_targets_by_pid/clusters).")
    # optional flags
    parser.add_argument("--trim_fasta", required=False, action="store_true", help="If present will trim the input FASTA to corresponding lenths based on the  MSA.")
    parser.add_argument("--get_representatives", required=False, action="store_true", help="Whether to save representative targets to JSON file (optional).")
    parser.add_argument("--plot_tree", required=False, action="store_true", help="Whether to make a dendrogram plot (optional).")
    parser.add_argument("--plot_distributions", required=False,action="store_true",  help="Whether to make PID distribution plots (optional).")
    parser.add_argument("--get_consensuses", required=False, action="store_true", help="Whether to make consensus sequences for each cluster (default: False).")
    parser.add_argument("--allowed_gap_size", required=False, type=int, default=3, help="When trimming FASTA to MSA, allow up to this many gaps in a row in the MSA without trimming (default: 3).")
    # args for removing targets for testing
    parser.add_argument("--remove_target_for_testing", required=False, type=str, help="Target ID to remove for testing purposes (optional).")
    parser.add_argument("--remove_thresholds", required=False, type=float, nargs='+', help="List of PID thresholds for removing targets (optional).")
    parser.add_argument("--remove_output_json_prefix", required=False, type=str, help="Prefix for output JSON files when removing targets (optional).")
    return parser.parse_args()

def main():
    args = parse_args()
    Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)

    # trim overhang - if one or two sequences has long left or right hang, cut edges
    trim_overhang(args.msa, args.msa, n_samples_to_distinguish_overhang=15)

    # trim fasta if specified
    if args.trim_fasta:
        print("Trimming input FASTA...")
        trimmed_fasta_path = args.msa.replace(".fasta", "_trimmed.fasta")
        trim_fasta_to_no_gaps_in_msa(args.msa, trimmed_fasta_path, allow_gap=args.allowed_gap_size)
        args.msa = trimmed_fasta_path
    else:
        print("Not trimming input FASTA...")

    # make clusters
    cluster_df = process_fasta_and_cluster(args.msa, args.pid_threshold)
    to_save = cluster_df.sort(by=["cluster_id", "is_best_representative"], descending=[False, True])
    to_save.write_csv(args.output_prefix + ".clusters.tsv", separator="\t")
    print(f"Clusters at {args.pid_threshold * 100}% identity saved to {args.output_prefix + '.clusters.tsv'}.\n")

    # make a json of representatives if specified
    if args.get_representatives:
        output_json_path = args.output_prefix + ".representatives.json"
        if args.trim_fasta:
            json_key = f"clustered_{args.pid_threshold}_pid_trimmed"
        else:
            json_key = f"clustered_{args.pid_threshold}_pid"
        generate_representative_targets_json(cluster_df, output_json_path, json_key=json_key)
        print(f"Representative targets JSON saved to {output_json_path}.\n")

    # make tree plot if specified
    if args.plot_tree:
        tree_plot_path = args.output_prefix + ".tree.png"
        build_tree_plot_from_fasta(cluster_df, args.msa, args.pid_threshold, tree_plot_path)
        print(f"Pairwise distance tree plot saved to {tree_plot_path}.\n")

    # make pid distribution plots if specified
    if args.plot_distributions:
        distributions_plot_path = args.output_prefix + ".pid_distributions.png"
        build_pid_distribution_plots_from_fasta(args.msa, distributions_plot_path)
        print(f"PID distribution plots saved to {distributions_plot_path}.\n")

    # make consensus sequences if specified
    if args.get_consensuses:
        consensus_fasta_path = args.output_prefix + ".consensus_sequences.fasta"
        consensus_json_path = args.output_prefix + ".consensus_sequences.json"
        if args.trim_fasta:
            json_key = f"consensus_clusters_{args.pid_threshold}_pid_trimmed"
        else:
            json_key = f"consensus_clusters_{args.pid_threshold}_pid"
        build_consensus_sequence_for_clusters(cluster_df, args.pid_threshold, args.msa, consensus_fasta_path, consensus_json_path, json_key=json_key, trimmed=args.trim_fasta)
        print(f"Consensus sequences saved to {consensus_fasta_path}.\n")
    
    # remove target for testing if specified
    if args.remove_target_for_testing:
        print(f"Removing {args.remove_target_for_testing} alongside targets at {','.join(map(str, args.remove_thresholds))} percent identity thresholds for testing purposes.\n")
        if not args.remove_thresholds:
            args.remove_thresholds = [0.8, 0.90, 0.95, 0.98]
        if not args.remove_output_json_prefix:
            args.remove_output_json_prefix = str(RESULTS_DIR / "cluster_targets_by_pid" / f"representative_targets_with_{args.remove_target_for_testing}_removed")
        remove_targets_by_identity(cluster_df, args.pid_threshold, args.msa, args.remove_target_for_testing, np.array(args.remove_thresholds), args.remove_output_json_prefix)

if __name__ == "__main__":
    main()