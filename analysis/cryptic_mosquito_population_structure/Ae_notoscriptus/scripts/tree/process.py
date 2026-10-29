from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.determine_subtypes.process import msa_with_mafft
from analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.process import parse_fasta_file, compute_pairwise_differences_between_sequences, cluster_sequences_by_pid
from analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.plots import plot_tree
from analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.determine_subtypes.process import (find_Cx_annulirostris_frag_position_in_ref,
                                                                  get_ids_from_fasta, 
                                                                  extract_Cx_annulirostris_amplicon_region_from_reference,
                                                                  find_unique_kmers_for_subtype)
import os
import polars as pl

def build_Ae_notoscriptus_tree(fasta_path:str, output_path:str, tmp_dir: str) -> None:
    # make msa
    msa_output = os.path.join(tmp_dir, "Ae_notoscriptus_msa.fasta")
    msa_with_mafft(fasta_path, msa_output)
    # plot tree
    sequence_df = parse_fasta_file(msa_output)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(sequence_df)
    plot_tree(dist_matrix, ids, pid_threshold=0.95, output_path=output_path)

def minimise_fasta_if_identical_sequences(fasta_path: str, output_path: str) -> None:
    sequence_df = parse_fasta_file(fasta_path)
    unique_sequences = sequence_df.unique(subset="sequence", keep="first")
    current_sequence_n = len(unique_sequences)
    with open(output_path, 'w') as out_f:
        for row in unique_sequences.iter_rows(named=True):
            out_f.write(f">{row['id']}\n{row['sequence']}\n")
    print(f"Minimised fasta from {fasta_path} to {output_path}, reduced from {len(sequence_df)} to {current_sequence_n} sequences")

def build_Ae_notoscriptus_tree_including_frags(fasta_path: str, fragments_path: str, output_path: str, tmp_dir: str) -> None:
    # find positions of fragments in each contig in fasta path and put in new fasta file
    contigs = get_ids_from_fasta(fasta_path)
    fasta_dir = os.path.dirname(fasta_path)
    fragments_in_refs_path = os.path.join(fasta_dir, "Ae_notoscriptus_fragments_in_references.fasta")
    with open(fragments_in_refs_path, 'w') as out_f:
        for contig in contigs:
            positions = find_Cx_annulirostris_frag_position_in_ref(fragments_path, fasta_path, contig, 31, tmp_dir)
            start = positions["start"]["start"]
            end = positions["end"]["end"]
            fragment_seq = extract_Cx_annulirostris_amplicon_region_from_reference(positions, fasta_path, contig)
            out_f.write(f">{contig}\n{fragment_seq}\n")
    # minimise fragments fasta if identical sequences
    minimised_fragments_path = os.path.join(tmp_dir, "Ae_notoscriptus_fragments_minimised.fasta")
    minimise_fasta_if_identical_sequences(fragments_path, minimised_fragments_path)
    # combine with fragments fasta
    combined_fasta_path = os.path.join(tmp_dir, "Ae_notoscriptus_fragments_combined.fasta")
    combine_two_fastas(fragments_in_refs_path, minimised_fragments_path, combined_fasta_path)
    # make msa
    msa_output = os.path.join(tmp_dir, "Ae_notoscriptus.COI_fragments_msa.fasta")
    msa_with_mafft(combined_fasta_path, msa_output)
    # plot tree
    sequence_df = parse_fasta_file(msa_output)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(sequence_df)
    plot_tree(dist_matrix, ids, pid_threshold=0.97, output_path=output_path, bold_leaves=contigs)
    # cluster all sequences (not minimised) by pid at 0.97
    combined_fasta_path = os.path.join(tmp_dir, "Ae_notoscriptus_fragments_combined.all.fasta")
    combine_two_fastas(fragments_path, fragments_in_refs_path, combined_fasta_path)
    msa_output = os.path.join(tmp_dir, "Ae_notoscriptus_fragments_combined.all.msa.fasta")
    msa_with_mafft(combined_fasta_path, msa_output)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(parse_fasta_file(combined_fasta_path))
    cluster_df = cluster_sequences_by_pid(ids, dist_matrix, pid_threshold=0.97)
    results_dir = os.path.dirname(output_path)
    cluster_output = os.path.join(results_dir, "Ae_notoscriptus_fragments_and_refs_combined.COI_cluster_at_0.97.tsv")
    cluster_df.write_csv(cluster_output, separator='\t')

def combine_two_fastas(fasta1: str, fasta2: str, output_path: str) -> None:
    with open(output_path, 'w') as out_f:
        with open(fasta1, 'r') as f1:
            for line in f1:
                out_f.write(line)
        with open(fasta2, 'r') as f2:
            for line in f2:
                out_f.write(line)

def find_kmers_for_each_Ae_notoscriptus_cluster(cluster_info: str, fragment_fasta: str, reference_fasta: str, kmer_length: int, output_path: str) -> None:
    # parse cluster df
    cluster_df = pl.read_csv(cluster_info, separator='\t')
    # combined fragments and references fasta
    combined_fasta_path = os.path.join(os.path.dirname(output_path), "Ae_notoscriptus_fragments_combined.all.fasta")
    combine_two_fastas(fragment_fasta, reference_fasta, combined_fasta_path)
    # rename cluster col to subtype 
    cluster_df = cluster_df.rename({"cluster_id": "subtype", "seq_id": "fasta_id"})
    # find unique kmers for each cluster
    unique_kmers_df = find_unique_kmers_for_subtype(cluster_df, combined_fasta_path, kmer_length)
    unique_kmers_df.write_csv(output_path, separator='\t')

