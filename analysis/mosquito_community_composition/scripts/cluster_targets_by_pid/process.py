from Bio import SeqIO
import polars as pl
import numpy as np
import scipy.cluster.hierarchy as sch
import json
import os
from analysis.common.mosquito_ids import return_species_from_target_id

def parse_fasta_file(fasta_path: str) -> pl.DataFrame:
    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_id = record.id.replace(" ", "_").replace(";", "").replace("(", "").replace(")", "").replace(",", "") # remove spaces and special characters from seq_id to avoid issues later
        sequence = str(record.seq)
        records.append({"id": seq_id, "sequence": sequence})
    
    df = pl.DataFrame(records)
    return df

def trim_fasta_to_no_gaps_in_msa(fasta_path: str, trimmed_fasta_path: str, allow_gap: int) -> None:
    # turn seqs in fasta into matrix of letters
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_id = record.id.replace(" ", "_").replace(";", "").replace("(", "").replace(")", "").replace(",", "") # remove spaces and special characters from seq_id to avoid issues later
        ids.append(seq_id)
        sequences.append(list(str(record.seq)))
    seq_matrix = np.array(sequences)
    # identify longest stretch of columns with no '-'
    gap_free_cols = np.all(seq_matrix != '-', axis=0)
    while np.sum(gap_free_cols) < 500: # this just means it likely won't trim it to very small
        # raise ValueError("No gap-free columns found in the MSA.")
        # try removing sequences with most gaps and check again
        gap_counts = np.sum(seq_matrix == '-', axis=1)
        most_gappy_index = np.argmax(gap_counts)
        print(f"No gap-free columns found. Removing sequence {ids[most_gappy_index]}.")
        seq_matrix = np.delete(seq_matrix, most_gappy_index, axis=0)
        ids.pop(most_gappy_index)
        gap_free_cols = np.all(seq_matrix != '-', axis=0)
    # find indices of stretches with no gaps
    gap_free_indices = np.where(gap_free_cols)[0]
    # find longest continuous stretch
    max_len = 0
    best_start = 0
    best_end = 0
    current_start = gap_free_indices[0]
    current_end = gap_free_indices[0]
    current_len = 1
    for idx in gap_free_indices[1:]:
        if idx <= current_end + allow_gap:
            current_end = idx
            current_len += 1
        else:
            if current_len > max_len:
                max_len = current_len
                best_start = current_start
                best_end = current_end
            current_start = idx
            current_end = idx
            current_len = 1
    # final check at end of loop
    if current_len > max_len:
        max_len = current_len
        best_start = current_start
        best_end = current_end
    # extract no-gap matrix
    no_gap_matrix = seq_matrix[:, best_start:best_end + 1]
    trimmed_len = best_end - best_start + 1
    print(f"Trimmed sequences to positions {best_start} to {best_start + max_len} (length {trimmed_len}) with few gaps. Gaps of size {allow_gap} were permitted.\n")
    # write trimmed sequences to new fasta
    with open(trimmed_fasta_path, 'w') as out_fasta:
        with open(trimmed_fasta_path.replace(".fasta", "_no_gaps.fasta"), 'w') as out_fasta_no_gaps:
            for i, seq_id in enumerate(ids):
                full_sequence = seq_matrix[i, :]
                start_in_sequence = 1
                for j in range(len(full_sequence)):
                    if full_sequence[j] != '-':
                        start_in_sequence += 1
                    if j == best_start - 1:
                        break
                end_in_sequence = start_in_sequence
                for j in range(best_start, best_end + 1):
                    if full_sequence[j] != '-':
                        end_in_sequence += 1
                    if j == best_end:
                        break
                trimmed_seq = ''.join(no_gap_matrix[i, :])
                no_gap_trimmed_seq = trimmed_seq.replace('-', '')
                out_fasta_no_gaps.write(f">{seq_id}_trimmed_{start_in_sequence}-{end_in_sequence}\n{no_gap_trimmed_seq}\n")
                out_fasta.write(f">{seq_id}_trimmed_{start_in_sequence}-{end_in_sequence}\n{trimmed_seq}\n")
                
def compute_pairwise_differences_between_sequences(df: pl.DataFrame) -> np.ndarray:
    ids = df["id"].to_list()
    sequences = df["sequence"].to_list()
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences[i]
            seq2 = sequences[j]
            shared_sequenced_length = sum((c1 not in ['-', 'N'] and c2 not in ['-', 'N']) for c1, c2 in zip(seq1, seq2))
            differences = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
            dist_matrix[i, j] = float(differences) / shared_sequenced_length if shared_sequenced_length > 0 else 1.0
            dist_matrix[j, i] = float(differences) / shared_sequenced_length if shared_sequenced_length > 0 else 1.0
    return ids, dist_matrix

def cluster_sequences_by_pid(ids: list, dist_matrix: np.ndarray, pid_threshold: float) -> pl.DataFrame:
    n = len(dist_matrix)
    assert dist_matrix.shape == (n, n)
    assert np.allclose(np.diag(dist_matrix), 0)
    dist_condensed = sch.distance.squareform(dist_matrix)
    Z = sch.linkage(dist_condensed, method='complete')
    cluster_labels = sch.fcluster(Z, t=1 - pid_threshold, criterion='distance')
    cluster_df = pl.DataFrame({
        "seq_id": ids,
        "cluster_id": cluster_labels.tolist()
    })
    return cluster_df

def identitify_best_representative_within_clusters(cluster_df: pl.DataFrame, dist_matrix: np.ndarray) -> pl.DataFrame:
    cluster_groups = np.unique(cluster_df["cluster_id"].to_numpy())
    best_reps = []
    for cluster in cluster_groups:
        members = cluster_df.filter(pl.col("cluster_id") == cluster)["seq_id"].to_list()
        if len(members) == 1:
            best_reps.append(members[0])
        elif len(members) == 2:
            # take first arbitrarily
            best_reps.append(members[0])
        else:
            member_indicies = [cluster_df['seq_id'].to_list().index(m) for m in members]
            sub_dist_matrix = dist_matrix[np.ix_(member_indicies, member_indicies)]
            avg_distances = sub_dist_matrix.mean(axis=1)
            best_index = np.argmin(avg_distances)
            best_reps.append(members[best_index])
    cluster_df = cluster_df.with_columns(
        pl.col("seq_id").is_in(best_reps).alias("is_best_representative")
    )
    return cluster_df

def process_fasta_and_cluster(fasta_path: str, pid_threshold: float) -> pl.DataFrame:
    df = parse_fasta_file(fasta_path)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(df)
    cluster_df = cluster_sequences_by_pid(ids, dist_matrix, pid_threshold)
    best_rep_df = identitify_best_representative_within_clusters(cluster_df, dist_matrix)
    print(f"Identified {best_rep_df.select(pl.col('cluster_id').n_unique()).item()} clusters for {len(df)} sequences at {pid_threshold * 100}% identity.\n")
    return best_rep_df

def generate_representative_targets_json(cluster_df: pl.DataFrame, output_json_path: str, json_key="targets") -> None:
    representative_targets = {}
    representative_targets['target_sets'] = {}
    representative_targets['target_sets'][json_key] = []
    for seq_id in cluster_df.filter(pl.col("is_best_representative") == True)["seq_id"].to_list():
        representative_targets['target_sets'][json_key].append(seq_id)
    with open(output_json_path, 'w') as json_file:
        json.dump(representative_targets, json_file, indent=4)

def remove_targets_by_identity(cluster_df: pl.DataFrame, cluster_threshold: float, fasta_path: str, target_to_remove: str, thresholds: np.ndarray, output_json_prefix: str) -> None:
    df = parse_fasta_file(fasta_path)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(df)
    thresholds = np.sort(thresholds) # ascending order
    target_to_remove_index = [id.split("_trimmed")[0] for id in ids].index(target_to_remove) # handle trimmed ids
    target_to_remove_full = ids[target_to_remove_index]
    for threshold in thresholds:
        to_remove = [target_to_remove_full]
        for i, id in enumerate(ids):
            if i != target_to_remove_index:
                if 1 - dist_matrix[target_to_remove_index, i] >= threshold:
                    to_remove.append(id)
        if len(to_remove) == len(ids):
            print(f"All other targets would be removed at >{int(threshold*100)}% identity with {target_to_remove}. Skipping this threshold.\n")
            continue
        if len(to_remove) == 1:
            print(f"No targets to remove at >{int(threshold*100)}% identity with {target_to_remove}.\n")
            continue
        to_remove_filter = np.isin(ids, to_remove)
        filtered_cluster_df = cluster_df.filter(to_remove_filter == False)
        filtered_cluster_df = filtered_cluster_df.drop("is_best_representative")
        filtered_dist_matrix = dist_matrix[to_remove_filter == False, :][:, to_remove_filter == False]
        assert filtered_dist_matrix.shape[0] == filtered_dist_matrix.shape[1] == len(filtered_cluster_df)
        output_json_path = f"{output_json_prefix}_targets_removed_at_gt_{threshold}_pid.json"
        filtered_cluster_df = identitify_best_representative_within_clusters(filtered_cluster_df, filtered_dist_matrix)
        generate_representative_targets_json(filtered_cluster_df.filter(pl.col("is_best_representative") == True), output_json_path, json_key=os.path.basename(output_json_path).replace(".json", ""))

def build_consensus_sequence_for_clusters(cluster_df: pl.DataFrame, pid_threshold: float, fasta_path: str, output_fasta_path: str, output_json_path: str, json_key="consensus_clusters", trimmed=False) -> None:
    sequences_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}
    cluster_groups = cluster_df.select(["cluster_id"]).unique().sort("cluster_id")["cluster_id"].to_list()
    cluster_ids_for_json = {}
    cluster_ids_for_json["target_sets"] = {}
    cluster_ids_for_json["target_sets"][json_key] = []
    with open(output_fasta_path, 'w') as out_fasta:
        for cluster in cluster_groups:
            members = cluster_df.filter(pl.col("cluster_id") == cluster)["seq_id"].to_list()
            cluster_id = give_cluster_a_useful_name(cluster, pid_threshold, cluster_df, trimmed=trimmed)
            if len(members) == 1:
                seq = sequences_dict[members[0]]
                seq_no_gaps = seq.replace('-', '')
                out_fasta.write(f">{cluster_id}\n{seq_no_gaps}\n")
            else:
                member_seqs = [sequences_dict[m] for m in members]
                seqs_matrix = np.array([list(seq) for seq in member_seqs])
                # remove cols with all gaps
                non_gap_cols = ~np.all(seqs_matrix == '-', axis=0)
                seqs_matrix = seqs_matrix[:, non_gap_cols]
                # build consensus
                consensus_seq = ""
                for pos in range(seqs_matrix.shape[1]):
                    bases_at_pos = seqs_matrix[:, pos]
                    unique_bases, counts = np.unique(bases_at_pos, return_counts=True)
                    filtered_bases = [(b, c) for b, c in zip(unique_bases, counts) if b not in ['-', 'N']] # remove gaps and N from consideration
                    if filtered_bases:
                        consensus_base = max(filtered_bases, key=lambda x: x[1])[0]
                    else:
                        consensus_base = 'N'  # if all are gaps or N
                    consensus_seq += consensus_base
                consensus_seq_no_gaps = consensus_seq.replace('-', '')
                out_fasta.write(f">{cluster_id}\n{consensus_seq_no_gaps}\n")
            cluster_ids_for_json["target_sets"][json_key].append(cluster_id)
    with open(output_json_path, 'w') as json_file:
        json.dump(cluster_ids_for_json, json_file, indent=4)

def give_cluster_a_useful_name(cluster_id: int, pid_threshold: float, cluster_df: pl.DataFrame, trimmed=False) -> str:
    members = cluster_df.filter(pl.col("cluster_id") == cluster_id)["seq_id"].to_list()
    full_species_in_cluster = [return_species_from_target_id(id_) for id_ in members]
    genera, species = zip(*(i.split(" ") for i in full_species_in_cluster))
    genera_in_cluster = list(set(genera))
    species_in_cluster = list(set(species))
    if trimmed:
        cluster_id = f"{'/'.join(genera_in_cluster)}_{'/'.join(species_in_cluster)}_pid_{pid_threshold}_trimmed_consensus_cluster{cluster_id}"
    else:
        cluster_id = f"{'/'.join(genera_in_cluster)}_{'/'.join(species_in_cluster)}_pid_{pid_threshold}_consensus_cluster{cluster_id}"
    return cluster_id

def trim_overhang(fasta_path: str, output_fasta_path: str, n_samples_to_distinguish_overhang: int = 2) -> None:
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        ids.append(record.id)
        sequences.append(list(str(record.seq)))
    seq_matrix = np.array(sequences)
    # identify columns where at least n_samples_to_distinguish_overhang have a gap
    potential_overhang_cols = np.sum(seq_matrix != '-', axis=0) <= n_samples_to_distinguish_overhang
    # trim columns from left until we hit a non-gap column
    left_trim_index = 0
    for i in range(seq_matrix.shape[1]):
        if potential_overhang_cols[i]:
            left_trim_index = i + 1
        else:
            break
    # trim columns from right until we hit a non-gap column
    right_trim_index = seq_matrix.shape[1]
    for i in range(seq_matrix.shape[1] - 1, -1, -1):
        if potential_overhang_cols[i]:
            right_trim_index = i
        else:
            break
    trimmed_seq_matrix = seq_matrix[:, left_trim_index:right_trim_index]
    with open(output_fasta_path, 'w') as out_fasta:
        for i, seq_id in enumerate(ids):
            trimmed_seq = ''.join(trimmed_seq_matrix[i, :])
            out_fasta.write(f">{seq_id}\n{trimmed_seq}\n")
