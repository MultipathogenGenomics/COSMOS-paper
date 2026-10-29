from matplotlib import pyplot as plt
import seaborn as sns
import polars as pl
import numpy as np
import scipy.cluster.hierarchy as sch
from analysis.common.mosquito_ids import shortened_mosquito_target_id, return_species_from_target_id
from analysis.mosquito_community_composition.scripts.cluster_targets_by_pid.process import compute_pairwise_differences_between_sequences, parse_fasta_file

def plot_tree(distance_matrix: np.ndarray, ids: list, pid_threshold: float, output_path: str, colour_leaves: list = None, bold_leaves: list = None, shade_out_leaves: list = None, annotate_pid_to_target: str = None) -> None:
    dist_condensed = sch.distance.squareform(distance_matrix)
    Z = sch.linkage(dist_condensed, method='complete')
    fig, ax = plt.subplots(figsize=(8, max(6, len(ids) * 0.1)))
    tree_labels = [shortened_mosquito_target_id(id_) for id_ in ids]
    dn = sch.dendrogram(Z, labels=tree_labels, leaf_rotation=0, orientation='left', ax=ax)
    if colour_leaves:
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            if lbl.get_text() in [shortened_mosquito_target_id(id_) for id_ in colour_leaves]:
                lbl.set_color('tab:red')
    if bold_leaves:
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            if lbl.get_text() in [shortened_mosquito_target_id(id_) for id_ in bold_leaves]:
                lbl.set_fontweight('bold')
    if shade_out_leaves:
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            if lbl.get_text() in [shortened_mosquito_target_id(id_) for id_ in shade_out_leaves]:
                lbl.set_color('lightgray')
    if annotate_pid_to_target:
        idx_dict = {shortened_mosquito_target_id(id_): idx for idx, id_ in enumerate(ids)}
        focal_idx = idx_dict[shortened_mosquito_target_id(annotate_pid_to_target)]
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            y = lbl.get_position()[1]
            target_id = lbl.get_text()
            target_idx = idx_dict.get(target_id, None)
            pid_value = 1 - distance_matrix[focal_idx, target_idx]
            ax.text(1, y, f"{pid_value:.2f}", transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize='x-small', color='tab:blue')
            ax.tick_params(axis='y', which='both', pad=15)
    if pid_threshold > 0:
        ax.axvline(x=1 - pid_threshold, color='k', linestyle='--', label=f'pid threshold: {pid_threshold}')
        ax.legend(frameon=False, loc='upper left', fontsize='small')
    ax.set_xlabel('')
    ax.set_ylabel('distance (1 - pid)')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.tick_params(axis='x', which='both', length=0)
    fig.tight_layout()
    sns.despine(left=True, bottom=True, right=True, top=True, ax=ax)
    fig.savefig(output_path)
    plt.close(fig)

def plot_distribution(input_arr: np.ndarray, ax: plt.Axes, title: str = None) -> plt.Axes:
    sns.histplot(input_arr, bins=50, kde=True, ax=ax, color='tab:blue')
    if title:
        ax.set_title(title)
    return ax

def plot_pid_distributions(dist_matrix: np.ndarray, ids: list, output_path: str) -> None:
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
    ax = ax.flatten()
    species = [return_species_from_target_id(id_) for id_ in ids]
    # plot within species pid distribution
    within_species_pids = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if species[i] == species[j]:
                pid = 1 - dist_matrix[i, j]
                within_species_pids.append(pid)
    plot_distribution(np.array(within_species_pids), ax[0], title='within species')
    # plot between species pid distribution
    between_species_pids = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if species[i] != species[j]:
                pid = 1 - dist_matrix[i, j]
                between_species_pids.append(pid)
    plot_distribution(np.array(between_species_pids), ax[1], title='between species')
    # plot all pid distribution
    all_pids = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            pid = 1 - dist_matrix[i, j]
            all_pids.append(pid)
    plot_distribution(np.array(all_pids), ax[2], title='all pairs')
    # plot pid with nearest non-species neighbour
    nearest_non_species_pids = []
    for specie in set(species):
        is_species = np.isin(np.array(species), specie)
        min_value = np.min(dist_matrix[is_species][:, ~is_species])
        nearest_non_species_pids.append(1 - min_value)
    plot_distribution(np.array(nearest_non_species_pids), ax[3], title='nearest non-species neighbour')
    sns.despine(fig=fig, top=True, right=True, trim=True)
    fig.tight_layout()
    fig.savefig(output_path)

def build_tree_plot_from_fasta(df: pl.DataFrame, fasta_path: str, pid_threshold: float, output_path: str) -> None:
    sequence_df = parse_fasta_file(fasta_path)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(sequence_df)
    plot_tree(dist_matrix, ids, pid_threshold, output_path, bold_leaves=df.filter(pl.col("is_best_representative") == True)["seq_id"].to_list())

def build_pid_distribution_plots_from_fasta(fasta_path: str, output_path: str) -> None:
    sequence_df = parse_fasta_file(fasta_path)
    ids, dist_matrix = compute_pairwise_differences_between_sequences(sequence_df)
    plot_pid_distributions(dist_matrix, ids, output_path)