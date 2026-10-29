from analysis.common.paths import RESULTS_DIR
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform
from analysis.common.mosquito_ids import return_species_from_target_id
from Bio import SeqIO


# Helpers retained from the original dynamics_poster modules.
def short_genus(species_name: str) -> str:
    if "Culex" in species_name:
        return species_name.replace("Culex", "Cx.")
    elif "Anopheles" in species_name:
        return species_name.replace("Anopheles", "An.")
    elif "Aedes" in species_name:
        return species_name.replace("Aedes", "Ae.")
    elif "Mansonia" in species_name:
        return species_name.replace("Mansonia", "Ma.")
    elif "Ochlerotatus" in species_name:
        return species_name.replace("Ochlerotatus", "Ae.")
    elif "Haemagogus" in species_name:
        return species_name.replace("Haemagogus", "Hg.")
    elif "Toxorhynchites" in species_name:
        return species_name.replace("Toxorhynchites", "Tx.")
    elif "Coquillettidia" in species_name:
        return species_name.replace("Coquillettidia", "Cq.")
    elif "Psorophora" in species_name:
        return species_name.replace("Psorophora", "Ps.")
    elif "Sabethes" in species_name:
        return species_name.replace("Sabethes", "Sa.")
    elif "Tripteroides" in species_name:
        return species_name.replace("Tripteroides", "Tr.")
    elif "Uranotaenia" in species_name:
        return species_name.replace("Uranotaenia", "Ur.")
    elif "Verallina" in species_name:
        return species_name.replace("Verallina", "Ve.")
    else:
        return species_name

def parse_fasta_file(fasta_file: str) -> dict:
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def get_kmers_from_sequence(seq: str, kmer_size: int) -> set:
    seq_upper = seq.upper()
    return set(seq_upper[i:i+kmer_size] for i in range(len(seq_upper) - kmer_size + 1))

def compute_jaccard_index(kmers1: set, kmers2: set) -> float:
    intersection = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    return intersection / union if union > 0 else 0.0

def compute_mash_distance(jaccard_index: float, kmer_size: int) -> float:
    if jaccard_index <= 0:
        return float(1.0)  # Maximum distance if no shared kmers
    elif jaccard_index >= 1:
        return float(0.0)  # Minimum distance if all kmers are shared
    else:
        mash_distance = -1.0 / kmer_size * np.log(2 * jaccard_index / (1 + jaccard_index))
        return mash_distance

def compute_mash_matrix(sequences: dict, kmer_size=31):
    matrix = np.zeros((len(sequences), len(sequences)))
    for i, (id1, seq1) in enumerate(sequences.items()):
        kmers1 = get_kmers_from_sequence(seq1, kmer_size)
        jaccard_idx = compute_jaccard_index(kmers1, kmers1)  # Jaccard index with itself should be 1
        matrix[i, i] = compute_mash_distance(jaccard_idx, kmer_size)
        for j, (id2, seq2) in enumerate(sequences.items()):
            if i < j:  # Compute only for upper triangle (matrix is symmetric)
                kmers2 = get_kmers_from_sequence(seq2, kmer_size)
                jaccard_idx = compute_jaccard_index(kmers1, kmers2)
                mash_dist = compute_mash_distance(jaccard_idx, kmer_size)
                matrix[i, j] = mash_dist
                matrix[j, i] = mash_dist  # Symmetric entry
    return matrix, list(sequences.keys())

def draw_mash_tree(matrix, labels, title=None):
    distance_matrix = np.asarray(matrix, dtype=float)
    condensed_distance = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_distance, method="complete", optimal_ordering=True) # check if should be complete
    # draw plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 14))
    dendrogram(
        linkage_matrix,
        labels=labels,
        ax=ax,
        leaf_rotation=0,
        leaf_font_size=8,
        orientation='left',
        above_threshold_color='k',
        color_threshold=0.0,
    )
    # Use symlog scale to spread out values near 0 for better visibility
    ax.set_xscale('symlog', linthresh=0.01)
    # Set custom tick labels in normal decimal notation instead of scientific
    tick_positions = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    tick_labels = ['0', '0.001', '0.01', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '1']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    if title is not None:
        ax.set_title(title, fontsize=12, weight="semibold")
    # if first axis, make the leaf labels italic
    ax.set_yticklabels(ax.get_yticklabels(), fontstyle="italic")
    ax.grid(axis="x", linestyle="--", linewidth=0.6, alpha=0.45)
    sns.despine(ax=ax, left=True, bottom=True, right=True, top=True)
    # make leaf labels large
    ax.tick_params(axis="y", labelsize=14)
    ax.set_xlabel("MASH Distance", fontsize=18)
    return fig, ax

def main():
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description="Draw a Mash-distance tree from a reference FASTA")
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=RESULTS_DIR / 'cluster_targets_by_pid/mash/reference_mash_tree.png')
    args = parser.parse_args()
    seqs = parse_fasta_file(args.fasta)
    matrix, labels = compute_mash_matrix(seqs)
    mosquito_species_labels = [return_species_from_target_id(i) for i in labels]
    for i, (label, sp) in enumerate(zip(labels, mosquito_species_labels)):
        if "pipiens" in sp:
            mosquito_species_labels[i] = sp + (" molestus" if "molestus" in label else " pipiens")
    fig, ax = draw_mash_tree(matrix, [short_genus(i) for i in mosquito_species_labels])
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
