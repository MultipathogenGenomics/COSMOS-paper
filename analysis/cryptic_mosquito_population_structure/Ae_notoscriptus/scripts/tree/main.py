from analysis.Ae_notoscriptus.tree.process import build_Ae_notoscriptus_tree, build_Ae_notoscriptus_tree_including_frags, find_kmers_for_each_Ae_notoscriptus_cluster
import os

def main():
    fasta = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/Ae_notoscriptus/Ae_notoscriptus_mt_genomes.fasta"
    output_path = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/tree/Ae_notoscriptus.whole_mt_genome_tree.png"
    tmp_dir = "scratch/Ae_notoscriptus"
    # make sure paths exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    # analysis
    build_Ae_notoscriptus_tree(fasta, output_path, tmp_dir)
    # build tree including fragments
    fragments_path = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/Ae_notoscriptus/Ae_notoscriptus_fragments.Endersby_et_al.fasta"
    output_path_fragments = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/tree/Ae_notoscriptus.CO1_fragments_tree.png"
    build_Ae_notoscriptus_tree_including_frags(fasta, fragments_path, output_path_fragments, tmp_dir)
    # find kmers which define each cluster at 0.97 pid threshold
    cluster_info = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/tree/Ae_notoscriptus_fragments_and_refs_combined.COI_cluster_at_0.97.tsv"
    output_path = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Ae_notoscriptus/tree/Ae_notoscriptus.0.97_clusters.COI_unique_kmers.tsv"
    fragments_in_refs_path = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/Ae_notoscriptus/Ae_notoscriptus_fragments_in_references.fasta"
    find_kmers_for_each_Ae_notoscriptus_cluster(cluster_info, fragments_path, fasta, 31, output_path)

if __name__ == "__main__":
    main()