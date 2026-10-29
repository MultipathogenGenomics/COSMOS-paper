from analysis.common.paths import DATA_DIR, RESULTS_DIR
from analysis.cryptic_mosquito_population_structure.Ae_notoscriptus.scripts.subtypes.process import find_Ae_notoscriptus_cluster_defining_snps, merge_cluster_defining_snps
import os

def main():
    msa_path = f"{RESULTS_DIR}/intermediates/Ae_notoscriptus/Ae_notoscriptus_fragments_combined.all.msa.fasta"
    fragment_metadata = f"{RESULTS_DIR}/Ae_notoscriptus/tree/Ae_notoscriptus_fragments_and_refs_combined.COI_cluster_at_0.97.tsv"
    reference_path = f"{DATA_DIR}/references/lineages/KM676218.1_Aedes_notoscriptus_MC04_COX1.fasta"
    reference_contig = "KM676218.1_Aedes_notoscriptus_isolate_MC04_mitochondrion_complete_genome"
    tmp_dir = f"{RESULTS_DIR}/intermediates/Ae_notoscriptus"
    ignore_clusters = None
    output_path = f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps.tsv"
    # make sure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    results_df = find_Ae_notoscriptus_cluster_defining_snps(msa_path, fragment_metadata, reference_path, reference_contig, tmp_dir, ignore_clusters)
    results_df.write_csv(output_path, separator='\t')
    # get SNPs for 1 vs 3 and then 1 vs 2 by ignoring the other cluster in each case
    output_path = f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_3.tsv"
    ignore_clusters = [2]
    results_df = find_Ae_notoscriptus_cluster_defining_snps(msa_path, fragment_metadata, reference_path, reference_contig, tmp_dir, ignore_clusters)
    results_df.write_csv(output_path, separator='\t')
    output_path = f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_2.tsv"
    ignore_clusters = [3]
    results_df = find_Ae_notoscriptus_cluster_defining_snps(msa_path, fragment_metadata, reference_path, reference_contig, tmp_dir, ignore_clusters)
    results_df.write_csv(output_path, separator='\t')
    # also do 2 vs 3
    output_path = f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_2_vs_3.tsv"
    ignore_clusters = [1]
    results_df = find_Ae_notoscriptus_cluster_defining_snps(msa_path, fragment_metadata, reference_path, reference_contig, tmp_dir, ignore_clusters)
    results_df.write_csv(output_path, separator='\t')

    # merge to get 1 vs 2 or 3
    merge_cluster_defining_snps(
        files_to_merge=[
            f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_3.tsv",
            f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_2.tsv"
        ],
        focal_subtype=1,
        output_file=f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_1_vs_2_or_3.tsv"
    )

    # run on 2 vs 3 as well to fix bases
    merge_cluster_defining_snps(
        files_to_merge=[
            f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_2_vs_3.tsv"
        ],
        focal_subtype=2,
        output_file=f"{RESULTS_DIR}/Ae_notoscriptus/subtypes/Ae_notoscriptus_cluster_defining_snps_2_vs_3_fixed.tsv"
    )


if __name__ == "__main__":
    main()