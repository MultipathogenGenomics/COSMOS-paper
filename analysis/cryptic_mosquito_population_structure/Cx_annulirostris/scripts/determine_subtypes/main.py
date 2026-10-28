from analysis.Cx_annulirostris.determine_subtypes.process import find_Cx_annulirostris_subtype_defining_snps, find_unique_kmers_between_Cx_annulirostris_subtypes
import os
import sys

def main():

    # input files
    fragment_fasta = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/Cx_annulirostris_Hemmerter_et_al/Cx_annulirostris_diversity.fasta"
    fragment_metadata = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/Cx_annulirostris_Hemmerter_et_al/Cx_annulirostris_fragment_metadata.tsv"
    reference_fasta = "/Users/afor0966/Documents/code/multipathogen_genomics/data/targets/mozz_2024-11-01_min1kb_AF_trial_update_plus_additional_COX1s.fasta"
    reference_contig = "mosquito_NC_054323.1-1430-2966-Culex-annulirostris-isolate-W01-mitochondrion"
    kmer_length = 31
    tmp_output_dir = "/Users/afor0966/Documents/code/multipathogen_genomics/scratch/Cx_annulirostris"
    outputs_dir = "/Users/afor0966/Documents/code/multipathogen_genomics/results/Cx_annulirostris"
    all_subtypes = ["ann PNG-1", "ann PNG-2", "ann S-AUS", "ann SI"]

    # make sure tmp output dir exists
    os.makedirs(tmp_output_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)

    # find subtype defining kmers
    unique_kmers_df = find_unique_kmers_between_Cx_annulirostris_subtypes(fragment_metadata, fragment_fasta, kmer_length)
    output_file =  os.path.join(outputs_dir, "Cx_annulirostris.subtype_unique_kmers.tsv")
    unique_kmers_df.select(["subtype", "kmer"]).write_csv(output_file, separator='\t')

    # find subtype private kmers which are present in at least some members of a subtype,
    # while remaining absent from all other subtypes
    min_subtype_fraction = 0.5
    private_kmers_df = find_unique_kmers_between_Cx_annulirostris_subtypes(
        fragment_metadata,
        fragment_fasta,
        kmer_length,
        min_subtype_fraction=min_subtype_fraction
    )
    output_file = os.path.join(
        outputs_dir,
        f"Cx_annulirostris.subtype_private_kmers.min_fraction_{min_subtype_fraction}.tsv"
    )
    private_kmers_df.write_csv(output_file, separator='\t')
    sys.exit()

    # find subtype defining SNPs
    output_file = os.path.join(outputs_dir, "Cx_annulirostris.SNPs_between_subtype_and_other_fragments.tsv")
    results = find_Cx_annulirostris_subtype_defining_snps(
        fragment_fasta, fragment_metadata, reference_fasta, reference_contig, kmer_length, tmp_output_dir
    )
    # save results to output file
    results.write_csv(output_file, separator='\t')
    # version just comparing AUS and S-AUS
    ignore_subtypes = ["ann PNG-1", "ann PNG-2", "ann SI"]
    results_ignore_subtypes = find_Cx_annulirostris_subtype_defining_snps(
        fragment_fasta, fragment_metadata, reference_fasta, reference_contig, kmer_length, tmp_output_dir, ignore_subtypes=ignore_subtypes
    )
    output_file_ignore_subtypes = os.path.join(outputs_dir, "Cx_annulirostris.SNPs_between_AUS_and_S-AUS.tsv")
    results_ignore_subtypes.write_csv(output_file_ignore_subtypes, separator='\t')
    # version with all but SI
    ignore_subtypes = ["ann SI"]
    results_ignore_subtypes = find_Cx_annulirostris_subtype_defining_snps(
        fragment_fasta, fragment_metadata, reference_fasta, reference_contig, kmer_length, tmp_output_dir, ignore_subtypes=ignore_subtypes
    )
    output_file_ignore_subtypes = os.path.join(outputs_dir, "Cx_annulirostris.SNPs_between_all_excluding_SI.tsv")
    results_ignore_subtypes.write_csv(output_file_ignore_subtypes, separator='\t')
    # version with AUS and S-AUS + each individually
    for subtype in ["ann PNG-1", "ann PNG-2", "ann SI"]:
        remaining_to_ignore = [s for s in all_subtypes if s not in [subtype, "ann S-AUS", "ann AUS"]]
        results_ignore_subtypes = find_Cx_annulirostris_subtype_defining_snps(
            fragment_fasta, fragment_metadata, reference_fasta, reference_contig, kmer_length, tmp_output_dir, ignore_subtypes=remaining_to_ignore
        )
        output_file_ignore_subtypes = os.path.join(outputs_dir, f"Cx_annulirostris.SNPs_between_AUS_S-AUS_and_{subtype.replace('ann ', '')}.tsv")
        results_ignore_subtypes.write_csv(output_file_ignore_subtypes, separator='\t')

if __name__ == "__main__":
    main()
