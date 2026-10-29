from analysis.common.BLAST_helpers import make_blast_db, blast_kmer, extract_best_blast_hit
import os
from Bio import SeqIO
import subprocess
import polars as pl
import numpy as np
from collections import Counter

def get_kmers_from_fasta(fasta_file: str, kmer_length: int) -> set:
    kmers = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i+kmer_length]
            kmers.append(kmer)
    return kmers

def find_Cx_annulirostris_frag_position_in_ref(fragment_fasta: str, reference_fasta: str, reference_contig: str, kmer_length: int, blast_output_dir: str) -> dict:
    # get start and end kmers from the fragment fasta
    fragment_kmers = get_kmers_from_fasta(fragment_fasta, kmer_length)
    start_kmer = fragment_kmers[0]
    end_kmer = fragment_kmers[-1]
    # Create BLAST database from reference fasta
    db = os.path.join(blast_output_dir, os.path.basename(reference_fasta))
    make_blast_db(reference_fasta, db)
    results = {}
    # Blast the start and end kmers against the reference database
    for skmer in [start_kmer, end_kmer]:
        output_blast_file = os.path.join(blast_output_dir, f"{skmer}_blast_output.txt")
        blast_kmer(skmer, db, output_blast_file)
        sseqid, best_start, best_end = extract_best_blast_hit(skmer, output_blast_file, reference_contig)
        label = "start" if skmer == start_kmer else "end"
        results[label] = {
            'kmer': skmer,
            'start': best_start,
            'end': best_end
        }
    return results

def extract_Cx_annulirostris_amplicon_region_from_reference(positions: dict, reference_fasta: str, reference_contig: str) -> str:
    start_pos = positions['start']['start'] # assuming the start kmer is on the forward strand
    end_pos = positions['end']['end']
    # extract_amplicon region from reference fasta
    for record in SeqIO.parse(reference_fasta, "fasta"):
        if record.id == reference_contig:
            amplicon_region = record.seq[start_pos - 1:end_pos] # -1 because fasta is 1-based and python is 0-based
            # check start with same letter as the start kmer and end with same letter as end kmer
            assert amplicon_region.startswith(positions['start']['kmer'][0]), f"Amplicon region in ref does not start with the same letter as the start of amplicon"
            assert amplicon_region.endswith(positions['end']['kmer'][-1]), f"Amplicon region in ref does not end with the same letter as the end of amplicon"
            return str(amplicon_region)
        
def msa_with_mafft(msa_input_fasta: str, msa_output_fasta: str) -> None:
    if os.path.exists(msa_output_fasta):
        print("WARNING: MSA Path exists, skipping MSA step. If you want to re-run MSA, delete the existing MSA output file.")
        return
    # run mafft to do msa
    cmd = f"mafft --auto {msa_input_fasta} > {msa_output_fasta}"
    subprocess.run(cmd, shell=True, check=True)
        
def msa_Cx_annulirostris_ref_and_amplicons(amplicon_seq: str, reference_contig: str, fragment_fasta: str, blast_output_dir: str) -> None:
    # write amplicon region to fasta file
    amplicon_fasta = os.path.join(blast_output_dir, "amplicon_region.fasta")
    with open(amplicon_fasta, "w") as f:
        f.write(f">{reference_contig}_amplicon_region\n")
        f.write(amplicon_seq)
    # concatenate amplicon fasta and fragment fasta into one fasta file for msa
    msa_input_fasta = os.path.join(blast_output_dir, "msa_input.fasta")
    with open(msa_input_fasta, "w") as f:
        for record in SeqIO.parse(fragment_fasta, "fasta"):
            SeqIO.write(record, f, "fasta")
        for record in SeqIO.parse(amplicon_fasta, "fasta"):
            SeqIO.write(record, f, "fasta")
    # run msa with mafft
    msa_output_fasta = os.path.join(blast_output_dir, "msa_output.fasta")
    msa_with_mafft(msa_input_fasta, msa_output_fasta)

def make_msa_seq_matrix(msa_fasta: str) -> tuple[np.ndarray, list]:
    sequences = []
    ids = []
    for record in SeqIO.parse(msa_fasta, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    seq_matrix = np.array([list(seq) for seq in sequences])
    return seq_matrix, ids

def parse_subtype_metadata(fragment_metadata: str, ids: list) -> pl.DataFrame:
    # parse and filter metadata
    df = pl.read_csv(fragment_metadata, separator='\t')
    df = df.filter(pl.col("PCR-RFLP (ITS1 rDNA)") == "Culex annulirostris") # only keep rows that are Culex annulirostris by PCR-RFLP
    df = df.filter(pl.col("Acc").is_not_null()) # only keep rows that have an Accession number
    df = df.with_columns(
        pl.col("mtDNA (COI)").str.split("annulirostris").list.get(1).str.strip_chars('"\n\t ').alias("subtype")
    )
    # add the fasta id to a column using the Acc column
    ids_dict = {}
    for genbank_acc in df["Acc"]:
        for fasta_id in ids:
            if genbank_acc in fasta_id:
                ids_dict[genbank_acc] = fasta_id
                break
    df = df.with_columns(
        pl.col("Acc").map_elements(lambda acc: ids_dict.get(acc, None)).alias("fasta_id")
    )
    return df

def find_differences_between_subtype_and_ref(ids: list, ref_seq_index: int, seq_matrix: np.ndarray) -> pl.DataFrame:
    differences = []
    ref_seq = seq_matrix[ref_seq_index, :]
    matrix_without_ref = np.delete(seq_matrix, ref_seq_index, axis=0)
    for pos in range(ref_seq.shape[0]):
        ref_base = ref_seq[pos]
        if ref_base == "-":
            continue # skip positions where the reference has a gap
        non_ref_bases = matrix_without_ref[:, pos]
        # if all non-ref bases are the same and different from the ref base, then this position is a difference between the subtype and the reference
        if len(set(non_ref_bases)) == 1 and non_ref_bases[0] != ref_base:
            row = {
                "position": pos,
                "ref_base": ref_base,
                "subtype_base": non_ref_bases[0]
            }
            differences.append(row)
    return pl.DataFrame(differences)

def remove_terminal_gaps_from_msa(seq_matrix: np.ndarray) -> np.ndarray:
    # where some sequences are shorter than others, trim down to length where all have sequence
    max_non_gap_pos = seq_matrix.shape[1]
    for seq in seq_matrix:
        while seq[max_non_gap_pos - 1] == "-":
            max_non_gap_pos -= 1
    return seq_matrix[:, :max_non_gap_pos]

def find_differences_between_subtype_and_all_others(subtype_ids: list, all_ids: list, ref_id: str, seq_matrix: np.ndarray) -> pl.DataFrame:
    differences = []
    subtype_indices = [all_ids.index(id) for id in subtype_ids]
    ref_index = all_ids.index(ref_id)
    non_subtype_indices = [i for i in range(len(all_ids)) if i not in subtype_indices and i != ref_index]
    # Count non-gap positions in reference as we go through the MSA
    ref_base_position = 0  # Position counter for non-gap bases in reference (0-based)
    for pos in range(seq_matrix.shape[1]):
        ref_base = seq_matrix[ref_index, pos]
        # Increment counter only for non-gap positions in reference
        if ref_base != "-":
            subtype_bases = set(seq_matrix[subtype_indices, pos])
            non_subtype_bases = set(seq_matrix[non_subtype_indices, pos])
            if (len(subtype_bases) == 1) and (list(subtype_bases)[0] not in list(non_subtype_bases)):
                row = {
                    "position": ref_base_position,  # Use the non-gap counter
                    "ref_base": ref_base,
                    "subtype_base": list(subtype_bases)[0],
                    "non_subtype_bases": ",".join(non_subtype_bases),
                    "n_subtype_seqs": len(subtype_indices),
                    "n_non_subtype_seqs": len(non_subtype_indices)
                }
                differences.append(row)
            ref_base_position += 1  # Increment after processing non-gap position
    return pl.DataFrame(differences)

def find_unique_kmers_for_subtype(
        metadata_df: pl.DataFrame,
        fragment_fasta: str,
        kmer_length: int,
        min_subtype_fraction: float = 1.0,
        min_subtype_count: int = None,
        max_other_subtype_count: int = 0
        ) -> pl.DataFrame:
    # Parse FASTA file once and create a dictionary mapping fasta_id to sequence
    fasta_dict = {}
    for record in SeqIO.parse(fragment_fasta, "fasta"):
        fasta_dict[record.id] = str(record.seq)
    # Get kmers for each fasta_id
    id_kmers = {}
    for fasta_id, sequence in fasta_dict.items():
        id_kmers[fasta_id] = set()
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i+kmer_length]
            id_kmers[fasta_id].add(kmer)
    # Get unique kmers for each subtype
    rows = []
    for subtype in metadata_df["subtype"].unique():
        fasta_ids = [
            fasta_id for fasta_id in
            metadata_df.filter(pl.col("subtype") == subtype)["fasta_id"].to_list()
            if fasta_id in id_kmers
        ]
        other_fasta_ids = [
            fasta_id for fasta_id in
            metadata_df.filter(pl.col("subtype") != subtype)["fasta_id"].to_list()
            if fasta_id in id_kmers
        ]
        if not fasta_ids:
            continue

        subtype_kmer_counts = Counter()
        for fasta_id in fasta_ids:
            subtype_kmer_counts.update(id_kmers[fasta_id])

        other_kmer_counts = Counter()
        for fasta_id in other_fasta_ids:
            other_kmer_counts.update(id_kmers[fasta_id])

        min_required_count = (
            min_subtype_count
            if min_subtype_count is not None
            else int(np.ceil(min_subtype_fraction * len(fasta_ids)))
        )

        n_subtype_kmers = 0
        for kmer, n_subtype_seqs_with_kmer in subtype_kmer_counts.items():
            n_other_seqs_with_kmer = other_kmer_counts.get(kmer, 0)
            if (
                n_subtype_seqs_with_kmer >= min_required_count
                and n_other_seqs_with_kmer <= max_other_subtype_count
            ):
                n_subtype_kmers += 1
                rows.append({
                    "subtype": subtype,
                    "kmer": kmer,
                    "n_subtype_seqs_with_kmer": n_subtype_seqs_with_kmer,
                    "n_subtype_seqs": len(fasta_ids),
                    "subtype_fraction": n_subtype_seqs_with_kmer / len(fasta_ids),
                    "n_other_seqs_with_kmer": n_other_seqs_with_kmer,
                    "n_other_seqs": len(other_fasta_ids),
                    "min_subtype_fraction": min_subtype_fraction,
                    "min_subtype_count": min_required_count,
                    "max_other_subtype_count": max_other_subtype_count
                })
        print(
            f"Found {n_subtype_kmers} unique kmers for subtype {subtype} "
            f"with min_subtype_fraction={min_subtype_fraction}, "
            f"min_subtype_count={min_required_count}, "
            f"max_other_subtype_count={max_other_subtype_count}"
        )
    return pl.DataFrame(
        rows,
        schema={
            "subtype": pl.String,
            "kmer": pl.String,
            "n_subtype_seqs_with_kmer": pl.Int64,
            "n_subtype_seqs": pl.Int64,
            "subtype_fraction": pl.Float64,
            "n_other_seqs_with_kmer": pl.Int64,
            "n_other_seqs": pl.Int64,
            "min_subtype_fraction": pl.Float64,
            "min_subtype_count": pl.Int64,
            "max_other_subtype_count": pl.Int64
        }
    )

def find_Cx_annulirostris_subtype_defining_snps(fragment_fasta: str, fragment_metadata: str, reference_fasta: str, reference_contig: str, kmer_length: int, tmp_output_dir: str, ignore_subtypes: list = None) -> None:
    # find position of fragment in reference
    positions = find_Cx_annulirostris_frag_position_in_ref(fragment_fasta, reference_fasta, reference_contig, kmer_length, tmp_output_dir)
    start = positions['start']['start'] # assuming the start kmer is on the forward strand
    end = positions['end']['end']
    # find the sequence of the amplicon region in the reference fasta
    amplicon_seq = extract_Cx_annulirostris_amplicon_region_from_reference(positions, reference_fasta, reference_contig)
    # run multiple sequence alignment with amplicon fragments from Hemmerter et al. - MSA will be in tmp_output_dir/msa_output.fasta
    msa_Cx_annulirostris_ref_and_amplicons(amplicon_seq, reference_contig, fragment_fasta, tmp_output_dir)
    # make a sequence matrix from the MSA output fasta
    msa_output_fasta = os.path.join(tmp_output_dir, "msa_output.fasta")
    # make a sequence matrix from the MSA output fasta and get the corresponding ids for each row of the sequence matrix
    seq_matrix, ids = make_msa_seq_matrix(msa_output_fasta)
    # parse subtypes fragment metadata to get subtype information for each fragment
    fragment_metadata = parse_subtype_metadata(fragment_metadata, ids)
    # remove subtypes to ignore from the metadata table and the sequence matrix and ids list
    if ignore_subtypes is not None:
        fragment_metadata = fragment_metadata.filter(~pl.col("subtype").is_in(ignore_subtypes))
        ids_to_keep = fragment_metadata["fasta_id"].to_list() + [f"{reference_contig}_amplicon_region"] # make sure to keep the reference sequence in the MSA
        indices_to_keep = [ids.index(id) for id in ids if id in ids_to_keep]
        seq_matrix = seq_matrix[indices_to_keep, :]
        ids = [ids[i] for i in indices_to_keep]
    # trim terminal gaps from the MSA sequence matrix
    previous_shape = seq_matrix.shape
    seq_matrix = remove_terminal_gaps_from_msa(seq_matrix)
    new_shape = seq_matrix.shape
    print(f"Trimmed amplicons to length {new_shape[1]} from original length {previous_shape[1]} to remove terminal gaps")
    # final dataframe to hold subtype defining SNPs
    results = None
    for subtype in fragment_metadata["subtype"].unique():
        # get genbank ids from the metadata table
        subtype_genbank_ids = fragment_metadata.filter(pl.col("subtype") == subtype)["Acc"].to_list()
        # extract the rows from the sequence matrix that correspond to the genbank ids for this subtype
        subtype_indices = [i for i, fasta_id in enumerate(ids) if any(genbank_id in fasta_id for genbank_id in subtype_genbank_ids)]
        subtype_ids = [ids[i] for i in subtype_indices]
        # find positions where all sequences for this subtype are the same and different from the other sequences.
        differences_df = find_differences_between_subtype_and_all_others(subtype_ids, ids, f"{reference_contig}_amplicon_region", seq_matrix)
        if differences_df.shape[0] == 0:
            continue # if there are no differences between this subtype and all others, move on
        differences_df = differences_df.with_columns(
            pl.lit(subtype).alias("subtype"),
            (pl.col("position") + start).alias("position_in_reference")
            )
        # concatenate the differences dataframe for this subtype with the results dataframe for all subtypes
        if results is None:
            results = differences_df
        else:
            results = pl.concat([results, differences_df], how="vertical")
    return results

def get_ids_from_fasta(fasta_file: str) -> list:
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids.append(record.id)
    return ids

def find_unique_kmers_between_Cx_annulirostris_subtypes(
        fragment_metadata: str,
        fragment_fasta: str,
        kmer_length: int,
        min_subtype_fraction: float = 1.0,
        min_subtype_count: int = None,
        max_other_subtype_count: int = 0
        ) -> pl.DataFrame:
    ids = get_ids_from_fasta(fragment_fasta)
    metadata_df = parse_subtype_metadata(fragment_metadata, ids)
    unique_kmers_df = find_unique_kmers_for_subtype(
        metadata_df,
        fragment_fasta,
        kmer_length,
        min_subtype_fraction=min_subtype_fraction,
        min_subtype_count=min_subtype_count,
        max_other_subtype_count=max_other_subtype_count
    )
    return unique_kmers_df
