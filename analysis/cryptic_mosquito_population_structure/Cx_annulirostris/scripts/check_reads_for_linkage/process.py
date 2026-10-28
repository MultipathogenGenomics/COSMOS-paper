import pysam, os
import polars as pl

def check_if_reads_support_subtypes_by_linkage(snp_ratio_file: str, bam_dir: str, abundance_results: str, focal_subtype: str, contig: str, 
                                               mapq_threshold: int=20, sample_list: list=[], max_check: int=500, chunk: int=None, total_chunks: int=None) -> None:
    snp_df = pl.read_csv(snp_ratio_file, separator='\t')
    if sample_list == []:
        samples = snp_df['sample'].unique().to_list()
    else:
        samples = [i for i in snp_df['sample'].unique().to_list() if i in sample_list ]
    results = []
    if chunk is not None:
        chunk_size = len(samples) // total_chunks + 1
        samples = samples[chunk*chunk_size:(chunk+1)*chunk_size]
    for sample in samples:
        print("Assessing sample:", sample)
        bam_file = os.path.join(bam_dir, sample, f"{sample}_kma_mapped_to_genome.sorted.bam")
        bam = pysam.AlignmentFile(bam_file, "rb")
        positions = snp_df.filter(pl.col("sample") == sample)['position_in_reference'].to_list()
        subtype_bases = snp_df.filter(pl.col("sample") == sample)[f'{focal_subtype}_subtype_bases'].to_list()
        non_subtype_bases = snp_df.filter(pl.col("sample") == sample)[f'{focal_subtype}_non_subtype_bases'].to_list()
        assessed_reads = []
        rows = []
        for i, pos in enumerate(positions[:-1]):
            # gather all the bases at this position for this sample
            A_bases = [base.upper() for base in subtype_bases[i].split(",")]
            a_bases = [base.upper() for base in non_subtype_bases[i].split(",")]
            next_pos = positions[i+1]
            B_bases = [base.upper() for base in subtype_bases[i+1].split(",")]
            b_bases = [base.upper() for base in non_subtype_bases[i+1].split(",")]
            # counts to store
            AB = 0
            Ab = 0
            aB = 0
            ab = 0
            # now let's iterate over reads at pos
            pos_0based = pos - 1
            next_pos_0based = next_pos - 1
            start = min(pos_0based, next_pos_0based)
            end = max(pos_0based, next_pos_0based) + 1
            reads_checked = 0
            for read in bam.fetch(contig, start=start, end=end):
                if reads_checked >= max_check:
                    print(f"Max number of reads checked (position={pos}, N={max_check}), moving to next SNP pair...")
                    break
                if read.query_name in assessed_reads:
                    continue
                try:
                    mate = bam.mate(read)
                except ValueError:
                    continue
                if read.is_paired and not read.mate_is_unmapped and not read.is_duplicate and read.mapping_quality >= mapq_threshold and mate.mapping_quality >= mapq_threshold:
                    # find base at pos and next_pos for this read
                    read_pos = read.get_aligned_pairs()
                    read_pos_mate = mate.get_aligned_pairs()
                    read_ref_to_query = {ref: query for query, ref in read_pos if ref is not None and query is not None}
                    mate_ref_to_query = {ref: query for query, ref in read_pos_mate if ref is not None and query is not None}
                    all_ref_pos = list(read_ref_to_query.keys()) + list(mate_ref_to_query.keys())
                    if pos_0based not in all_ref_pos or next_pos_0based not in all_ref_pos:
                        continue
                    # get base of pos in read
                    if pos_0based in read_ref_to_query:
                        pos_in_read = read_ref_to_query[pos_0based]
                        base_at_pos = read.query_sequence[pos_in_read].upper()
                    else:
                        pos_in_read = mate_ref_to_query.get(pos_0based)
                        base_at_pos = mate.query_sequence[pos_in_read].upper()
                    # get base of next_pos in read
                    if next_pos_0based in read_ref_to_query:
                        pos_in_read = read_ref_to_query[next_pos_0based]
                        base_at_next_pos = read.query_sequence[pos_in_read].upper()
                    else:
                        pos_in_read = mate_ref_to_query.get(next_pos_0based)
                        base_at_next_pos = mate.query_sequence[pos_in_read].upper()
                    # now check if base at pos is in A
                    if base_at_pos in A_bases and base_at_next_pos in B_bases:
                        AB += 1
                    elif base_at_pos in A_bases and base_at_next_pos in b_bases:
                        Ab += 1
                    elif base_at_pos in a_bases and base_at_next_pos in B_bases:
                        aB += 1
                    elif base_at_pos in a_bases and base_at_next_pos in b_bases:
                        ab += 1
                    assessed_reads.append(read.query_name)
                reads_checked += 1
            rows.append({
                "sample": sample,
                "locusA": pos,
                "baseA": ",".join(A_bases),
                "basea": ",".join(a_bases),
                "locusB": next_pos,
                "baseB": ",".join(B_bases),
                "baseb": ",".join(b_bases),
                "AB": AB,
                "Ab": Ab,
                "aB": aB,
                "ab": ab
            })
        results.extend(rows)
    print("Finito!")
    return pl.DataFrame(results)