from Bio import SeqIO
from Bio.Seq import Seq
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import os

def parse_alignment_file(alignment_file: str, keep_sequences: list = None) -> dict:
    sequences = {}
    for record in SeqIO.parse(alignment_file, "fasta"):
        if keep_sequences is None or record.id in keep_sequences:
            sequences[record.id] = str(record.seq)
    return sequences

def cut_to_segments(sequences: dict, ref_seq_id: str, deliniation_kmers: dict) -> dict:
    ref_seq = sequences[ref_seq_id]
    segments = {}
    for segment_name, (start_kmer, end_kmer) in deliniation_kmers.items():
        start_index = ref_seq.find(start_kmer)
        end_index = ref_seq.find(end_kmer)
        if start_index == -1 or end_index == -1:
            raise ValueError(f"Could not find k-mers for segment {segment_name}")
        segments[segment_name] = {seq_id: seq[start_index:end_index + len(end_kmer)] for seq_id, seq in sequences.items()}
    return segments

def remove_gaps_where_all_gaps(sequences: dict) -> dict:
    # Convert sequences to a list of lists for easier manipulation
    seq_list = [list(seq) for seq in sequences.values()]
    # Transpose to get columns
    transposed = list(zip(*seq_list))
    # Identify columns where all sequences are missing. Treat N like a gap here
    # because it represents an uncalled base rather than sequence coverage.
    columns_to_keep = [i for i, col in enumerate(transposed) if not all(base.upper() in ["-", "N"] for base in col)]
    # Create new sequences without the all-gap columns
    new_sequences = {seq_id: ''.join(seq[i] for i in columns_to_keep) for seq_id, seq in sequences.items()}
    return new_sequences

def map_ref_positions_to_alignment_positions(ref_seq: str) -> dict:
    ref_pos_to_aln_pos = {}
    ref_pos = 0
    for aln_pos, base in enumerate(ref_seq):
        if base.upper() in ["A", "C", "G", "T"]:
            ref_pos += 1
            ref_pos_to_aln_pos[ref_pos] = aln_pos
    return ref_pos_to_aln_pos

def get_nearest_alignment_position(ref_pos_to_aln_pos: dict, ref_pos: int) -> int:
    if ref_pos in ref_pos_to_aln_pos:
        return ref_pos_to_aln_pos[ref_pos]
    available_positions = sorted(ref_pos_to_aln_pos.keys())
    nearest_ref_pos = min(available_positions, key=lambda x: abs(x - ref_pos))
    return ref_pos_to_aln_pos[nearest_ref_pos]

def draw_gene_annotations(ax, segment: str, segment_seq: str, genes: dict) -> None:
    segment_genes = {gene: info for gene, info in genes.items() if info["segment"] == segment}
    ref_pos_to_aln_pos = map_ref_positions_to_alignment_positions(segment_seq)
    ymax = max(0.42, 0.16 * len(segment_genes) + 0.16)

    for gene_idx, (gene, gene_info) in enumerate(segment_genes.items()):
        start_pos = get_nearest_alignment_position(ref_pos_to_aln_pos, gene_info["start_pos"])
        end_pos = get_nearest_alignment_position(ref_pos_to_aln_pos, gene_info["end_pos"]) + 1
        y_pos = 0.08 + gene_idx * 0.16
        ax.annotate(
            "",
            xy=(end_pos, y_pos),
            xytext=(start_pos, y_pos),
            arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 0.8, "mutation_scale": 8},
        )
        ax.text(
            (start_pos + end_pos) / 2,
            y_pos + 0.025,
            gene.replace("_", " "),
            ha="center",
            va="bottom",
            fontsize=6,
        )

    ax.set_xlim(0, len(segment_seq))
    ax.set_ylim(0, ymax)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

def translate_codon(codon: str):
    codon = codon.upper()
    if len(codon) != 3 or any(base not in ["A", "C", "G", "T"] for base in codon):
        return None
    return str(Seq(codon).translate())

def get_codon_at_alignment_positions(seq: str, aln_positions: list) -> str:
    return "".join(seq[aln_pos] for aln_pos in aln_positions)

def make_aa_difference_matrix_and_summary(segments: dict, ref_seq_id: str, genes: dict) -> tuple:
    aa_difference_matrix = {
        segment: {seq_id: [0] * len(next(iter(segment_seqs.values()))) for seq_id in segment_seqs.keys()}
        for segment, segment_seqs in segments.items()
    }
    gene_summary_rows = []

    for gene, gene_info in genes.items():
        segment = gene_info["segment"]
        ref_segment_seq = segments[segment][ref_seq_id]
        ref_pos_to_aln_pos = map_ref_positions_to_alignment_positions(ref_segment_seq)

        for seq_id, seq_segment in segments[segment].items():
            aa_aligned_with_ref = 0
            aa_matching_ref = 0
            aa_different_to_ref = 0

            for codon_start in range(gene_info["start_pos"], gene_info["end_pos"] - 1, 3):
                codon_ref_positions = [codon_start, codon_start + 1, codon_start + 2]
                if any(ref_pos not in ref_pos_to_aln_pos for ref_pos in codon_ref_positions):
                    continue
                codon_aln_positions = [ref_pos_to_aln_pos[ref_pos] for ref_pos in codon_ref_positions]
                ref_codon = get_codon_at_alignment_positions(ref_segment_seq, codon_aln_positions)
                seq_codon = get_codon_at_alignment_positions(seq_segment, codon_aln_positions)
                ref_aa = translate_codon(ref_codon)
                seq_aa = translate_codon(seq_codon)
                if ref_aa is None or seq_aa is None:
                    continue

                aa_aligned_with_ref += 1
                if seq_aa == ref_aa:
                    aa_matching_ref += 1
                else:
                    aa_different_to_ref += 1
                    if seq_id != ref_seq_id:
                        for aln_pos in codon_aln_positions:
                            aa_difference_matrix[segment][seq_id][aln_pos] += 1

            gene_summary_rows.append(
                {
                    "seq_id": seq_id,
                    "gene": gene,
                    "segment": segment,
                    "aa_aligned_with_ref": aa_aligned_with_ref,
                    "aa_matching_ref": aa_matching_ref,
                    "aa_different_to_ref": aa_different_to_ref,
                }
            )

    return aa_difference_matrix, gene_summary_rows

def make_sequence_summary(segments: dict, ref_seq_id: str, gene_summary_rows: list) -> list:
    summary_rows = []
    valid_bases = ["A", "C", "G", "T"]

    for seq_id in segments[next(iter(segments.keys()))].keys():
        bases_aligned_with_ref = 0
        bases_matching_ref = 0
        for segment, segment_seqs in segments.items():
            ref_segment_seq = segment_seqs[ref_seq_id]
            seq_segment = segment_seqs[seq_id]
            for ref_base, seq_base in zip(ref_segment_seq, seq_segment):
                ref_base = ref_base.upper()
                seq_base = seq_base.upper()
                if ref_base in valid_bases and seq_base in valid_bases:
                    bases_aligned_with_ref += 1
                    if seq_base == ref_base:
                        bases_matching_ref += 1

        seq_gene_rows = [row for row in gene_summary_rows if row["seq_id"] == seq_id]
        aa_aligned_with_ref = sum(row["aa_aligned_with_ref"] for row in seq_gene_rows)
        aa_matching_ref = sum(row["aa_matching_ref"] for row in seq_gene_rows)
        aa_different_to_ref = sum(row["aa_different_to_ref"] for row in seq_gene_rows)
        summary_rows.append(
            {
                "seq_id": seq_id,
                "bases_aligned_with_ref": bases_aligned_with_ref,
                "bases_matching_ref": bases_matching_ref,
                "base_match_proportion": bases_matching_ref / bases_aligned_with_ref if bases_aligned_with_ref > 0 else None,
                "aa_aligned_with_ref": aa_aligned_with_ref,
                "aa_matching_ref": aa_matching_ref,
                "aa_different_to_ref": aa_different_to_ref,
                "aa_match_proportion": aa_matching_ref / aa_aligned_with_ref if aa_aligned_with_ref > 0 else None,
            }
        )

    return summary_rows

def write_LSV_summary_tables(summary_rows: list, gene_summary_rows: list, output_file: str) -> None:
    output_dir = os.path.dirname(output_file)
    output_prefix = os.path.splitext(output_file)[0]
    if output_dir != "":
        os.makedirs(output_dir, exist_ok=True)
    pl.DataFrame(summary_rows).write_csv(f"{output_prefix}_summary.tsv", separator="\t")
    pl.DataFrame(gene_summary_rows).write_csv(f"{output_prefix}_gene_summary.tsv", separator="\t")

def print_ref_white_regions(segments: dict, ref_seq_id: str, seqs_to_draw: list) -> None:
    for segment, segment_seqs in segments.items():
        ref_seq = segment_seqs[ref_seq_id]
        white_positions = [i for i, base in enumerate(ref_seq) if base.upper() not in ["A", "C", "G", "T"]]
        if len(white_positions) == 0:
            continue

        white_regions = []
        region_start = white_positions[0]
        previous_position = white_positions[0]
        for position in white_positions[1:]:
            if position != previous_position + 1:
                white_regions.append((region_start, previous_position + 1))
                region_start = position
            previous_position = position
        white_regions.append((region_start, previous_position + 1))

        print(f"\nReference white regions in {segment}: {len(white_regions)}")
        for region_idx, (start, end) in enumerate(white_regions, start=1):
            context_start = max(0, start - 10)
            context_end = min(len(ref_seq), end + 10)
            print(f"{segment} region {region_idx}: segment columns {start}-{end - 1}, length {end - start}")
            print(f"context columns {context_start}-{context_end - 1}")
            for seq_id in seqs_to_draw:
                print(f">{seq_id}")
                print(segment_seqs[seq_id][context_start:context_end])

def draw_LSV_diagram(alignment_file: str, seqs_to_draw_info: dict, ref_seq_id: str, deliniation_kmers: dict, genes: dict, output_file: str):
    seq_labels = {}
    for seq_id, seq_info in seqs_to_draw_info.items():
        short_seq_id = seq_id.split(" |")[0]  # SeqIO only uses bit before space
        seq_labels[short_seq_id] = seq_info.get("label", short_seq_id)
    seqs_to_draw = list(seq_labels.keys())
    ref_seq_id = ref_seq_id.split(" |")[0]  # SeqIO only uses bit before space
    sequences = parse_alignment_file(alignment_file, keep_sequences=seqs_to_draw)
    print(f"Parsed sequences: {list(sequences.keys())}")
    sequences = remove_gaps_where_all_gaps(sequences)
    segments = cut_to_segments(sequences, ref_seq_id, deliniation_kmers)
    aa_difference_matrix, gene_summary_rows = make_aa_difference_matrix_and_summary(segments, ref_seq_id, genes)
    summary_rows = make_sequence_summary(segments, ref_seq_id, gene_summary_rows)
    write_LSV_summary_tables(summary_rows, gene_summary_rows, output_file)
    segment_colours = {
        "segment_L": "#7B5EA7",
        "segment_M": "#E69F00",
        "segment_S": "#CC79A7",
    }
    coverage_colour = "#F6C7D9"
    aa_difference_colour = "#D55E00"
    overlapping_aa_difference_colour = "#7B1E3A"
    # figure set up - cols = segments, rows = gene annotations + sequences, widths = segment lengths
    segment_lengths = [len(next(iter(seg.values()))) for seg in segments.values()]
    # make sure ref sequence is the first row
    seqs_to_draw = [ref_seq_id] + [seq for seq in seqs_to_draw if seq != ref_seq_id]
    print_ref_white_regions(segments, ref_seq_id, seqs_to_draw)
    fig, axes = plt.subplots(
        nrows=len(seqs_to_draw) + 1,
        ncols=len(segments),
        figsize=(12, 0.38 * len(seqs_to_draw) + 0.8),
        gridspec_kw={"width_ratios": segment_lengths, "height_ratios": [0.45] + [0.45] + [0.35] * (len(seqs_to_draw) - 1)},
        squeeze=False,
    )

    for segment_idx, segment in enumerate(segments.keys()):
        ax = axes[0, segment_idx]
        draw_gene_annotations(ax, segment, segments[segment][ref_seq_id], genes)

    for row_idx, seq_id in enumerate(seqs_to_draw):
        for segment_idx, segment in enumerate(segments.keys()):
            ax = axes[row_idx + 1, segment_idx]
            seq_segment = segments[segment][seq_id]
            for col_idx, base in enumerate(seq_segment):
                base = base.upper()
                if base in ["A", "C", "G", "T"]:
                    if row_idx == 0:
                        base_colour = segment_colours[segment]
                    else:
                        aa_difference_count = aa_difference_matrix[segment][seq_id][col_idx]
                        if aa_difference_count > 1:
                            base_colour = overlapping_aa_difference_colour
                        elif aa_difference_count == 1:
                            base_colour = aa_difference_colour
                        else:
                            base_colour = coverage_colour
                else:
                    base_colour = "white"
                ax.add_patch(plt.Rectangle((col_idx, 0), 1, 1, color=base_colour))
            ax.set_xlim(0, len(seq_segment))
            ax.set_ylim(0, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            if row_idx == 0:
                ax.text(
                    len(seq_segment) / 2,
                    0.5,
                    segment.replace("_", " "),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white",
                )
            if segment_idx == 0:
                ax.set_ylabel(seq_labels[seq_id], rotation=0, ha="right", va="center", fontsize=7)
            for spine in ax.spines.values():
                spine.set_visible(False)

    # lets add a legend where pink = coverage, red = aa difference, dark red = overlapping aa difference
    legend_elements = [
        Patch(facecolor=coverage_colour, edgecolor=None, label="Covered bases"),
        Patch(facecolor=aa_difference_colour, edgecolor=None, label="AA difference"),
        # Patch(facecolor=overlapping_aa_difference_colour, edgecolor=None, label="Overlapping AA difference"),
    ]
    fig.legend(handles=legend_elements, loc="upper left", ncol=1, fontsize=8, frameon=False)
    plt.tight_layout()
    # plt.subplots_adjust(hspace=0.1, wspace=0.04)
    plt.savefig(output_file)
    plt.savefig(output_file.replace(".pdf", ".png"), dpi=300)


def draw_LSV_coverage_plot(coverage_df: pl.DataFrame, output_file: str) -> None:
    coverage_df = coverage_df.sort("position")
    plt.figure(figsize=(12, 4))
    plt.plot(coverage_df["position"], coverage_df["depth"], color="#1f77b4", linewidth=1)
    plt.fill_between(coverage_df["position"], coverage_df["depth"], color="#1f77b4", alpha=0.3)
    plt.yscale("log")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _: f"{int(value):,}" if value >= 1 else f"{value:g}"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.xlim(0, coverage_df["position"].max())
    plt.ylim(0, coverage_df["depth"].max() * 1.1)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()