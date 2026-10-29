from analysis.common.software import CLIQUESNV_JAR, TABULATE_ALIGNMENTS_EXECUTABLE
from analysis.common.paths import DATA_DIR, RESULTS_DIR
import polars as pl
import os
import json
import subprocess
from Bio import SeqIO, Phylo
import io
import seaborn as sns
from matplotlib import pyplot as plt
from statsmodels.stats.proportion import proportion_confint
from pathlib import Path

def parse_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id.replace(":", "_")] = str(record.seq)
    return sequences

def sort_seq_into_species(sequences, species_strings):
    species_sequences = {}
    for seq_id, seq in sequences.items():
        for species in species_strings:
            if species in seq_id:
                if species not in species_sequences:
                    species_sequences[species] = {}
                species_sequences[species][seq_id] = seq
                break
    return species_sequences

def mafft_align_sequences(fasta_file, output_file):
    mafft_command = ["mafft", "--auto", fasta_file]
    with open(output_file, "w") as out_f:
        subprocess.run(mafft_command, stdout=out_f, check=True)

def print_pc_identity_sharing_matrix_between_all_seqs(species_sequences):
    # Keep each sequence with its ID; accession prefixes do not identify species.
    all_sequences = [
        (seq_id, seq)
        for sequences in species_sequences.values()
        for seq_id, seq in sequences.items()
    ]
    all_seq_ids = [seq_id for seq_id, seq in all_sequences]
    n = len(all_seq_ids)
    identity_matrix = pl.DataFrame({"seq_id": all_seq_ids})
    for i in range(n):
        row = []
        seq_i = all_sequences[i][1]
        for j in range(n):
            seq_j = all_sequences[j][1]
            matches = sum(1 for a, b in zip(seq_i, seq_j) if a == b and a != '-')
            total_positions = sum(1 for a, b in zip(seq_i, seq_j) if a != '-' and b != '-')
            identity_percentage = (matches / total_positions * 100) if total_positions > 0 else 0
            row.append(identity_percentage)
        identity_matrix = identity_matrix.with_columns(pl.Series(f"seq_{i}", row))
    print(identity_matrix)

def determine_species_defining_snps(species_sequences):
    # make a dict with a position in the alignment, and the base for each species at that position. 
    # Only consider species defining if completely in that species and not in the other.
    species_defining_snps = {}
    species_list = list(species_sequences.keys())
    alignment_length = len(next(iter(species_sequences[species_list[0]].values())))
    for pos in range(alignment_length):
        # if either a gap, skip
        if any(seq[pos] == '-' for species in species_list for seq in species_sequences[species].values()):
            continue
        speciesA_bases = set(seq[pos] for seq in species_sequences[species_list[0]].values())
        speciesB_bases = set(seq[pos] for seq in species_sequences[species_list[1]].values())
        if len(speciesA_bases) == 1 and len(speciesB_bases) == 1 and speciesA_bases != speciesB_bases:
            species_defining_snps[pos] = {
                species_list[0]: next(iter(speciesA_bases)),
                species_list[1]: next(iter(speciesB_bases))
            }
    return species_defining_snps

def make_single_reference(sequences, id_for_new_ref, output_file):
    ref_seq = sequences[id_for_new_ref]
    with open(output_file, "w") as f:
        f.write(f">{id_for_new_ref}\n")
        f.write(f"{ref_seq}\n")
    # faidx the new reference
    subprocess.run(["samtools", "faidx", output_file], check=True)

def align_reads_to_reference_with_minimap2(fastq_r1, fastq_r2, reference_fasta, output_bam, threads=8):
    if threads < 1:
        raise ValueError("threads must be at least 1")
    # Samtools counts additional threads beyond its main thread.
    samtools_threads = str(threads - 1)
    sam_file = output_bam.replace(".bam", ".sam")
    minimap2_command = [
        "minimap2",
        "-t", str(threads),
        "-ax", "sr",
        reference_fasta,
        fastq_r1,
        fastq_r2
    ]
    with open(sam_file, "w") as sam_f:
        subprocess.run(minimap2_command, stdout=sam_f, check=True)
    # convert to bam and sort
    subprocess.run(["samtools", "view", "-@", samtools_threads, "-bS", sam_file, "-o", output_bam], check=True)
    subprocess.run(["samtools", "sort", "-@", samtools_threads, "-o", output_bam, output_bam], check=True)
    subprocess.run(["samtools", "index", "-@", samtools_threads, output_bam], check=True)
    # delete the intermediate sam file
    os.remove(sam_file)

def run_tabulate_alignments(sample_bam: str, position: int, tabulate_alignments_executable: str, target_id: str, mapq_threshold: int) -> pl.DataFrame:
    # SNP dictionary positions are zero-based; tabulate-alignments uses one-based positions.
    position += 1
    # use tabulate alignments to get pileup across kmer
    command = f"{tabulate_alignments_executable} -reads {sample_bam} -range {target_id}:{position}-{position} -mq-bins 0 {mapq_threshold}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    header = ["file","chromosome","position","mq>=0:A","mq>=0:C","mq>=0:G","mq>=0:T","mq>=0:-","mq>=0:+",f"mq>={mapq_threshold}:A",f"mq>={mapq_threshold}:C",f"mq>={mapq_threshold}:G",f"mq>={mapq_threshold}:T",f"mq>={mapq_threshold}:-",f"mq>={mapq_threshold}:+"] 
    df = pl.read_csv(io.StringIO(result.stdout), separator='\t', has_header=True)
    return df

def compute_snp_counts_at_position_in_bam(sample_bam: str, position: int, base_A: str, base_B: str, target_id: str, tabulate_alignments_executable: str=TABULATE_ALIGNMENTS_EXECUTABLE, mapq_threshold: int=10) -> float:
    df = run_tabulate_alignments(sample_bam, position, tabulate_alignments_executable, target_id, mapq_threshold)
    # capitalise bases
    base_A = base_A.upper()
    base_B = base_B.upper()
    # sum the counts for bases_A and bases_B
    count_A = df[f"mq>={mapq_threshold}:{base_A}"].sum()
    count_B = df[f"mq>={mapq_threshold}:{base_B}"].sum()
    total_count = 0
    for base in ['A', 'C', 'G', 'T']:
        total_count += df[f"mq>={mapq_threshold}:{base}"].sum()
    if total_count == 0:
        return None
    non_AB_count = total_count - (count_A + count_B)
    return count_A, count_B, non_AB_count

def make_table_of_snp_counts_for_samples(sample_bams: list[str], species_defining_snps: dict, target_id: str, tabulate_alignments_executable: str=TABULATE_ALIGNMENTS_EXECUTABLE, mapq_threshold: int=10) -> pl.DataFrame:
    rows = []
    species_names = list(next(iter(species_defining_snps.values())).keys())
    species_A, species_B = species_names
    for sample_bam in sample_bams:
        sample_name = os.path.basename(sample_bam).replace(".bam", "")
        for pos, bases in species_defining_snps.items():
            base_A = bases[list(bases.keys())[0]]
            base_B = bases[list(bases.keys())[1]]
            counts = compute_snp_counts_at_position_in_bam(sample_bam, pos, base_A, base_B, target_id, tabulate_alignments_executable, mapq_threshold)
            if counts is not None:
                count_A, count_B, non_AB_count = counts
                rows.append({
                    "sample": sample_name,
                    "position": pos,
                    f"base_{species_A}": base_A,
                    f"base_{species_B}": base_B,
                    f"count_{species_A}": count_A,
                    f"count_{species_B}": count_B,
                    f"other_bases_count": non_AB_count
                })
    return pl.DataFrame(rows)

def plot_forest_plot_of_snp_counts(snp_counts_df: pl.DataFrame, output_file: str):
    # A panel per mixture, matching L1/CulexAus1 and L2/CulexAus2 in the heatmap.
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 7), sharex=True)
    # Use all counted bases as the denominator for each allele fraction.
    total = (
        pl.col("count_pipiens")
        + pl.col("count_australicus")
        + pl.col("other_bases_count")
    )
    snp_counts_df = snp_counts_df.with_columns(
        (pl.col("count_pipiens") / total).alias("pipiens_ratio"),
        (pl.col("count_australicus") / total).alias("australicus_ratio"),
        (pl.col("other_bases_count") / total).alias("other_ratio"),
    )
    # Calculate 95% Wilson intervals on arrays, then add named columns.
    totals = snp_counts_df.select(total.alias("total")).to_series().to_numpy()
    for allele, count_column in [
        ("pipiens", "count_pipiens"),
        ("australicus", "count_australicus"),
        ("other", "other_bases_count"),
    ]:
        lower, upper = proportion_confint(
            snp_counts_df[count_column].to_numpy(),
            totals,
            alpha=0.05,
            method="wilson",
        )
        snp_counts_df = snp_counts_df.with_columns(
            pl.Series(f"{allele}_ci_lower", lower),
            pl.Series(f"{allele}_ci_upper", upper),
        )
    # plot scatter of pipiens ratio with error bars and other ratio with error bars, offset from each other horizontally
    positions = snp_counts_df["position"].unique().sort().to_list()
    position_to_tick = {position: index for index, position in enumerate(positions)}
    mixture_titles = {
        "mixture1_aligned": "Mixture 1",
        "mixture2_aligned": "Mixture 2",
    }
    reference_fractions = {"mixture1_aligned": 0.2, "mixture2_aligned": 0.1}
    for i, (mixture, title) in enumerate(mixture_titles.items()):
        ax = axes[i]
        ax.hlines(reference_fractions[mixture], xmin=0, xmax=len(positions) - 1,
                  color="0.35", linestyle="--", linewidth=1, zorder=0)
        mixture_df = snp_counts_df.filter(pl.col("sample") == mixture).sort("position")
        x = pl.Series([position_to_tick[position] for position in mixture_df["position"]])
        ax.errorbar(x - 0.12, mixture_df["pipiens_ratio"], yerr=[mixture_df["pipiens_ratio"] - mixture_df["pipiens_ci_lower"], mixture_df["pipiens_ci_upper"] - mixture_df["pipiens_ratio"]], fmt='o', markersize=4, label='Pipiens-associated allele', color='#0072B2')
        ax.errorbar(x + 0.12, mixture_df["other_ratio"], yerr=[mixture_df["other_ratio"] - mixture_df["other_ci_lower"], mixture_df["other_ci_upper"] - mixture_df["other_ratio"]], fmt='s', markersize=4, label='Non-reference alleles', color='#D55E00')
        ax.set_title(title)
        ax.set_axisbelow(True)
        ax.xaxis.grid(True, color='0.9', linewidth=0.6)
        ax.set_xticks(range(len(positions)), labels=[str(position) for position in positions])
        ax.tick_params(axis="x", labelbottom=True, labelrotation=90, labelsize=10)
        ax.set_ylabel('Fraction of reads in pileup', fontsize=10)
        sns.despine(ax=ax, top=True, right=True, trim=True)
    axes[0].legend(loc="upper left", frameon=False)
    axes[-1].set_xlabel("Species-defining polymorphism")
    # save the figure
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

def run_cliquesnv(
    bam_file: str,
    output_dir: str | Path,
    cliquesnv_jar: str,
    mode: str,
    threads: int,
) -> tuple[Path, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_stem = Path(bam_file).stem
    fasta_output = output_dir / f"{output_stem}.fasta"
    json_output = output_dir / f"{output_stem}.json"
    if all(path.is_file() and path.stat().st_size > 0 for path in (fasta_output, json_output)):
        try:
            report = json.loads(json_output.read_text())
        except json.JSONDecodeError:
            report = None
        if isinstance(report, dict) and report.get("error") in (None, "", "none"):
            return fasta_output, json_output

    # CliqueSNV's Illumina mode expects SAM input, including the header.
    sam_file = output_dir / f"{output_stem}.sam"
    subprocess.run(
        ["samtools", "view", "-@", str(max(threads - 1, 0)), "-h",
         "-o", str(sam_file), bam_file],
        check=True,
    )
    command = [
        "java",
        "-jar",
        cliquesnv_jar,
        "-m",
        mode,
        "-in",
        str(sam_file),
        "-outDir",
        str(output_dir),
        "-threads",
        str(threads),
    ]
    subprocess.run(command, check=True)
    if not all(path.is_file() and path.stat().st_size > 0 for path in (fasta_output, json_output)):
        raise RuntimeError(f"CliqueSNV did not produce complete outputs in {output_dir}")
    report = json.loads(json_output.read_text())
    if report.get("error") not in (None, "", "none"):
        raise RuntimeError(f"CliqueSNV reported: {report['error']}")
    return fasta_output, json_output

def combine_fasta_files(fasta_files: list[str], output_file: str) -> None:
    records = {}
    for fasta_file in fasta_files:
        mixture = Path(fasta_file).parent.name
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id.replace(":", "_")
            if mixture in ("mixture1", "mixture2"):
                name = f"{mixture}_{name}"
            if name in records:
                if records[name].seq == record.seq:
                    continue
                raise ValueError(f"Different sequences share FASTA ID: {name}")
            record.id = name
            record.name = name
            record.description = ""
            records[name] = record
    SeqIO.write(records.values(), output_file, "fasta")


def run_iqtree_on_fasta(fasta_file: str, output_prefix: str, threads: int=4) -> None:
    subprocess.run(
        ["iqtree", "-s", fasta_file, "--prefix", output_prefix,
         "-m", "GTR+F", "-B", "1000", "-T", str(threads)],
        check=True,
    )


def plot_haplotype_tree(tree_file: str, output_file: str) -> None:
    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    colours = {
        tip.name: ("#0072B2" if tip.name.startswith("mixture1_") else
                   "#D55E00" if tip.name.startswith("mixture2_") else "black")
        for tip in tree.get_terminals()
    }
    fig, ax = plt.subplots(figsize=(14, max(4, len(colours) * 0.35)))
    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_func=lambda clade: clade.name if clade.is_terminal() else None,
        label_colors=colours, show_confidence=True,
    )
    for text in ax.texts:
        text.set_fontsize(8)
    ax.set_xlabel("Substitutions/site")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.grid(False)
    sns.despine(ax=ax, left=True, trim=True)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    (RESULTS_DIR / "pip_australicus_mixtures_check").mkdir(parents=True, exist_ok=True)
    fasta_file = f"{DATA_DIR}/references/Culex_pipiens_australicus_COX1_six_references.fasta"
    species_strings = ["pipiens", "australicus"]

    aligned_output_file = f"{RESULTS_DIR}/pip_australicus_mixtures_check/pipiens_australicus_res_aligned.fasta"
    mafft_align_sequences(fasta_file, aligned_output_file)

    sequences = parse_fasta(aligned_output_file)
    species_sequences = sort_seq_into_species(sequences, species_strings)

    print_pc_identity_sharing_matrix_between_all_seqs(species_sequences)

    species_defining_snps = determine_species_defining_snps(species_sequences)

    target_id = "MN389460.1_1446-2982_Culex_pipiens_molestus_isolate_W82_COX1"
    make_single_reference(sequences, target_id, f"{RESULTS_DIR}/pip_australicus_mixtures_check/pipiens_ref.fasta")
    seq_data = {
        "mixture1": {
            "R1": f"{DATA_DIR}/reads/trimmed/24-Cosmos-CulexAus1_1_trimmed.fastq.gz",
            "R2": f"{DATA_DIR}/reads/trimmed/24-Cosmos-CulexAus1_2_trimmed.fastq.gz",
        },
        "mixture2": {
            "R1": f"{DATA_DIR}/reads/trimmed/24-Cosmos-CulexAus2_1_trimmed.fastq.gz",
            "R2": f"{DATA_DIR}/reads/trimmed/24-Cosmos-CulexAus2_2_trimmed.fastq.gz",
        }
    }
    for mixture_name, read_files in seq_data.items():
        # Use new BAMs so cached alignments with colon-containing reference names are not reused.
        output_bam = f"{RESULTS_DIR}/pip_australicus_mixtures_check/{mixture_name}_aligned.bam"
        if not os.path.exists(output_bam):
            align_reads_to_reference_with_minimap2(read_files["R1"], read_files["R2"], f"{RESULTS_DIR}/pip_australicus_mixtures_check/pipiens_ref.fasta", output_bam, threads=16)

    output_file = f"{RESULTS_DIR}/abundance-validation/tables/Cx_australicus_pipiens_snp_counts_table.tsv"
    sample_bams = [f"{RESULTS_DIR}/pip_australicus_mixtures_check/{mixture}_aligned.bam" for mixture in seq_data.keys()]
    if not os.path.exists(output_file):
        snp_table = make_table_of_snp_counts_for_samples(sample_bams, species_defining_snps, target_id)
        snp_table.write_csv(output_file, separator='\t')

    results_df = pl.read_csv(output_file, separator='\t')
    plot_output_file = f"{RESULTS_DIR}/abundance-validation/plots/Cx_australicus_pipiens_mixtures_snp_counts_forest_plot.png"
    plot_forest_plot_of_snp_counts(results_df, plot_output_file)

    # run CliqueSNV on the BAMs to see if it can detect the two species in the mixtures
    cliquesnv_jar = CLIQUESNV_JAR
    for mixture_name in seq_data.keys():
        output_dir = Path(f"{RESULTS_DIR}/pip_australicus_mixtures_check/{mixture_name}")
        bam_file = f"{RESULTS_DIR}/pip_australicus_mixtures_check/{mixture_name}_aligned.bam"
        run_cliquesnv(bam_file, output_dir, cliquesnv_jar, mode="snv-illumina", threads=16)

    # combine clique fastas with the original reference sequences for alignment and SNP analysis
    clique_fasta_files = [fasta_file, f"{RESULTS_DIR}/pip_australicus_mixtures_check/mixture1/mixture1_aligned.fasta", f"{RESULTS_DIR}/pip_australicus_mixtures_check/mixture2/mixture2_aligned.fasta"]
    combined_fasta_output = f"{RESULTS_DIR}/pip_australicus_mixtures_check/combined_clique_and_reference_labelled.fasta"
    combine_fasta_files(clique_fasta_files, combined_fasta_output)
    aligned_output = f"{RESULTS_DIR}/pip_australicus_mixtures_check/combined_clique_and_reference_labelled_aligned.fasta"
    if not os.path.exists(aligned_output):
        mafft_align_sequences(combined_fasta_output, aligned_output)
    iqtree_output_prefix = f"{RESULTS_DIR}/pip_australicus_mixtures_check/combined_clique_and_reference_labelled_aligned"
    if not os.path.exists(f"{iqtree_output_prefix}.treefile"):
        run_iqtree_on_fasta(aligned_output, iqtree_output_prefix, threads=14)
    plot_haplotype_tree(
        f"{iqtree_output_prefix}.treefile",
        f"{RESULTS_DIR}/abundance-validation/plots/Cx_australicus_pipiens_haplotype_tree.png",
    )
        
    


if __name__ == "__main__":
    main()
