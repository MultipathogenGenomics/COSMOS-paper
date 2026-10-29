# Prepare mosquito sequence alignments and percentage-difference clusters.

from pathlib import Path
import re
import subprocess
import os
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance


UNAMBIGUOUS_BASES = set("ACGT")


def read_fasta(fasta_file: str, source: str) -> list[dict]:
    records = []
    header = None
    sequence = []
    with open(fasta_file) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(
                        {
                            "header": header,
                            "sequence": "".join(sequence).upper(),
                            "source": source,
                        }
                    )
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence found before a header in {fasta_file}")
                sequence.append(line)
    if header is not None:
        records.append(
            {
                "header": header,
                "sequence": "".join(sequence).upper(),
                "source": source,
            }
        )
    return records


def remove_duplicate_headers(
    records: list[dict],
) -> tuple[list[dict], pd.DataFrame]:
    retained = []
    duplicate_rows = []
    seen_headers = {}
    for record in records:
        if record["header"] in seen_headers:
            duplicate_rows.append(
                {
                    "header": record["header"],
                    "retained_source": seen_headers[record["header"]],
                    "excluded_source": record["source"],
                }
            )
        else:
            retained.append(record)
            seen_headers[record["header"]] = record["source"]
    return retained, pd.DataFrame(
        duplicate_rows,
        columns=["header", "retained_source", "excluded_source"],
    )


def filter_records_by_species(records: list[dict], species_name: str) -> list[dict]:
    species_text = f"_{species_name.replace(' ', '_')}_"
    species_records = [record for record in records if species_text in record["header"]]
    if not species_records:
        raise ValueError(f"No {species_name} sequences were found")
    return species_records


def collapse_identical_sequences(records: list[dict]) -> tuple[list[dict], pd.DataFrame]:
    representative_by_sequence = {}
    representative_records = []
    membership_rows = []
    for record in records:
        if record["sequence"] not in representative_by_sequence:
            representative_by_sequence[record["sequence"]] = record["header"]
            representative_records.append(record)
        membership_rows.append(
            {
                "header": record["header"],
                "representative_header": representative_by_sequence[record["sequence"]],
            }
        )
    membership = pd.DataFrame(membership_rows)
    leaf_sizes = membership["representative_header"].value_counts()
    membership["n_sequences_in_leaf"] = membership["representative_header"].map(leaf_sizes)
    return representative_records, membership


def shorten_sequence_header(header: str) -> tuple[str, str]:
    accession = header.split("_", 1)[0]
    if "isolate_" not in header:
        raise ValueError(f"Could not find isolate name in {header}")
    isolate = header.split("isolate_", 1)[1]
    isolate = re.split("_cytochrome|_haplogroup|_mitochondrion", isolate, maxsplit=1)[0]
    return accession, isolate


def load_cx_annulirostris_types(metadata_file: str) -> dict:
    metadata = pd.read_csv(metadata_file, sep="\t").dropna(subset=["H", "mtDNA (COI)"])
    type_by_isolate = {}
    for _, row in metadata.iterrows():
        coi_type = row["mtDNA (COI)"]
        match = re.search(r"ann (AUS|S-AUS|PNG-1|PNG-2|SI)", coi_type)
        if match is not None:
            type_by_isolate[row["H"]] = match.group(1)
    return type_by_isolate


def add_leaf_labels(
    sequence_metadata: pd.DataFrame,
    cluster_ids: np.ndarray,
    linkage_matrix: np.ndarray,
    leaf_sizes: dict,
    config: dict,
) -> pd.DataFrame:
    sequence_metadata = sequence_metadata.copy()
    shortened_headers = sequence_metadata["header"].apply(shorten_sequence_header)
    sequence_metadata["accession"] = shortened_headers.str[0]
    sequence_metadata["isolate"] = shortened_headers.str[1]

    if config["leaf_annotation"] == "cx_type":
        type_by_isolate = load_cx_annulirostris_types(config["metadata_file"])
        sequence_metadata["group"] = sequence_metadata["isolate"].map(type_by_isolate)
        if sequence_metadata["group"].isna().any():
            missing = sequence_metadata.loc[sequence_metadata["group"].isna(), "isolate"].tolist()
            raise ValueError(f"No Culex annulirostris type found for: {missing}")
    else:
        unique_clusters = sorted(set(cluster_ids))
        if len(unique_clusters) != 3:
            raise ValueError(
                f"Expected three Aedes notoscriptus clusters at "
                f"{config['pid_threshold']} identity, found {len(unique_clusters)}"
            )
        sequence_metadata["group"] = [f"population {cluster_id}" for cluster_id in cluster_ids]

    sequence_metadata["leaf_label"] = (
        sequence_metadata["accession"]
        + " | "
        + sequence_metadata["isolate"]
    )
    sequence_metadata["n_sequences_in_leaf"] = sequence_metadata["header"].map(leaf_sizes)
    sequence_metadata.loc[
        sequence_metadata["n_sequences_in_leaf"] > 1,
        "leaf_label",
    ] += "*"
    return sequence_metadata


def make_supplementary_tree_table(
    records: list[dict],
    membership: pd.DataFrame,
    sequence_metadata: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    rows = []
    type_by_isolate = None
    if config["leaf_annotation"] == "cx_type":
        type_by_isolate = load_cx_annulirostris_types(config["metadata_file"])
    group_by_header = sequence_metadata.set_index("header")["group"].to_dict()
    leaf_label_by_header = sequence_metadata.set_index("header")["leaf_label"].to_dict()
    representative_by_header = membership.set_index("header")["representative_header"].to_dict()

    for record in records:
        accession, isolate = shorten_sequence_header(record["header"])
        representative_header = representative_by_header[record["header"]]
        subtype = (
            type_by_isolate[isolate]
            if type_by_isolate is not None
            else group_by_header[representative_header]
        )
        rows.append(
            {
                "species": config["species_name"],
                "ncbi_accession": accession,
                "sample_id": isolate,
                "subtype": subtype,
                "leaf_label": leaf_label_by_header[representative_header],
            }
        )
    return pd.DataFrame(rows)


def write_fasta(records: list[dict], output_file: str) -> pd.DataFrame:
    metadata_rows = []
    with open(output_file, "w") as handle:
        for index, record in enumerate(records, start=1):
            sequence_id = f"sequence_{index:04d}"
            handle.write(f">{sequence_id}\n{record['sequence']}\n")
            metadata_rows.append(
                {
                    "sequence_id": sequence_id,
                    "header": record["header"],
                    "source": record["source"],
                    "sequence_length": len(record["sequence"]),
                }
            )
    return pd.DataFrame(metadata_rows)


def run_mafft(input_fasta: str, output_fasta: str, executable: str = "mafft") -> None:
    with open(output_fasta, "w") as output_handle:
        subprocess.run(
            [executable, "--auto", input_fasta],
            stdout=output_handle,
            check=True,
        )


def compute_pairwise_differences(
    aligned_records: list[dict],
    min_shared_sites: int,
) -> tuple[list[str], np.ndarray, pd.DataFrame]:
    sequence_ids = [record["header"] for record in aligned_records]
    sequences = [record["sequence"] for record in aligned_records]
    if len(set(len(sequence) for sequence in sequences)) != 1:
        raise ValueError("MAFFT output sequences do not all have the same aligned length")
    difference_matrix = np.zeros((len(sequences), len(sequences)), dtype=float)
    pairwise_rows = []

    for first_index, first_sequence in enumerate(sequences):
        for second_index in range(first_index + 1, len(sequences)):
            second_sequence = sequences[second_index]
            shared_sites = [
                (first_base, second_base)
                for first_base, second_base in zip(first_sequence, second_sequence)
                if first_base in UNAMBIGUOUS_BASES and second_base in UNAMBIGUOUS_BASES
            ]
            if len(shared_sites) < min_shared_sites:
                raise ValueError(
                    f"{sequence_ids[first_index]} and {sequence_ids[second_index]} "
                    f"share only {len(shared_sites)} unambiguous sites; "
                    f"at least {min_shared_sites} are required"
                )
            n_differences = sum(first_base != second_base for first_base, second_base in shared_sites)
            proportion_difference = n_differences / len(shared_sites)
            difference_matrix[first_index, second_index] = proportion_difference
            difference_matrix[second_index, first_index] = proportion_difference
            pairwise_rows.append(
                {
                    "sequence_id_1": sequence_ids[first_index],
                    "sequence_id_2": sequence_ids[second_index],
                    "shared_sites": len(shared_sites),
                    "n_differences": n_differences,
                    "proportion_difference": proportion_difference,
                    "percent_difference": proportion_difference * 100,
                }
            )
    return sequence_ids, difference_matrix, pd.DataFrame(pairwise_rows)


def prepare_sequence_tree(config: dict, output_dir: str) -> dict:
    species_dir = Path(output_dir) / config["output_prefix"]
    species_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for fasta_file in config["sequence_fastas"]:
        records.extend(read_fasta(fasta_file, Path(fasta_file).name))
    records = filter_records_by_species(records, config["species_name"])
    records, duplicate_headers = remove_duplicate_headers(records)
    all_records = records
    records, sequence_membership = collapse_identical_sequences(records)

    combined_fasta = species_dir / "combined.fasta"
    aligned_fasta = species_dir / "combined.aligned.fasta"
    sequence_metadata = write_fasta(records, str(combined_fasta))
    duplicate_headers.to_csv(species_dir / "duplicate_headers_excluded.tsv", sep="\t", index=False)

    run_mafft(str(combined_fasta), str(aligned_fasta), config.get("mafft_executable", "mafft"))
    aligned_records = read_fasta(str(aligned_fasta), "aligned")
    sequence_ids, difference_matrix, pairwise_df = compute_pairwise_differences(
        aligned_records,
        config.get("min_shared_sites", 100),
    )
    pairwise_df.to_csv(species_dir / "pairwise_differences.tsv", sep="\t", index=False)

    linkage_matrix = hierarchy.linkage(
        distance.squareform(difference_matrix),
        method="complete",
    )
    pid_threshold = config.get("pid_threshold")
    if pid_threshold is None:
        if config["leaf_annotation"] == "ae_population":
            raise ValueError("pid_threshold is required to assign Aedes populations")
        cluster_ids = np.full(len(sequence_ids), pd.NA)
    else:
        cluster_ids = hierarchy.fcluster(
            linkage_matrix,
            t=1 - pid_threshold,
            criterion="distance",
        )
    sequence_metadata = add_leaf_labels(
        sequence_metadata,
        cluster_ids,
        linkage_matrix,
        sequence_membership.set_index("representative_header")["n_sequences_in_leaf"].to_dict(),
        config,
    )
    sequence_metadata.to_csv(species_dir / "sequence_metadata.tsv", sep="\t", index=False)
    cluster_df = sequence_metadata.merge(
        pd.DataFrame({"sequence_id": sequence_ids, "cluster_id": cluster_ids}),
        on="sequence_id",
        how="inner",
        validate="one_to_one",
    )
    cluster_df.to_csv(species_dir / "cluster_assignments.tsv", sep="\t", index=False)
    supplementary_table = make_supplementary_tree_table(
        all_records,
        sequence_membership,
        sequence_metadata,
        config,
    )

    return {
        **config,
        "sequence_ids": sequence_ids,
        "sequence_metadata": sequence_metadata,
        "difference_matrix": difference_matrix,
        "linkage_matrix": linkage_matrix,
        "cluster_assignments": cluster_df,
        "supplementary_table": supplementary_table,
    }

def run_iqtree_on_fasta(fasta_file: str, output_prefix: str, threads: int=4) -> None:
    command = f"iqtree -s {fasta_file} --prefix {output_prefix} -m GTR+F -B 1000 -T {threads}" # +R6 to check model seleciton?
    subprocess.run(command, shell=True)

def make_phylogenetic_trees_for_species(config: dict, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    for species in config:
        name = species["species_name"]
        output_prefix = species["output_prefix"]
        fastas = species["sequence_fastas"]
        # combine fastas into one file
        combined_fasta = os.path.join(output_dir, output_prefix, f"{output_prefix}_combined.fasta")
        with open(combined_fasta, "w") as outfile:
            for fasta in fastas:
                with open(fasta, "r") as infile:
                    outfile.write(infile.read())
        # run mafft to align sequences
        aligned_fasta = os.path.join(output_dir, output_prefix, f"{output_prefix}_aligned.fasta")
        run_mafft(combined_fasta, aligned_fasta)
        # run iqtree to build phylogenetic tree
        run_iqtree_on_fasta(aligned_fasta, os.path.join(output_dir, output_prefix, f"{output_prefix}_phylogenetic_tree"), threads=4)
