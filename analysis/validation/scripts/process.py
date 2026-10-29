# Prepare data for abundance validation plots.

import re

import numpy as np
import polars as pl

from analysis.common.abundance import compute_mixture_abundance_error_records, summarise_mixture_abundance_errors
from analysis.common.mosquito_ids import return_species_from_target_id


GENUS_ORDER = {
    "Mansonia": 0,
    "Anopheles": 1,
    "Coquillettidia": 2,
    "Aedes": 3,
    "Culex": 4,
}

SPECIES_ORDER_OVERRIDES = {
    "aedes vigilax": 10_000,
    "culex annulirostris": -1,
    "culex australicus": 10_000,
    "coquillettidia linealis": 10_000,
}


def load_control_data(abundance_file: str, annotations_file: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    abundance_df = pl.read_csv(abundance_file, separator="\t")
    annotations_df = pl.read_csv(annotations_file, separator="\t").rename({"sampleid": "sample"})
    return abundance_df, annotations_df


def load_control_annotations(annotations_file: str) -> pl.DataFrame:
    return pl.read_csv(annotations_file, separator="\t").rename({"sampleid": "sample"})


def get_group_data_and_labels(abundance_df: pl.DataFrame, ordered_samples: list[str], sample_to_label: dict[str, str]) -> tuple:
    ordered_sample_df = pl.DataFrame({"sample": ordered_samples, "sample_order": np.arange(len(ordered_samples))})
    group_df = ordered_sample_df.join(abundance_df, on="sample", how="left").sort("sample_order")
    labels = [sample_to_label.get(sample, sample) for sample in group_df["sample"].to_list()]
    data = group_df.drop(["sample", "sample_order"]).to_pandas()
    return data, labels


def get_genus(species_name: str) -> str:
    genus = species_name.strip().split()[0]
    if genus == "Ochlerotatus":
        return "Aedes"
    return genus


def normalise_species_name(species_name: str) -> str:
    normalised_name = species_name.strip().replace("_", " ")
    if normalised_name.startswith("Ochlerotatus "):
        normalised_name = normalised_name.replace("Ochlerotatus ", "Aedes ", 1)
    return normalised_name


def letter_to_number(letter: str) -> int:
    l = letter.strip().upper()
    if l == "L1":
        return 12
    if l == "L2":
        return 13
    if len(l) != 1 or not ("A" <= l <= "Z"):
        raise ValueError("Expected A-Z, L1, or L2")
    return ord(l) - ord("A") + 1


def convert_species_letter_to_number(species_name: str, species_to_order: dict[str, str]) -> int:
    letter = species_to_order[normalise_species_name(species_name).lower()]
    if not letter:
        return 10_000
    return letter_to_number(letter)


def get_species_plot_label(species_name: str, species_to_order: dict[str, str]) -> str:
    return species_to_order[normalise_species_name(species_name).lower()]


def get_species_sort_order(species_name: str, species_to_order: dict[str, str]) -> int:
    species_key = normalise_species_name(species_name).lower()
    return SPECIES_ORDER_OVERRIDES.get(
        species_key,
        convert_species_letter_to_number(species_name, species_to_order),
    )


def format_species_sample_label(species_name: str, species_label: str) -> str:
    italic_name = normalise_species_name(species_name).replace(" ", r"\ ")
    return f"$\\it{{{italic_name}}}$"


def fix_COSMOS6_label(label: str) -> str:
    cleaned = str(label).replace("\n", " ")
    matches = re.findall(r"([0-9]*\.?[0-9]+)\s*%\s*([A-Za-z0-9]+)", cleaned)
    label_dict = {name.upper(): float(value) for value, name in matches}
    culex_label = "I" if "I" in label_dict else "L1"
    if "B" not in label_dict or culex_label not in label_dict:
        return label

    l1_adjusted = label_dict[culex_label] * 2.5
    total = label_dict["B"] + l1_adjusted
    b_frac = label_dict["B"] / total
    l1_frac = l1_adjusted / total
    new_label = f"{b_frac*100:.1f}% B, {l1_frac*100:.1f}% {culex_label}"
    return new_label


def format_ratio_label(label: str) -> str:
    percentages = re.findall(r"([0-9]*\.?[0-9]+)\s*%", str(label))
    if len(percentages) >= 2:
        return f"{percentages[0]}:{percentages[1]}"
    return str(label)


def format_mixture_label(samples: list[str], proportions: list[float], species_to_order: dict[str, str]) -> str:
    label_parts = []
    for sample, proportion in zip(samples, proportions):
        species_label = species_to_order[normalise_species_name(sample).lower()]
        label_parts.append((letter_to_number(species_label), f"{proportion * 100:g}% {species_label}"))
    label_parts = sorted(label_parts, key=lambda label_part: label_part[0])
    return ", ".join(label for _, label in label_parts)


def select_representative_mixture_ids(mixture_df: pl.DataFrame, n_mixtures: int) -> list[int]:
    sorted_mixtures = mixture_df.sort("kma_rmse").to_dicts()
    if len(sorted_mixtures) <= n_mixtures:
        return [row["mixture_id"] for row in sorted_mixtures]

    selected_ids = []
    candidate_indices = np.linspace(0, len(sorted_mixtures) - 1, n_mixtures)
    for index in candidate_indices:
        mixture_id = sorted_mixtures[int(round(index))]["mixture_id"]
        if mixture_id not in selected_ids:
            selected_ids.append(mixture_id)

    selected_species = set()
    for row in sorted_mixtures:
        if row["mixture_id"] in selected_ids:
            selected_species.update([row[f"sample_{index}"] for index in range(1, 5)])

    all_species = {
        row[f"sample_{index}"]
        for row in sorted_mixtures
        for index in range(1, 5)
    }
    for species in sorted(all_species - selected_species):
        replacement = next(row for row in sorted_mixtures if species in [row[f"sample_{index}"] for index in range(1, 5)])
        if replacement["mixture_id"] not in selected_ids:
            selected_ids[-1] = replacement["mixture_id"]
            selected_species.update([replacement[f"sample_{index}"] for index in range(1, 5)])

    return selected_ids[:n_mixtures]


def get_abundance_value(row: dict, abundance_column_lookup: dict[str, str], species_name: str) -> float:
    abundance_column = abundance_column_lookup.get(normalise_species_name(species_name).lower())
    if abundance_column is None:
        return 0.0
    value = row.get(abundance_column)
    return 0.0 if value is None else float(value)


def prepare_in_silico_mixture_df(
    abundance_df: pl.DataFrame,
    metadata_file: str,
) -> pl.DataFrame:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    error_df = compute_mixture_abundance_error_records(
        abundance_df,
        metadata_df,
        tool="abundance_within_filter_pass_reads",
        relative_error_margin=0.10,
    )
    mixture_rmse_df = (
        error_df
        .group_by("mixture_id")
        .agg(pl.col("mixture_rmse").first().alias("kma_rmse"))
    )
    mixture_df = (
        abundance_df
        .with_columns(pl.col("sample").str.extract(r"mixture_([0-9]+)", 1).cast(pl.Int64).alias("mixture_id"))
        .join(metadata_df, on="mixture_id", how="inner")
        .join(mixture_rmse_df, on="mixture_id", how="left")
    )
    return mixture_df


def prepare_in_silico_mixture_heatmap_group(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    abundance_columns: list[str],
    species_to_order: dict[str, str],
    n_mixtures: int = 18,
) -> tuple:
    mixture_df = prepare_in_silico_mixture_df(abundance_df, metadata_file)
    abundance_column_lookup = {
        normalise_species_name(column).lower(): column
        for column in abundance_df.columns
        if column != "sample"
    }
    selected_ids = select_representative_mixture_ids(mixture_df, n_mixtures)
    selected_mixtures = mixture_df.filter(pl.col("mixture_id").is_in(selected_ids))
    mixture_order_df = pl.DataFrame({"mixture_id": selected_ids, "mixture_order": list(range(len(selected_ids)))})
    selected_mixtures = selected_mixtures.join(mixture_order_df, on="mixture_id", how="left").sort("mixture_order")

    rows = []
    labels = []
    for row in selected_mixtures.to_dicts():
        true_proportions = [float(row[f"proportion_{index}"]) for index in range(1, 5)]
        samples = [row[f"sample_{index}"] for index in range(1, 5)]
        heatmap_row = {
            column: get_abundance_value(row, abundance_column_lookup, column)
            for column in abundance_columns
        }
        rows.append(heatmap_row)
        labels.append(format_mixture_label(samples, true_proportions, species_to_order))

    return pl.DataFrame(rows).select(abundance_columns).to_pandas(), labels


def get_mixture_id_for_components(metadata_df: pl.DataFrame, components: list[str]) -> int:
    target_components = {normalise_species_name(component).lower() for component in components}
    matches = []
    for row in metadata_df.to_dicts():
        mixture_components = {
            normalise_species_name(str(row[f"sample_{index}"])).lower()
            for index in range(1, 5)
        }
        if mixture_components == target_components:
            matches.append(int(row["mixture_id"]))

    if len(matches) != 1:
        raise ValueError(f"Expected exactly one simulated mixture for {components}, found {len(matches)}")
    return matches[0]


def prepare_selected_in_silico_mixture_heatmap_group(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    abundance_columns: list[str],
    selected_mixtures: list[dict],
) -> tuple:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    mixture_ids = [
        get_mixture_id_for_components(metadata_df, mixture["components"])
        for mixture in selected_mixtures
    ]
    mixture_df = prepare_in_silico_mixture_df(abundance_df, metadata_file)
    mixture_order_df = pl.DataFrame({"mixture_id": mixture_ids, "mixture_order": list(range(len(mixture_ids)))})
    selected_mixture_df = (
        mixture_df
        .filter(pl.col("mixture_id").is_in(mixture_ids))
        .join(mixture_order_df, on="mixture_id", how="inner")
        .sort("mixture_order")
    )
    abundance_column_lookup = {
        normalise_species_name(column).lower(): column
        for column in abundance_df.columns
        if column != "sample"
    }

    rows = []
    for row in selected_mixture_df.to_dicts():
        heatmap_row = {
            column: get_abundance_value(row, abundance_column_lookup, column)
            for column in abundance_columns
        }
        rows.append(heatmap_row)

    labels = [mixture["label"] for mixture in selected_mixtures]
    return pl.DataFrame(rows).select(abundance_columns).to_pandas(), labels


def prepare_in_silico_mixture_boxplot_data(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    abundance_columns: list[str],
    species_to_order: dict[str, str],
) -> pl.DataFrame:
    plot_column_lookup = {
        normalise_species_name(column).lower(): column
        for column in abundance_columns
    }
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    error_df = compute_mixture_abundance_error_records(
        abundance_df,
        metadata_df,
        tool="abundance_within_filter_pass_reads",
        relative_error_margin=0.10,
    )
    rows = []
    for row in error_df.to_dicts():
        species_key = normalise_species_name(row["species"]).lower()
        if species_key not in plot_column_lookup:
            continue
        plot_species_name = normalise_species_name(plot_column_lookup[species_key])
        rows.append(
            {
                "mixture_id": row["mixture_id"],
                "species": plot_species_name,
                "species_label": species_to_order[species_key],
                "expected_abundance": row["expected_abundance"],
                "estimated_abundance": row["estimated_abundance"],
                "kma_rmse": row["mixture_rmse"],
            }
        )
    return pl.DataFrame(rows)


def prepare_kmer_matching_error_data(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    kmer_matching_file: str,
) -> pl.DataFrame:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    error_by_species = (
        compute_mixture_abundance_error_records(
            abundance_df,
            metadata_df,
            tool="abundance_within_filter_pass_reads",
        )
        .with_columns(
            pl.col("species")
            .map_elements(lambda value: normalise_species_name(str(value)).lower(), return_dtype=pl.Utf8)
            .alias("species_key")
        )
        .group_by(["species_key", "species"])
        .agg(pl.col("error_estimated_minus_expected").mean().alias("mean_signed_error"))
    )

    kmer_matching_df = (
        pl.read_csv(kmer_matching_file, separator="\t")
        .with_columns(
            pl.col("species")
            .map_elements(lambda value: normalise_species_name(str(value)).lower(), return_dtype=pl.Utf8)
            .alias("species_key"),
            (
                pl.col("correct_only_weighted_proportion")
                + pl.col("shared_correct_and_incorrect_weighted_proportion")
            ).alias("correct_weighted_proportion"),
        )
    )

    return (
        error_by_species
        .join(kmer_matching_df, on="species_key", how="inner")
        .with_columns(pl.col("species").alias("plot_label"))
        .sort("species")
    )


def prepare_non_confident_assignment_error_data(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    single_species_abundance_file: str,
) -> pl.DataFrame:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    error_by_species = (
        compute_mixture_abundance_error_records(
            abundance_df,
            metadata_df,
            tool="abundance_within_filter_pass_reads",
        )
        .with_columns(
            pl.col("species")
            .map_elements(lambda value: normalise_species_name(str(value)).lower(), return_dtype=pl.Utf8)
            .alias("species_key")
        )
        .group_by(["species_key", "species"])
        .agg(pl.col("error_estimated_minus_expected").mean().alias("mean_signed_error"))
    )
    non_confident_by_species = (
        pl.read_csv(single_species_abundance_file, separator="\t")
        .select(["sample", "non_confidently_assigned"])
        .with_columns(
            pl.col("sample")
            .map_elements(lambda value: normalise_species_name(str(value)).lower(), return_dtype=pl.Utf8)
            .alias("species_key")
        )
    )
    return (
        error_by_species
        .join(non_confident_by_species, on="species_key", how="inner")
        .with_columns(pl.col("species").alias("plot_label"))
        .sort("species")
    )


def prepare_in_silico_mixture_error_summary(
    abundance_df: pl.DataFrame,
    metadata_file: str,
    tool: str = "abundance_within_filter_pass_reads",
    relative_error_margin: float = 0.10,
) -> pl.DataFrame:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    error_df = compute_mixture_abundance_error_records(
        abundance_df,
        metadata_df,
        tool=tool,
        relative_error_margin=relative_error_margin,
    )
    return summarise_mixture_abundance_errors(error_df)


def write_in_silico_mixture_error_summary(summary_df: pl.DataFrame, output_file: str) -> None:
    summary_df.write_csv(output_file, separator="\t")


def prepare_in_silico_mixture_violin_data(
    joined_assignment_file: str,
    metadata_file: str,
) -> pl.DataFrame:
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    n_components = len([column for column in metadata_df.columns if column.startswith("sample_")])
    joined_df = (
        pl.read_csv(joined_assignment_file, separator="\t")
        .with_columns(
            pl.col("target_id")
            .map_elements(lambda value: normalise_species_name(return_species_from_target_id(value)), return_dtype=pl.Utf8)
            .alias("species"),
            pl.col("kma:readCount").fill_null(0).cast(pl.Float64),
            pl.col("kallisto:est_counts_mean").fill_null(0).cast(pl.Float64),
            pl.col("sample").str.extract(r"mixture_([0-9]+)", 1).cast(pl.Int64).alias("mixture_id"),
        )
    )

    rows = []
    for tool, count_column in [
        ("KMA", "kma:readCount"),
        ("kallisto", "kallisto:est_counts_mean"),
    ]:
        count_df = (
            joined_df
            .group_by(["sample", "mixture_id", "species"])
            .agg(pl.col(count_column).sum().alias("count"))
        )
        abundance_df = count_df.join(metadata_df, on="mixture_id", how="inner")

        for row in abundance_df.group_by(["sample", "mixture_id"]).agg(pl.all()).to_dicts():
            true_species = [
                normalise_species_name(str(row[f"sample_{index}"][0])).lower()
                for index in range(1, n_components + 1)
            ]
            true_species_set = set(true_species)
            species_to_count = {
                normalise_species_name(species).lower(): count
                for species, count in zip(row["species"], row["count"])
            }
            true_counts = [species_to_count.get(species, 0.0) for species in true_species]
            true_count_total = sum(true_counts)
            estimated_true_abundances = [count / true_count_total for count in true_counts] if true_count_total > 0 else [0.0] * len(true_counts)
            expected_abundances = [float(row[f"proportion_{index}"][0]) for index in range(1, n_components + 1)]
            rmse = float(np.sqrt(np.mean([(estimated - expected) ** 2 for estimated, expected in zip(estimated_true_abundances, expected_abundances)])))
            total_count = sum(species_to_count.values())
            off_target_count = sum(count for species, count in species_to_count.items() if species not in true_species_set)
            proportion_off_target = off_target_count / total_count if total_count > 0 else 0.0
            rows.append(
                {
                    "tool": tool,
                    "sample": row["sample"],
                    "mixture_id": row["mixture_id"],
                    "rmse": rmse,
                    "proportion_off_target": proportion_off_target,
                }
            )
    return pl.DataFrame(rows).sort(["tool", "mixture_id"])


def remove_empty_heatmap_columns(heatmap_data: dict, groups: list[str]) -> dict:
    column_metadata = heatmap_data["column_metadata"]
    columns_to_keep = []
    for column_index, metadata in enumerate(column_metadata):
        keep_column = any(
            heatmap_data[group][0].iloc[:, column_index].fillna(0).sum() > 0
            for group in groups
        )
        if not keep_column and "simulated_mixture_boxplot_data" in heatmap_data:
            species = metadata["species"]
            keep_column = (
                heatmap_data["simulated_mixture_boxplot_data"]
                .filter(pl.col("species") == species)
                .height
                > 0
            )
        if keep_column:
            columns_to_keep.append(column_index)

    filtered_heatmap_data = {}
    for key, value in heatmap_data.items():
        if key == "column_metadata":
            filtered_heatmap_data[key] = [column_metadata[index] for index in columns_to_keep]
        elif key == "simulated_mixture_boxplot_data":
            species_to_keep = [column_metadata[index]["species"] for index in columns_to_keep]
            filtered_heatmap_data[key] = value.filter(pl.col("species").is_in(species_to_keep))
        elif key in groups:
            data, labels = value
            filtered_heatmap_data[key] = (data.iloc[:, columns_to_keep], labels)
        else:
            filtered_heatmap_data[key] = value
    return filtered_heatmap_data


def prepare_control_heatmap_data(
    abundance_df: pl.DataFrame,
    annotations_df: pl.DataFrame,
    in_silico_abundance_df: pl.DataFrame | None = None,
    in_silico_metadata_file: str | None = None,
    n_in_silico_mixtures: int = 18,
    exclude_sample_ids: list[str] | None = None,
    exclude_species_names: list[str] | None = None,
    selected_simulated_mixtures: list[dict] | None = None,
    drop_empty_columns: bool = True,
) -> dict:
    if exclude_sample_ids is not None:
        annotations_df = annotations_df.filter(~pl.col("sample").is_in(exclude_sample_ids))
    species_annotations = annotations_df.filter(pl.col("type") == "species")
    species_to_order = {}
    for content, label in zip(species_annotations["content"], species_annotations["label"]):
        species_name = normalise_species_name(str(content)).lower()
        if species_name not in species_to_order:
            species_to_order[species_name] = str(label).strip()
    if "culex pipiens" not in species_to_order and "culex australicus" in species_to_order:
        species_to_order["culex pipiens"] = species_to_order["culex australicus"]

    abundance_cols = [col for col in abundance_df.columns if col != "sample"]
    if exclude_species_names is not None:
        excluded_species_names = {normalise_species_name(species).lower() for species in exclude_species_names}
        abundance_cols = [
            col for col in abundance_cols
            if normalise_species_name(col).lower() not in excluded_species_names
        ]
    abundance_cols_ordered = sorted(
        abundance_cols,
        key=lambda col: (
            GENUS_ORDER.get(get_genus(col), len(GENUS_ORDER)),
            get_species_sort_order(col, species_to_order),
        ),
    )
    abundance_df = abundance_df.select(["sample"] + abundance_cols_ordered)
    column_metadata = [
        {
            "species": normalise_species_name(col),
            "genus": get_genus(col),
            "label": get_species_plot_label(col, species_to_order),
        }
        for col in abundance_cols_ordered
    ]

    annotations_df = annotations_df.with_columns(
        pl.when(pl.col("sample").is_in(["50-50-MultiAust_S44", "75-25-MultiAust_S45", "95-05-MultiAust_S46"]))
        .then(pl.col("label").map_elements(fix_COSMOS6_label, return_dtype=pl.Utf8))
        .otherwise(pl.col("label"))
        .alias("label")
    )
    annotations_df = annotations_df.with_columns(
        pl.when(pl.col("type") == "species")
        .then(
            pl.struct(["content", "label"]).map_elements(
                lambda row: format_species_sample_label(str(row["content"]), str(row["label"])),
                return_dtype=pl.Utf8,
            )
        )
        .otherwise(pl.col("label"))
        .alias("plot_label")
    )

    species_rows = sorted(
        annotations_df
        .filter((pl.col("type") == "species") & (~pl.col("label").is_in(["L1", "L2"])))
        .to_dicts(),
        key=lambda row: (
            GENUS_ORDER.get(get_genus(str(row["content"])), len(GENUS_ORDER)),
            get_species_sort_order(str(row["content"]), species_to_order),
        ),
    )
    morphology_rows = sorted(
        annotations_df
        .filter((pl.col("type") == "species") & (pl.col("label").is_in(["L1", "L2"])))
        .to_dicts(),
        key=lambda row: letter_to_number(str(row["label"])),
    )
    species_samples = [row["sample"] for row in species_rows]
    morphology_samples = [row["sample"] for row in morphology_rows]
    groups = {
        "simulated mixture": annotations_df.filter(pl.col("type") == "simulated mixture")["sample"].to_list(),
        "species": species_samples,
        "morphology test": morphology_samples,
        "laboratory mixture": annotations_df.filter(pl.col("type") == "laboratory mixture")["sample"].to_list(),
        "laboratory serial dilution": annotations_df.filter(pl.col("type") == "laboratory serial dilution")["sample"].to_list(),
    }

    sample_to_label = dict(zip(annotations_df["sample"], annotations_df["plot_label"]))
    for index, sample in enumerate(groups["morphology test"], start=1):
        sample_to_label[sample] = f"Mixture {index}"
    labels = [r"$\it{Ae.proc/Ae.mult/Cx.siti/Cx.pipi}$", r"$\it{Cq.line/Ae.noto/Ae.vigi/Cx.annu}$"]
    for sample, label in zip(groups["laboratory mixture"], labels):
        sample_to_label[sample] = label

    raw_sample_to_label = dict(zip(annotations_df["sample"], annotations_df["label"]))
    for sample in groups["laboratory serial dilution"]:
        sample_to_label[sample] = format_ratio_label(raw_sample_to_label.get(sample, sample))

    heatmap_data = {
        group: get_group_data_and_labels(abundance_df, samples, sample_to_label)
        for group, samples in groups.items()
    }
    if in_silico_abundance_df is not None and in_silico_metadata_file is not None:
        if selected_simulated_mixtures is not None:
            heatmap_data["simulated mixture"] = prepare_selected_in_silico_mixture_heatmap_group(
                in_silico_abundance_df,
                in_silico_metadata_file,
                abundance_cols_ordered,
                selected_simulated_mixtures,
            )
        else:
            heatmap_data["simulated mixture"] = prepare_in_silico_mixture_heatmap_group(
                in_silico_abundance_df,
                in_silico_metadata_file,
                abundance_cols_ordered,
                species_to_order,
                n_mixtures=n_in_silico_mixtures,
            )
        heatmap_data["simulated_mixture_boxplot_data"] = prepare_in_silico_mixture_boxplot_data(
            in_silico_abundance_df,
            in_silico_metadata_file,
            abundance_cols_ordered,
            species_to_order,
        )
    heatmap_data["column_metadata"] = column_metadata
    if drop_empty_columns:
        heatmap_data = remove_empty_heatmap_columns(heatmap_data, list(groups.keys()))
    return heatmap_data
