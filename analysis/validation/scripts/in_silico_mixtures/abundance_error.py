# compute RMSE and proportion off-target for KMA and kallisto outputs

import json
import numpy as np
import polars as pl

def compute_RMSE(true_proportions: list, estimated_proportions: list) -> float:
    true_proportions = np.array(true_proportions)
    estimated_proportions = np.array(estimated_proportions)
    rmse = np.sqrt(np.mean((true_proportions - estimated_proportions) ** 2))
    return rmse

def compute_proportion_reads_off_target(df: pl.DataFrame, true_species: list, count_column: str) -> float:
    total_reads = df[count_column].sum()
    off_target_reads = df.filter(~pl.col("species").is_in(true_species))[count_column].sum()
    proportion_off_target = off_target_reads / total_reads if total_reads > 0 else 0
    return proportion_off_target

def compute_abundance_error(mixture_info_file: str, abundance_estimates_file: str, output_file: str) -> None:
    # parse files
    mixture_df = pl.read_csv(mixture_info_file, separator='\t')
    abundance_estimates_df = pl.read_csv(abundance_estimates_file, separator='\t')
    n_components = len([i for i in mixture_df.columns if "proportion" in i])
    mixture_ids = mixture_df['mixture_id'].to_list()
    results = []
    for mixture in mixture_ids:
        true_proportions = [mixture_df.filter(pl.col("mixture_id") == mixture)[f"proportion_{i}"].item() for i in range(1, n_components + 1)]
        true_species = [mixture_df.filter(pl.col("mixture_id") == mixture)[f"sample_{i}"].item() for i in range(1, n_components + 1)]
        kallisto_estimated_counts = [abundance_estimates_df.filter(pl.col("mixture_id") == mixture).filter(pl.col("species") == species)['kallisto_estimated_counts'].sum() for i, species in enumerate(true_species, start=1)]
        kma_estimated_counts = [abundance_estimates_df.filter(pl.col("mixture_id") == mixture).filter(pl.col("species") == species)['kma_read_count'].sum() for i, species in enumerate(true_species, start=1)]
        kallisto_estimated_proportions = list(np.array(kallisto_estimated_counts) / np.sum(np.array(kallisto_estimated_counts)))
        kma_estimated_proportions = list(np.array(kma_estimated_counts) / np.sum(np.array(kma_estimated_counts)))
        kallisto_rmse = compute_RMSE(true_proportions, kallisto_estimated_proportions)
        kma_rmse = compute_RMSE(true_proportions, kma_estimated_proportions)
        mixture_estimates = abundance_estimates_df.filter(pl.col("mixture_id") == mixture)
        kma_proportion_off_target = compute_proportion_reads_off_target(mixture_estimates, true_species, "kma_read_count")
        kallisto_proportion_off_target = compute_proportion_reads_off_target(mixture_estimates, true_species, "kallisto_estimated_counts")
        results.append({
            "mixture_id": mixture,
            "true_proportions": json.dumps(true_proportions),
            "kma_estimated_proportions": json.dumps(kma_estimated_proportions),
            "kallisto_estimated_proportions": json.dumps(kallisto_estimated_proportions),
            "kma_total_reads": abundance_estimates_df.filter(pl.col("mixture_id") == mixture)['kma_read_count'].sum(),
            "kallisto_total_reads": abundance_estimates_df.filter(pl.col("mixture_id") == mixture)['kallisto_estimated_counts'].sum(),
            "kma_on_target_reads": abundance_estimates_df.filter(pl.col("mixture_id") == mixture).filter(pl.col("species").is_in(true_species))['kma_read_count'].sum(),
            "kallisto_on_target_reads": abundance_estimates_df.filter(pl.col("mixture_id") == mixture).filter(pl.col("species").is_in(true_species))['kallisto_estimated_counts'].sum(),
            "kma_rmse": kma_rmse,
            "kallisto_rmse": kallisto_rmse,
            "kma_proportion_off_target": kma_proportion_off_target,
            "kallisto_proportion_off_target": kallisto_proportion_off_target
        })
    results_df = pl.DataFrame(results)
    results_df.write_csv(output_file, separator='\t')
