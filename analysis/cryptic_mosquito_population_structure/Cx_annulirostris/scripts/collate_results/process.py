import polars as pl
import numpy as np

def compute_R(row):
    AB = row["AB"]
    aB = row["aB"]
    Ab = row["Ab"]
    ab = row["ab"]
    numerator = AB * ab - aB * Ab
    denominator = ((AB + aB) * (Ab + ab) * (AB + Ab) * (aB + ab)) ** 0.5
    if denominator == 0:
        return np.nan
    else:
        return numerator / denominator
    
def compute_D_prime(row):
    AB = row["AB"]
    aB = row["aB"]
    Ab = row["Ab"]
    ab = row["ab"]
    if sum([AB, aB, Ab, ab]) == 0:
        return np.nan
    pA = (AB + Ab) / (AB + aB + Ab + ab)
    pB = (AB + aB) / (AB + aB + Ab + ab)
    D = (AB / (AB + aB + Ab + ab)) - (pA * pB)
    if D >= 0:
        D_max = min(pA * (1 - pB), (1 - pA) * pB)
    else:
        D_max = min(pA * pB, (1 - pA) * (1 - pB))
    if D_max == 0:
        return np.nan
    else:
        return D / D_max

def collate_count_and_linkage_results_by_sample(counts_file: str, linkage_file: str, output_file: str, count_threshold_for_averages: int=20) -> None:
    # parse data
    counts_df = pl.read_csv(counts_file, separator='\t')
    linkage_df = pl.read_csv(linkage_file, separator='\t')
    # get samples
    samples = set(counts_df["sample"].unique().to_list()) | set(linkage_df["sample"].unique().to_list())
    # add pearson's R and D prime to linkage results
    Dprimes = []
    Rs = []
    for row in linkage_df.iter_rows(named=True):
        Dprimes.append(compute_D_prime(row))
        Rs.append(compute_R(row))
    linkage_df = linkage_df.with_columns(
        pl.Series("D_prime", Dprimes),
        pl.Series("R", Rs)
    )
    # add column to linkage_df for whether to accept row for averages based on count threshold in counts_df
    accept_row = generate_correlation_filter(counts_df, linkage_df, count_threshold_for_averages)
    linkage_df = linkage_df.with_columns(pl.Series("accept_row", accept_row))
    # collate results for each sample
    results = []
    for sample in samples:
        sample_counts_df = counts_df.filter(pl.col("sample") == sample)
        sample_linkage_df = linkage_df.filter(pl.col("sample") == sample)
        # extract info from counts
        subtype_name = [ i for i in sample_counts_df.columns if "_subtype_bases" in i ][0].replace("_subtype_bases", "")
        min_count = sample_counts_df["total_counts_at_position"].min()
        max_count = sample_counts_df["total_counts_at_position"].max()
        min_subtype_count = sample_counts_df[f"{subtype_name}_subtype_count"].min()
        max_subtype_count = sample_counts_df[f"{subtype_name}_subtype_count"].max()
        min_non_subtype_count = sample_counts_df[f"{subtype_name}_non_subtype_count"].min()
        max_non_subtype_count = sample_counts_df[f"{subtype_name}_non_subtype_count"].max()
        avg_subtype_ratio = sample_counts_df.filter(pl.col("total_counts_at_position") >= count_threshold_for_averages)[f"{subtype_name}_ratio"].mean()
        std_subtype_ratio = sample_counts_df.filter(pl.col("total_counts_at_position") >= count_threshold_for_averages)[f"{subtype_name}_ratio"].std()
        n_snps_which_meet_threshold = sample_counts_df.filter(pl.col("total_counts_at_position") >= count_threshold_for_averages).shape[0]
        list_snps_which_meet_threshold = ",".join(sample_counts_df.filter(pl.col("total_counts_at_position") >= count_threshold_for_averages)["position_in_reference"].cast(pl.String).to_list())
        # extract info from linkage
        total_AB = sample_linkage_df["AB"].sum()
        total_aB = sample_linkage_df["aB"].sum()
        total_Ab = sample_linkage_df["Ab"].sum()
        total_ab = sample_linkage_df["ab"].sum()
        avg_D_prime = sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("D_prime").is_finite()))["D_prime"].mean()
        std_D_prime = sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("D_prime").is_finite()))["D_prime"].std()
        avg_R = sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("R").is_finite()))["R"].mean()
        std_R = sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("R").is_finite()))["R"].std()
        R_values = ",".join(sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("R").is_finite()))["R"].round(3).cast(pl.String).to_list())
        D_prime_values = ",".join(sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("D_prime").is_finite()))["D_prime"].round(3).cast(pl.String).to_list())
        n_snps_computed_linkage_for = sample_linkage_df.filter((pl.col("accept_row") == True) & (pl.col("D_prime").is_finite())).shape[0]
        # add to results
        row = {
            "sample": sample,
            "subtype_name": subtype_name,
            "min_count": min_count,
            "max_count": max_count,
            "min_subtype_count": min_subtype_count,
            "max_subtype_count": max_subtype_count,
            "min_non_subtype_count": min_non_subtype_count,
            "max_non_subtype_count": max_non_subtype_count,
            "avg_subtype_ratio": avg_subtype_ratio,
            "std_subtype_ratio": std_subtype_ratio,
            "n_snps_which_meet_threshold": n_snps_which_meet_threshold,
            "list_snps_which_meet_threshold": list_snps_which_meet_threshold,
            "total_AB": total_AB,
            "total_aB": total_aB,
            "total_Ab": total_Ab,
            "total_ab": total_ab,
            "avg_D_prime": avg_D_prime,
            "std_D_prime": std_D_prime,
            "avg_R": avg_R,
            "std_R": std_R,
            "R_values": R_values,
            "D_prime_values": D_prime_values,
            "n_pairs_computed_linkage_between": n_snps_computed_linkage_for
        }
        results.append(row)
    results = pl.DataFrame(results)
    # round floats to 3dp
    results = results.with_columns(
        pl.col("avg_subtype_ratio").round(3),
        pl.col("std_subtype_ratio").round(3),
        pl.col("avg_D_prime").round(3),
        pl.col("std_D_prime").round(3),
        pl.col("avg_R").round(3),
        pl.col("std_R").round(3)
    )
    results.write_csv(output_file, separator='\t')

def generate_correlation_filter(counts_df: pl.DataFrame, linkage_df: pl.DataFrame, count_threshold_for_averages: int) -> pl.DataFrame:
    accept_row = []
    for row in linkage_df.iter_rows(named=True):
        sample = row["sample"]
        sample_counts_df = counts_df.filter(pl.col("sample") == sample)
        locusA = row["locusA"]
        locusB = row["locusB"]
        count_A = sample_counts_df.filter(pl.col("position_in_reference") == int(locusA))["total_counts_at_position"].unique()[0]
        count_B = sample_counts_df.filter(pl.col("position_in_reference") == int(locusB))["total_counts_at_position"].unique()[0]
        if count_A >= count_threshold_for_averages and count_B >= count_threshold_for_averages:
            accept_row.append(True)
        else:
            accept_row.append(False)
    return np.array(accept_row)