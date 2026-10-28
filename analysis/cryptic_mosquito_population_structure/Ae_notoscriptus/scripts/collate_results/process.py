import polars as pl

def get_Ae_notoscriptus_cluster1_to_3_ratios(cluster1_vs_2or3_file: str, cluster2_vs_3_file: str, outpath: str) -> pl.DataFrame:
    cluster1_vs_2or3_df = pl.read_csv(cluster1_vs_2or3_file, separator='\t')
    cluster2_vs_3_df = pl.read_csv(cluster2_vs_3_file, separator='\t')
    cluster1_vs_2or3_df = cluster1_vs_2or3_df.rename({"avg_subtype_ratio": "proportion_cluster1"})
    combined = cluster1_vs_2or3_df.join(
        cluster2_vs_3_df.select(["sample", "avg_subtype_ratio"]),
        on="sample",
        how="left"
    )
    combined = combined.with_columns(
        ratio_cluster2_given_not1=pl.col("avg_subtype_ratio").fill_null(0.0),
        remainder_after_cluster1=(1 - pl.col("proportion_cluster1"))
    ).with_columns(
        proportion_cluster2=pl.col("ratio_cluster2_given_not1") * pl.col("remainder_after_cluster1"),
        proportion_cluster3=(1 - pl.col("ratio_cluster2_given_not1")) * pl.col("remainder_after_cluster1")
    ).select(["sample", "proportion_cluster1", "proportion_cluster2", "proportion_cluster3"])
    # round all to 3dp
    combined = combined.with_columns(
        pl.col("proportion_cluster1").round(3),
        pl.col("proportion_cluster2").round(3),
        pl.col("proportion_cluster3").round(3)
    )
    combined.write_csv(outpath, separator='\t')
    return combined
    
def combine_chunks_of_linkage_results(input_files: list) -> pl.DataFrame:
    dfs = []
    for i in input_files:
        try:
            dfs.append(pl.read_csv(i, separator='\t'))
        except pl.exceptions.NoDataError:
            continue
    combined_df = pl.concat(dfs)
    # get rid of perfectly duplicated rows
    combined_df = combined_df.unique()
    return combined_df