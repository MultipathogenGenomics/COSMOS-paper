from importlib import import_module
import pickle
import os
import polars as pl
import numpy as np

assess_enrichment = import_module("analysis.assess-enrichment.process")
abundance_validation = import_module("analysis.abundance-validation.process")
from ..common.abundance import read_abundance_for_sample_sets
from ..common.mosquito_ids import shorted_genus


def shorten_italic_species_label(label: str) -> str:
    label = str(label)
    mathtext_prefix = r"$\it{"
    if not (label.startswith(mathtext_prefix) and label.endswith("}$")):
        return label
    species_name = label[len(mathtext_prefix):-2].replace(r"\ ", " ")
    shortened_name = shorted_genus(species_name)
    return mathtext_prefix + shortened_name.replace(" ", r"\,") + "}$"

def make_mt_genome_coverage_file(out_path: str):
    # inputs
    genomes_config = "config/mt_genomes.json"
    samples_config = "config/single_species_controls.trimmed.json"
    genus_defaults = {
        "Aedes": "data/mt_genomes/Aedes_vigilax_mt_genome.fasta",
        "Culex": "data/mt_genomes/Culex_annulirostris_mt_genome.fasta",
        "Anopheles": "data/mt_genomes/Anopheles_annulipes_mt_genome.fasta",
        "Mansonia": "data/mt_genomes/Mansonia_uniformis_mt_genome.fasta",
        "Coquillettidia": "data/mt_genomes/Coquillettidia_linealis_mt_genome.fasta",
        "Ochlerotatus": "data/mt_genomes/Aedes_vigilax_mt_genome.fasta"
    }
    # collate data
    coverage_df = assess_enrichment.concat_coverage_data(samples_config, genomes_config, "results/assess-enrichment/coverage", genus_defaults)
    coverage_df = assess_enrichment.add_alignment_positions_to_coverage_data(coverage_df, "results/assess-enrichment/alignments/references_aligned.fasta")
    # write the coverage data to a file
    coverage_df.write_csv(out_path)

def prepare_heatmap_data(out_path: str):
    estimate_abundance_dir = "results/estimate-abundance/COSMOS-plus-SS-pools"
    control_sample_sets = ["COSMOS2_controls", "COSMOS5_controls", "COSMOS6_controls"]
    in_silico_sample_sets = ["in_silico_mixtures"]
    abundance_filename = "abundance_within_filter_pass_reads.tsv"
    joined_assignment_filename = "kallisto_kma_outputs_joined.tsv"
    annotations_file = "data/COSMOS-controls/COSMOS_controls_pubfig_annotations.tsv"
    in_silico_metadata_file = "results/in-silico-mixtures/mixture_fastq/mixture_metadata.tsv"
    single_species_fastq_config = "config/single_species_controls.mosquito_cox1_extracted_fastq.json"
    assignment_reference_file = (
        "results/reference-set/COSMOS-plus-SS-pools/"
        "COSMOS-plus-SS-pools.representatives.fasta"
    )
    single_species_abundance_file = (
        f"{estimate_abundance_dir}/single-species-controls-mosquito-cox1-extracted/"
        "abundance_within_all_assigned_reads.tsv"
    )
    selected_simulated_mixtures = [
        {
            "label": r"$\it{Ae.proc/Ae.mult/Cx.siti/Cx.pipi}$", 
            "components": ["Aedes procax", "Aedes multiplex", "Culex sitiens", "Culex pipiens"],
        },
        {
            "label": r"$\it{Cq.line/Ae.noto/Ae.vigi/Cx.annu}$",
            "components": ["Coquillettidia linealis", "Aedes notoscriptus", "Aedes vigilax", "Culex annulirostris"],
        },
    ]

    abundance_df = read_abundance_for_sample_sets(
        estimate_abundance_dir,
        control_sample_sets,
        abundance_filename,
    ).drop("batch")
    in_silico_abundance_df = read_abundance_for_sample_sets(
        estimate_abundance_dir,
        in_silico_sample_sets,
        abundance_filename,
    ).drop("batch")
    annotations_df = abundance_validation.load_control_annotations(annotations_file)

    main_heatmap_data = abundance_validation.prepare_control_heatmap_data(
        abundance_df,
        annotations_df,
        in_silico_abundance_df,
        in_silico_metadata_file,
        n_in_silico_mixtures=18,
        exclude_sample_ids=[
            "50-50-MultiAust_S44",
            "75-25-MultiAust_S45",
            "95-05-MultiAust_S46",
        ],
        selected_simulated_mixtures=selected_simulated_mixtures,
    )
    for metadata in main_heatmap_data["column_metadata"]:
        metadata["display_species"] = shorted_genus(metadata["species"])
    species_data, species_labels = main_heatmap_data["species"]
    main_heatmap_data["species"] = (
        species_data,
        [shorten_italic_species_label(label) for label in species_labels],
    )

    with open(out_path, "wb") as handle:
        pickle.dump(main_heatmap_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return main_heatmap_data

def combine_cryptic_subpop_data(cryptic_dir: str, output_dir: str, datasets: list[str]) -> pl.DataFrame:
    subpop_data = []
    for dataset in ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]:
        subpop_file = os.path.join(cryptic_dir, f"{dataset}.Cx_annulirostris_collated_results.tsv")
        subpop_df = pl.read_csv(subpop_file, separator="\t")
        subpop_df = subpop_df.with_columns(pl.lit(dataset).alias("dataset"))
        subpop_data.append(subpop_df)
    combined_subpop_df = pl.concat(subpop_data, how="diagonal")
    return combined_subpop_df

def order_results_data_for_plotting(df: pl.DataFrame) -> pl.DataFrame:
    # make a col of Rural, Urban, MEt_fringe and Met
    MET_CODES = [1, 2, 3]
    URBAN_CODES = [4, 5]
    FRINGE_CODES = [6, 7]
    RURAL_CODES = [8, 9, 10, 11]
    df = df.with_columns(
        pl.when(pl.col("council_olg_code").is_in(RURAL_CODES))
        .then(pl.lit("Rural"))
        .when(pl.col("council_olg_code").is_in(FRINGE_CODES))
        .then(pl.lit("Met_Fringe"))
        .when(pl.col("council_olg_code").is_in(URBAN_CODES))
        .then(pl.lit("Urban"))
        .when(pl.col("council_olg_code").is_in(MET_CODES))
        .then(pl.lit("Met"))
        .otherwise(None)
        .alias("urban_rural")
    )
    # combine coastal and urban_rural into a single column
    df = df.with_columns(
        pl.concat_str(
            [
                pl.col("Coastal").fill_null("Unknown"),
                pl.col("urban_rural").fill_null("Unknown"),
            ],
            separator=" - ",
        ).alias("Coastal_urban_rural")
    )
    # Order the categories by their average longitude within each coastal
    # group, then order individual samples from west to east.
    df = df.with_columns(
        pl.col("longitude")
        .mean()
        .over(["Coastal", "urban_rural"])
        .alias("avg_longitude")
    )
    df = df.sort(
        ["Coastal", "avg_longitude", "urban_rural", "longitude"],
        descending=[False, False, False, False],
    )
    # print average longitude for each coastal_urban_rural group
    avg_longitudes = (
        df.select(["Coastal", "urban_rural", "avg_longitude"])
        .unique()
        .sort(["Coastal", "avg_longitude"], descending=[False, False])
    )
    print("Average longitude for each Coastal_urban_rural group:")
    print(avg_longitudes)
    n_blocks_of_urban_rural = (
        df.select(
            pl.col("urban_rural").ne_missing(pl.col("urban_rural").shift()).sum()
        )
        .item()
    )
    print(f"Number of blocks of sorted samples: {n_blocks_of_urban_rural}")
    return df

def add_whether_east_or_west_of_GDR(df):
    gdr_coords = np.array([
        [152.0310561857285, -29.05884438306106],
        [151.8419028204735, -29.74480814986759],
        [150.9466485608887, -32.0576776944991],
        [150.3174422786903, -33.71480730074187],
        [149.7800386106203, -34.75207284189212],
        [149.2961935232750, -35.36848469567541],
        [148.8017118226620, -36.21266565329890],
        [148.5659720513401, -36.82561305143643],
    ])
    # np.interp requires its x values—latitude here—to increase.
    order = np.argsort(gdr_coords[:, 1])
    gdr_latitudes = gdr_coords[order, 1]
    gdr_longitudes = gdr_coords[order, 0]
    sample_latitudes = df["latitude"].cast(pl.Float64).to_numpy()
    # Estimate the GDR longitude at each sample's latitude.
    gdr_longitude_at_sample = np.interp(
        sample_latitudes,
        gdr_latitudes,
        gdr_longitudes,
        left=np.nan,               # south of the supplied GDR coordinates
        right=gdr_longitudes[-1],  # north: use the northernmost longitude
    )
    df = df.with_columns(
        pl.Series(
            name="GDR_longitude",
            values=gdr_longitude_at_sample,
        )
    )
    df = df.with_columns(
        pl.when(
            pl.col("longitude").is_null()
            | pl.col("latitude").is_null()
            | pl.col("GDR_longitude").is_nan()
        )
        .then(pl.lit(None, dtype=pl.String))
        .when(pl.col("longitude") > pl.col("GDR_longitude"))
        .then(pl.lit("E"))
        .otherwise(pl.lit("W"))
        .alias("E_or_W_of_GDR")
    )
    n_blocks_of_EW = (
        df.select(
            pl.col("E_or_W_of_GDR").ne_missing(pl.col("E_or_W_of_GDR").shift()).sum()
        )
        .item()
    )
    print(f"Number of blocks of E/W samples: {n_blocks_of_EW}")
    return df

def combine_data_for_results_figure(
        metadata_file: str,
        SDI_results_file: str,
        mosquito_abundance_file: str,
        virus_nc2_file: str,
        map_file: str
) -> pl.DataFrame:
    virus_acronyms = {
        "stratfordvirus": "STRV",
        "rossrivervirus": "RRV",
        "murrayvalleyencephalitisvirus": "MVEV",
        "kunjinvirus": "KUN",
        "kokoberavirus": "KOK",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "barmahforestvirus": "BFV",
        "sindbisvirus": "SINDV",
        "umbrevirus": "LSV",
    }
    metadata_df = pl.read_csv(metadata_file, separator="\t")
    metadata_df = metadata_df.select([col for col in metadata_df.columns if col not in virus_acronyms.values()])
    SDI_df = pl.read_csv(SDI_results_file, separator="\t")
    mosquito_abundance_df = pl.read_csv(mosquito_abundance_file)
    virus_nc2_df = pl.read_csv(virus_nc2_file, separator="\t")
    virus_nc2_df = virus_nc2_df.select(["sampleid"] + list(virus_acronyms.keys()))
    virus_nc2_df = virus_nc2_df.rename(virus_acronyms).rename({"sampleid": "sample"})
    combined_df = metadata_df.join(SDI_df, on="sample", how="left")
    combined_df = combined_df.join(mosquito_abundance_df, on="sample", how="left")
    combined_df = combined_df.join(virus_nc2_df, on="sample", how="left")
    combined_df = combined_df.with_columns(
        pl.when(pl.col("distance_to_coast_km") < 20)
        .then(pl.lit(1))
        .when(pl.col("distance_to_coast_km") >= 20)
        .then(pl.lit(0))
        .otherwise(None)
        .alias("Coastal")
    )
    combined_df = combined_df.rename(
        {"Culex annulirostris": "Culex annulirostris proportion"}
    )

    # going to make divide Cx annulirostris into two cols, S-AUS and AUS
    cryptic_process = import_module("analysis.cryptic-subpopulations.process")

    panel_config = {
        "label": "Culex annulirostris",
        "species_abundance_column": "Culex annulirostris",
        "file_template": "data/cryptic_subpops/{dataset}.Cx_annulirostris_collated_results.tsv",
        "ratio_columns": ["avg_subtype_ratio", "non_subtype_ratio"],
        "legend_labels": ["S-AUS", "AUS"],
        "min_count_threshold": 20,
        "complement_ratio_column": "avg_subtype_ratio",
        "complement_output_column": "non_subtype_ratio",
        "draw_type": "gradient"
    }
    subpopulation_data = cryptic_process.prepare_cryptic_subpopulation_data(
        map_file,
        metadata_file,
        [panel_config],
        ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"],
    )
    cx_panel_points = subpopulation_data["panels"][0]["points"]
    cx_subtype_ratios_df = pl.from_pandas(
        cx_panel_points[
            ["sample", "avg_subtype_ratio", "non_subtype_ratio"]
        ].drop_duplicates(subset=["sample"])
    ).select(
        pl.col("sample").cast(pl.String),
        pl.col("avg_subtype_ratio").cast(pl.Float64),
        pl.col("non_subtype_ratio").cast(pl.Float64),
    )
    combined_df = combined_df.join(cx_subtype_ratios_df, on="sample", how="left")
    combined_df = combined_df.with_columns(
        (
            pl.col("Culex annulirostris proportion")
            * pl.col("avg_subtype_ratio")
        ).alias("Culex annulirostris S-AUS"),
        (
            pl.col("Culex annulirostris proportion")
            * pl.col("non_subtype_ratio")
        ).alias("Culex annulirostris AUS"),
    ).drop("non_subtype_ratio")

    print(f"Combined data for results figure: {combined_df.shape[0]} samples, {combined_df.shape[1]} columns")
    ordered_df = order_results_data_for_plotting(combined_df)
    ordered_df = add_whether_east_or_west_of_GDR(ordered_df)
    return ordered_df

def parse_annulirostris_mt_coverage(output_file_path: str, metadata_file: str = None) -> pl.DataFrame:
    coverage_data = []
    for sample_set in ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"]:
        coverage_dir = f"/Users/afor0966/Documents/code/multipathogen_genomics/results/{sample_set}/aligned/Cx_annulirostris_mt"
        for sample in os.listdir(coverage_dir):
            sample_dir = os.path.join(coverage_dir, sample)
            for filename in os.listdir(sample_dir):
                if not filename.endswith("_coverage.tsv.gz"):
                    continue
                sample_name = filename.replace("_coverage.tsv.gz", "")
                file_path = os.path.join(sample_dir, filename)
                df = pl.read_csv(file_path, separator='\t', has_header=False)
                df.columns = ["contig", "position", "coverage"]
                df = df.with_columns(
                pl.lit(sample_name).alias("sample"),
                pl.lit(sample_set).alias("sample_set"),
                pl.lit("Culex annulirostris").alias("ref_genome"),
                pl.lit(True).alias("ref_is_species"),
            )
            cox1_start = 1430
            cox1_end = 2966
            df = df.with_columns(
                pl.when(
                    (pl.col("position") >= cox1_start)
                    & (pl.col("position") <= cox1_end)
                )
                .then(pl.lit(True))
                .otherwise(pl.lit(False))
                .alias("is_cox1")
            )
            coverage_data.append(df)
    mito_df = pl.concat(coverage_data, how="diagonal")
    # now filter to just those in metadata
    if metadata_file is not None:
        metadata_df = pl.read_csv(metadata_file, separator="\t")
        mito_df = mito_df.join(metadata_df.select(["sample"]), on="sample", how="inner")
    mito_df.write_csv(output_file_path)
    # also compute median and mean
    median_coverage = assess_enrichment.compute_median_mean_in_cox1(mito_df)
    median_coverage.write_csv(output_file_path.replace(".csv", "_median_mean_in_cox1.csv"))
