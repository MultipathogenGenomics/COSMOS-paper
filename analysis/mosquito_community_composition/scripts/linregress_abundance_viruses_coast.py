from analysis.common.paths import RESULTS_DIR
from analysis.common.regression import fit_linear_regression, fit_logistic_regression
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova as skbio_permanova
from analysis.mosquito_community_composition.scripts.process import centre_log_ratio_transform

def parse_data():
    # abundance data
    mosquito_abundance_file = f"{RESULTS_DIR}/publication_figures/data/combined_abundance_data.tsv"
    abundance_df = pl.read_csv(mosquito_abundance_file)
    # virus data
    virus_data_file = f"{RESULTS_DIR}/final_castanet_data/COSMOS1-4.final_castanet_nc2.tsv"
    virus_df = pl.read_csv(virus_data_file, separator='\t')
    # metadata
    metadata = "data/COSMOS_metadata/COSMOS1-4.traps_combined_metadata.tsv"
    metadata_df = pl.read_csv(metadata, separator='\t')
    return abundance_df, virus_df, metadata_df

def add_species_count(abundance_df):
    mosquito_cols = [col for col in abundance_df.columns if col not in ['sample', 'batch']]
    return abundance_df.with_columns(
        pl.when(pl.all_horizontal([pl.col(col).is_null() for col in mosquito_cols]))
        .then(pl.lit(None))
        .otherwise(pl.sum_horizontal([pl.col(col) > 0 for col in mosquito_cols]))
        .alias('species_count_in_trap')
    )

def combine_data():

    abundance_df, virus_df, metadata_df = parse_data()

    # add a column called "species_count_in_trap", which is sum of species > 0 in each row
    abundance_df = add_species_count(abundance_df)

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
    # filter virus_df to only include columns in virus_acronyms
    virus_df = virus_df.select(['sampleid'] + [col for col in virus_df.columns if col in virus_acronyms.keys()])
    # rename virus_df columns to acronyms
    virus_df = virus_df.rename({col: virus_acronyms[col] for col in virus_df.columns if col in virus_acronyms.keys()})
    # where value in virus_df is >= 0.05, is 1, otherwise 0, across all virus columns
    virus_df = virus_df.with_columns([
        pl.when(pl.col(col) >= 0.05).then(1).otherwise(0).alias(col) for col in virus_df.columns if col not in ['sampleid']
    ])
    # metadata 

    metadata_df = metadata_df.select(['sample', 'distance_to_coast_km','council_olg_code'])
    # define location codes
    MET_CODES = [1, 2, 3]
    URBAN_CODES = [4, 5]
    FRINGE_CODES = [6, 7]
    RURAL_CODES = [8, 9, 10, 11]
    # define a new column in metadata_df called 'location_type' based on the 'council_olg_code' column
    metadata_df = metadata_df.with_columns([
        pl.when(pl.col('council_olg_code').is_in(MET_CODES))
        .then(pl.lit('Metropolitan'))
        .when(pl.col('council_olg_code').is_in(URBAN_CODES))
        .then(pl.lit('Urban'))
        .when(pl.col('council_olg_code').is_in(FRINGE_CODES))
        .then(pl.lit('Met_Fringe'))
        .when(pl.col('council_olg_code').is_in(RURAL_CODES))
        .then(pl.lit('Rural'))
        .otherwise(pl.lit('Unknown')).alias('location_type')
    ])
    # merge abundance_df, virus_df, and metadata_df on 'sampleid'
    merged_df = metadata_df.join(abundance_df, on='sample', how='left', validate='1:1').join(
        virus_df, left_on='sample', right_on='sampleid', how='inner', validate='1:1'
    )
    return merged_df

def run_linear_regression(merged_df):
    # run linear regression of every column apart from sample and batch with distance_to_coast_km and save to table
    results = []
    for col in merged_df.columns:
        if col not in ['sample', 'batch', 'sampleid', 'distance_to_coast_km', 'council_olg_code', 'location_type']:
            print(f"Running linear regression for {col} vs distance_to_coast_km")
            y = merged_df.select('distance_to_coast_km').to_numpy().flatten()
            x = merged_df.select(col).to_numpy().flatten()
            # if all values in x are the same, skip this column
            if len(set(x)) == 1:
                print(f"Skipping {col} because all values are the same")
                continue
            fit = fit_linear_regression(x, y)
            results.append({
                'variable': col,
                'slope': fit.slope,
                'intercept': fit.intercept,
                'r_value': fit.rvalue,
                'p_value': fit.pvalue,
                'std_err': fit.stderr
            })
    results_df = pl.DataFrame(results)
    return results_df.sort(by='p_value')

def plot_linear_regression(merged_df, results_df, outfile):
    # Make a multi-panel plot of the linear regression results for each variable
    # remove the species_count_in_trap variable from results_df
    results_df = results_df.filter(pl.col('variable') != 'species_count_in_trap')
    ncols = 4
    nrows = (len(results_df) + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 5 * nrows))
    axs = axs.flatten()
    # Plot the linear regression results
    for index, row in enumerate(results_df.iter_rows()):

        ax = axs[index]
        variable = row[0]
        r_value = row[3]
        p_value = row[4]

        y = merged_df.select('distance_to_coast_km').to_numpy().flatten()
        x = merged_df.select(variable).to_numpy().flatten()

        fit = fit_linear_regression(x, y)
        line_x, mean, lower, upper = fit.predict()
        ax.plot(line_x, mean, color='grey', linewidth=2)
        ax.fill_between(line_x, lower, upper, color='grey', alpha=0.2)
        if variable == 'species_count_in_trap':
            sns.scatterplot(
                data=merged_df.select(variable, 'distance_to_coast_km', 'location_type').to_pandas(),
                x=variable,
                y='distance_to_coast_km',
                hue='location_type',
                alpha=0.5,
                ax=ax,
            )
        else:
            ax.scatter(x, y, color='black', alpha=0.5)
        ax.set_xlabel(variable)
        ax.set_ylabel('Distance to Coast (km)')
        ax.set_title(f'Distance to Coast vs {variable}')
        if 0 < p_value < float('inf'):
            mantissa, exponent = f'{p_value:.3e}'.split('e')
            p_label = rf'{mantissa} \times 10^{{{int(exponent)}}}'
        else:
            p_label = str(p_value)
        ax.text(0.05, 0.95, rf'R² = {r_value**2:.3f}, $p = {p_label}$', transform=ax.transAxes, verticalalignment='top')
    for ax in axs[len(results_df):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

def species_count_analysis(merged_df, outfile):
    # run linear regression with species count as outcome and distance to coast as predictor
    plot_df = merged_df.filter(
        pl.col('distance_to_coast_km').is_finite()
        & pl.col('species_count_in_trap').is_finite()
    ).with_columns(
        pl.when(pl.col('distance_to_coast_km') <= 20)
        .then(pl.lit('Coastal'))
        .otherwise(pl.lit('Inland'))
        .alias('coastal'),
        pl.col('location_type').replace({'Metropolitan': 'Met', 'Met_Fringe': 'Met-fringe'}),
    ).to_pandas()
    x = plot_df['distance_to_coast_km'].to_numpy()
    y = plot_df['species_count_in_trap'].to_numpy()
    fit = fit_linear_regression(x, y)
    location_order = ['Met', 'Urban', 'Met-fringe', 'Rural']
    location_palette = dict(zip(location_order, ['#0072B2', '#E69F00', '#CC79A7', '#009E73']))
    fig, (reg_ax, box_ax) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    sns.scatterplot(data=plot_df, x='distance_to_coast_km', y='species_count_in_trap',
                    hue='location_type', hue_order=location_order,
                    palette=location_palette, alpha=0.6, ax=reg_ax)
    line_x, mean, lower, upper = fit.predict()
    reg_ax.plot(line_x, mean, color='grey', linewidth=2)
    reg_ax.fill_between(line_x, lower, upper, color='grey', alpha=0.2)
    reg_ax.set_xlabel('Distance to coast (km)')
    reg_ax.set_ylabel('Species count in trap')
    reg_ax.set_title('Species count vs distance to coast')
    reg_ax.legend(title='Location type')
    if 0 < fit.pvalue < float('inf'):
        mantissa, exponent = f'{fit.pvalue:.3e}'.split('e')
        p_label = rf'{mantissa} \times 10^{{{int(exponent)}}}'
    else:
        p_label = str(fit.pvalue)
    reg_ax.text(0.05, 0.95, rf'R² = {fit.rvalue**2:.3f}, $p = {p_label}$',
                transform=reg_ax.transAxes, verticalalignment='top')
    box_order = [f'{coast} - {location}'
                 for coast in ['Inland', 'Coastal']
                 for location in ['Rural', 'Urban', 'Met-fringe', 'Met']]
    plot_df['box_group'] = plot_df['coastal'] + ' - ' + plot_df['location_type']
    sns.boxplot(data=plot_df, x='box_group', y='species_count_in_trap',
                order=box_order, hue='location_type', dodge=False, width=0.6,
                palette=location_palette, saturation=1, legend=False,
                showfliers=False, ax=box_ax)
    sns.stripplot(data=plot_df, x='box_group', y='species_count_in_trap',
                  order=box_order, color='black', jitter=0.18, size=3,
                  alpha=0.5, zorder=3, ax=box_ax)
    box_ax.set_xticks(range(8), ['Rural', 'Urban', 'Met-fringe', 'Met'] * 2,
                     rotation=45, ha='right')
    box_ax.axvline(3.5, color='grey', linewidth=0.8, linestyle='--')
    for position, label in [(1.5, 'Inland (>20 km)'), (5.5, 'Coastal (≤20 km)')]:
        box_ax.text(position, -0.25, label, transform=box_ax.get_xaxis_transform(),
                    ha='center', va='top')
    box_ax.set_xlabel('')
    box_ax.set_ylabel('Species count in trap')
    box_ax.set_title('Species count by location and coast')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    # save the fit parameters to a tsv file
    fit_params = {
        'slope': fit.slope,
        'intercept': fit.intercept,
        'r_value': fit.rvalue,
        'p_value': fit.pvalue,
        'std_err': fit.stderr
    }
    fit_params_df = pl.DataFrame([fit_params])
    fit_params_df.write_csv(f"{RESULTS_DIR}/abundance-geospatial/species_count_linear_regression_results.tsv", separator='\t')

def prepare_community_analysis_data():
    # desired - a combined season+week variable (temporal signal), distance from coast in km, type of landscape (urban/rural/fringe); proportion of each major mosquito species; prop Cx-AUS (or split culex lineages); the nc2 coverage of each virus, the log10(dedup reads) for each virus
    abundance_df, virus_df, metadata_df = parse_data()
    abundance_df = add_species_count(abundance_df)
    diversity_df = pl.read_csv(
        f'{RESULTS_DIR}/abundance-geospatial/tables/simpsons_diversity_index.tsv',
        separator='\t',
    ).select('sample', pl.col('Simpsons_Diversity_Index').alias('simpson_diversity_index'))
    # define location codes
    MET_CODES = [1, 2, 3]
    URBAN_CODES = [4, 5]
    FRINGE_CODES = [6, 7]
    RURAL_CODES = [8, 9, 10, 11]
    # define a new column in metadata_df called 'location_type' based on the 'council_olg_code' column
    metadata_df = metadata_df.with_columns([
        pl.when(pl.col('council_olg_code').is_in(MET_CODES))
        .then(pl.lit('Metropolitan'))
        .when(pl.col('council_olg_code').is_in(URBAN_CODES))
        .then(pl.lit('Urban'))
        .when(pl.col('council_olg_code').is_in(FRINGE_CODES))
        .then(pl.lit('Met_Fringe'))
        .when(pl.col('council_olg_code').is_in(RURAL_CODES))
        .then(pl.lit('Rural'))
        .otherwise(pl.lit('Unknown')).alias('location_type')
    ])
    # Keep every metadata sample; join the selected NC2 columns once below.
    merged_df = metadata_df.join(abundance_df, on='sample', how='left', validate='1:1')
    # make a combined season+week variable
    merged_df = merged_df.with_columns([
        (pl.col('season').cast(pl.Utf8) + '_' + pl.col('weeknum').cast(pl.Utf8)).alias('season_week')
    ])
    cx_annu_props = f"{RESULTS_DIR}/publication_tables/Cx_annu_lineage_informative_snp_summary.csv"
    read_counts = f"{RESULTS_DIR}/final_castanet_data/COSMOS1-4.final_castanet_reads_dedup.tsv"
    cx_annu_props_df = pl.read_csv(cx_annu_props)
    read_counts_df = pl.read_csv(read_counts, separator='\t')
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
    # filter virus df and read_counts_df to only include columns in virus_acronyms
    virus_df = virus_df.select(['sampleid'] + [col for col in virus_df.columns if col in virus_acronyms.keys()])
    read_counts_df = read_counts_df.select(['sampleid'] + [col for col in read_counts_df.columns if col in virus_acronyms.keys()])
    # rename virus_df and read_counts_df columns to acronyms
    virus_df = virus_df.rename({col: virus_acronyms[col] for col in virus_df.columns if col in virus_acronyms.keys()})
    read_counts_df = read_counts_df.rename({col: virus_acronyms[col] for col in read_counts_df.columns if col in virus_acronyms.keys()})
    # add "_nc2" to the virus_df columns and "_dedup_reads" to the read_counts_df columns
    virus_df = virus_df.rename({col: col + '_nc2' for col in virus_df.columns if col not in ['sampleid']})
    read_counts_df = read_counts_df.rename({col: col + '_dedup_reads' for col in read_counts_df.columns if col not in ['sampleid']})
    # Log10-transform measured read counts; nulls remain null and zeros become -inf.
    read_counts_df = read_counts_df.with_columns([
        pl.col(col).log10().alias(col) for col in read_counts_df.columns if col not in ['sampleid']
    ])
    # Label the read-count columns with the transform used.
    read_counts_df = read_counts_df.rename({col: col + '_log10' for col in read_counts_df.columns if col not in ['sampleid']})
    # merge merged_df with cx_annu_props_df, virus_df, and read_counts_df on 'sample' and 'sampleid'
    merged_df = (
        merged_df.join(cx_annu_props_df, on='sample', how='left', validate='1:1')
        .join(diversity_df, on='sample', how='left', validate='1:1')
        .join(virus_df, left_on='sample', right_on='sampleid', how='left', validate='1:1')
        .join(read_counts_df, left_on='sample', right_on='sampleid', how='left', validate='1:1')
    )
    assert merged_df.shape[0] == metadata_df.shape[0], "Merged dataframe has different number of rows than metadata dataframe"

    # select only the columns we want to keep
    mosquito_cols = [col for col in abundance_df.columns
                     if col not in ['sample', 'batch', 'species_count_in_trap']]
    virus_order = ['SINDV', 'STRV', 'RRV', 'KOK', 'BFV', 'EHV', 'KUN', 'MVEV', 'JEV', 'LSV']
    cols_to_keep = (
        ['sample', 'season_week', 'distance_to_coast_km', 'location_type', 'council_olg_code', 'species_count_in_trap',
         'simpson_diversity_index',
         'location-site', 'location', 'longitude', 'latitude']
        + mosquito_cols
        + ['mean_S-AUS_subtype_ratio']
        + [f'{virus}_nc2' for virus in virus_order if f'{virus}_nc2' in virus_df.columns]
        + [f'{virus}_dedup_reads_log10' for virus in virus_order
           if f'{virus}_dedup_reads_log10' in read_counts_df.columns]
    )
    merged_df = merged_df.select(cols_to_keep).rename({
        col: 'ab_' + col.replace(' ', '_') for col in mosquito_cols
    })
    merged_df = merged_df.with_columns(
        (pl.col('ab_Culex_annulirostris') * pl.col('mean_S-AUS_subtype_ratio'))
        .alias('ab_S-AUS_Cx_annu')
    ).with_columns(
        (pl.col('ab_Culex_annulirostris') - pl.col('ab_S-AUS_Cx_annu'))
        .alias('ab_AUS_Cx_annu')
    )

    # add count of viruses in each trap >= 0.05
    virus_cols = [f'{virus}_nc2' for virus in virus_order if f'{virus}_nc2' in virus_df.columns]
    merged_df = merged_df.with_columns(
        pl.sum_horizontal([pl.col(col) >= 0.05 for col in virus_cols]).alias('virus_count_in_trap')
    )

    # add a column for any virus detected as a binary
    merged_df = merged_df.with_columns(
        pl.when(pl.any_horizontal([pl.col(col) >= 0.05 for col in virus_cols]))
        .then(pl.lit(1)).otherwise(pl.lit(0)).alias('any_virus_detected')
    )

    return merged_df

def logistic_regression_analysis(merged_df):
    # run a logistic of all mosquito species against any virus
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
    # get all the values
    virus_cols = [virus_acronyms[virus] for virus in virus_acronyms]
    merged_df = merged_df.with_columns(pl.any_horizontal([pl.col(col) > 0 for col in virus_cols]).alias('any_virus'))
    # Fit each mosquito abundance separately; virus indicators define the outcome only.
    excluded_cols = {
        'sample', 'batch', 'sampleid', 'distance_to_coast_km', 'council_olg_code',
        'location_type', 'any_virus', 'species_count_in_trap',
    } | set(virus_cols)
    results = []
    for col in merged_df.columns:
        if col not in excluded_cols:
            print(f"Running logistic regression for {col} vs any_virus")
            y = merged_df.select('any_virus').to_numpy().flatten()
            x = merged_df.select(col).to_numpy().flatten()
            try:
                fit = fit_logistic_regression(x, y)
            except ValueError as error:
                print(f"Skipping {col}: {error}")
                continue
            results.append({
                'variable': col,
                'coef': fit.slope,
                'intercept': fit.intercept,
                'p_value': fit.pvalue,
                'odds_ratio': fit.odds_ratio,
                'n': fit.n,
            })
    results_df = pl.DataFrame(results, schema={
        'variable': pl.String, 'coef': pl.Float64, 'intercept': pl.Float64,
        'p_value': pl.Float64, 'odds_ratio': pl.Float64, 'n': pl.Int64,
    })
    return results_df.sort(by='p_value')



def main():
    merged_df = combine_data()
    # save merged_df to csv
    merged_df.write_csv(f"{RESULTS_DIR}/abundance-geospatial/merged_data.csv", separator=',')
    results_df = run_linear_regression(merged_df)
    # save results_df to csv
    results_df.write_csv(f"{RESULTS_DIR}/abundance-geospatial/linear_regression_results.tsv", separator='\t')
    # plot linear regression results
    plot_linear_regression(merged_df, results_df, f"{RESULTS_DIR}/abundance-geospatial/linear_regression_plots.png")
    species_count_analysis(merged_df, f"{RESULTS_DIR}/abundance-geospatial/species_count_linear_regression_plot.png")
    # run logistic regression analysis
    logistic_results_df = logistic_regression_analysis(merged_df)
    # save logistic_results_df to csv
    logistic_results_df.write_csv(f"{RESULTS_DIR}/abundance-geospatial/logistic_regression_results.tsv", separator='\t')



if __name__ == "__main__":
    main()
