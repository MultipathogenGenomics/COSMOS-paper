from analysis.common.paths import RESULTS_DIR
from analysis.common.regression import fit_linear_regression, fit_logistic_regression
from analysis.common.mosquito_ids import shorted_genus
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy.stats import kruskal
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova as skbio_permanova
from analysis.mosquito_community_composition.scripts.process import centre_log_ratio_transform
from analysis.mosquito_community_composition.scripts.linregress_abundance_viruses_coast import prepare_community_analysis_data
import statsmodels.formula.api as smf
import statsmodels.api as sm

HEATMAP_CMAP = 'PuOr'  # ColourBrewer diverging orange-white-purple, centred on zero.
OUTLINE_COLOUR = 'red'  # Shared by the heatmap boxes and MVEV trap outlines.


def heatmap_label(name):
    labels = {
        'distance_to_coast_km_log10': "Distance to coast\n" + r"(km; $\log_{10}$)",
        'simpson_diversity_index': r"Simpson's diversity ($D_1$)",
        'virus_count_in_trap': 'Virus richness',
        'council_olg_code': 'Location type',
        'location_score': 'Location type (ordered)',
        'location_Urban': 'Location: urban',
        'location_Met': 'Location: metropolitan',
        'location_Met_fringe': 'Location: metropolitan fringe',
        'location_Rural': 'Location: rural',
        'season': 'Season (numeric trend)',
        'week': 'Week',
    }
    if name in labels:
        return labels[name]
    if name.startswith('season_') and name.removeprefix('season_').isdigit():
        code = name.removeprefix('season_')
        return f'Season: 20{code[:2]}-{code[2:]}'
    if name == 'ab_S-AUS_Cx_annu':
        return r'$\it{Cx.\ annulirostris}$ (S-AUS)'
    if name == 'ab_AUS_Cx_annu':
        return r'$\it{Cx.\ annulirostris}$ (AUS)'
    if name.startswith('ab_'):
        species = shorted_genus(name.removeprefix('ab_').replace('_', ' '))
        return r'$\it{' + species.replace(' ', r'\ ') + '}$'
    return name.removesuffix('_detected')


def kruskal_wallis_on_D(df, grouping_variable, value_variable):
    groups = (
        df.drop_nulls([grouping_variable, value_variable])
        .group_by(grouping_variable)
        .agg(pl.col(value_variable))
        .get_column(value_variable)
        .to_list()
    )
    result = kruskal(*groups)
    print(f"\nKruskal-Wallis test for {value_variable} by {grouping_variable}:\nH={result.statistic}, p={result.pvalue}")

def run_permanova(df, grouping_variable):
    # get out the mosquito abundance columns
    mosquito_cols = [col for col in df.columns if col.startswith('ab_')]
    # remove the Cx annulirostris column
    mosquito_cols.remove('ab_Culex_annulirostris')
    # get groups
    df = df.drop_nulls(mosquito_cols + ['sample', grouping_variable])
    ids = df.select("sample").to_numpy().flatten()
    groups = df.get_column(grouping_variable).to_list()
    # filter df
    species_df = df.select(mosquito_cols)
    # make a matrix of the mosquito abundance columns
    dist_mat = beta_diversity(metric='braycurtis', counts=species_df.to_numpy(), ids=ids)
    perm_res = skbio_permanova(dist_mat, groups, permutations=999)
    pseudo_f = perm_res['test statistic']
    n = perm_res['sample size']
    k = perm_res['number of groups']
    r_squared = pseudo_f * (k - 1) / (pseudo_f * (k - 1) + n - k)
    print(f"\nPERMANOVA test for mosquito abundance by {grouping_variable}:\nF={pseudo_f}, p={perm_res['p-value']}, R²={r_squared}")

def log_transform_desired_columns(df, columns, pseudocount=1e-10):
    for col in columns:
        df = df.with_columns(
            pl.when(pl.col(col).is_null())
            .then(pl.lit(None))
            .otherwise((pl.col(col) + pseudocount).log10())
            .alias(col + '_log10')
        )
    return df

def standardise_by_z_score(df, columns):
    for col in columns:
        if not df.schema[col].is_numeric():
            raise ValueError(f"Cannot z-score non-numeric column {col!r}")
        mean = df.select(pl.col(col)).mean().item()
        std = df.select(pl.col(col)).std().item()
        if std is None or std == 0:
            raise ValueError(f"Cannot z-score column {col!r}: insufficient values or zero variance")
        df = df.with_columns(
            pl.when(pl.col(col).is_null())
            .then(pl.lit(None))
            .otherwise((pl.col(col) - mean) / std)
            .alias(col)
        )
    return df

def corr_mosquito_abundance_and_vars(df, variable_cols, abundance_transform='clr'):
    mosquito_cols = [col for col in df.columns
                     if col.startswith('ab_') and col != 'ab_Culex_annulirostris']
    matrix = corr_all_with_all(
        df.select(mosquito_cols + variable_cols), variable_cols,
        abundance_transform=abundance_transform,
    )
    mosquito_cols = [col for col in matrix.columns if col.startswith('ab_')]
    indices = [matrix.columns.index(col) for col in mosquito_cols]
    corr_cols = [f'corr_with_{col}' for col in variable_cols if col in matrix.columns]
    drivers_df = pl.DataFrame({
        f'corr_with_{col}': matrix[col].gather(indices)
        for col in variable_cols if col in matrix.columns
    }).with_columns(pl.Series('mosquito', mosquito_cols))
    return drivers_df.with_columns(
        pl.sum_horizontal([pl.col(col).abs() for col in corr_cols]).alias('Absolute_Influence')
    ).sort('Absolute_Influence', descending=True)


def CLR_transform_df(df, cols_to_transform, pseudo_count=1e-9):
    # CLR transform the specified columns in the DataFrame
    for col in cols_to_transform:
        if not df.schema[col].is_numeric():
            raise ValueError(f"Cannot CLR transform non-numeric column {col!r}")
    df = df.with_columns(
        pl.DataFrame(
            centre_log_ratio_transform(
                df.select(cols_to_transform).to_numpy().astype(float),
                pseudo_count=pseudo_count,
            ),
            schema=cols_to_transform,
        )
    )
    return df

def get_corr(x, y):
    # Calculate Pearson correlation coefficient between two arrays
    if len(x) != len(y):
        raise ValueError("Arrays must be of the same length to compute correlation")
    if len(x) < 2:
        raise ValueError("At least two data points are required to compute correlation")
    return pl.DataFrame({'x': x, 'y': y}).select(pl.corr('x', 'y')).item()

def corr_all_with_all(df, variable_cols, abundance_transform='clr'):
    if abundance_transform not in {'clr', 'zscore'}:
        raise ValueError("abundance_transform must be 'clr' or 'zscore'")
    # assuming all mosquito cols and virus cols and defined list of variable cols
    all_cols = list(dict.fromkeys(
        [col for col in df.columns if col.startswith('ab_') or col.endswith('_detected')]
        + variable_cols
    ))
    # remove the Cx annulirostris column
    all_cols = [col for col in all_cols if col != 'ab_Culex_annulirostris']
    if not all_cols:
        raise ValueError('No columns selected for correlation')
    for col in all_cols:
        if not df.schema[col].is_numeric():
            raise ValueError(f'Correlation requires a numeric column: {col!r}')
    df = df.select(all_cols).with_columns(
        [pl.when(pl.col(col).is_finite()).then(pl.col(col)).otherwise(None).alias(col)
         for col in all_cols]
    )
    # CLR transform the mosquito abundance columns.
    mosquito_cols = [col for col in all_cols if col.startswith('ab_')]
    if mosquito_cols:
        # CLR requires a complete composition, but virus/metadata gaps need not
        # discard a sample from every other correlation.
        df = df.drop_nulls(mosquito_cols)
        if abundance_transform == 'clr':
            min_non_zero_mosquito_ab = df.select(
                pl.min_horizontal([
                    pl.col(col).filter(pl.col(col) > 0).min()
                    for col in mosquito_cols
                ])
            ).item()
            if min_non_zero_mosquito_ab is None:
                raise ValueError('CLR requires at least one positive mosquito abundance in complete samples')
            df = CLR_transform_df(df, mosquito_cols, pseudo_count=min_non_zero_mosquito_ab / 2)
    unusable_cols = [col for col in all_cols if df[col].drop_nulls().n_unique() < 2]
    if unusable_cols:
        print(f'Omitting columns with insufficient data or no variation: {unusable_cols}')
    all_cols = [col for col in all_cols if col not in unusable_cols]
    if len(all_cols) < 2:
        raise ValueError('Correlation requires at least two columns with observed variation')
    if abundance_transform == 'zscore':
        df = standardise_by_z_score(df, [col for col in mosquito_cols if col in all_cols])
    # Location scores and detection columns are never explicitly standardised.
    # Pairwise complete observations; undefined correlations remain NaN.
    return pl.DataFrame({
        col2: [df.select(pl.corr(col1, col2)).item() for col1 in all_cols]
        for col2 in all_cols
    })

def cluster_map(matrix, outfile):
    labels = [heatmap_label(col).replace('\n', ' ') for col in matrix.columns]
    values = matrix.to_numpy()
    vmin, vmax = np.nanmin(values), np.nanmax(values)
    if np.isfinite(values).all():
        from scipy.cluster.hierarchy import to_tree
        grid = sns.clustermap(values, cmap=HEATMAP_CMAP, center=0, vmin=vmin, vmax=vmax,
                              linewidths=0.35, linecolor='white', figsize=(15, 13),
                              xticklabels=labels, yticklabels=labels)
        # Outline the first major dendrogram branch in displayed order.
        root = to_tree(grid.dendrogram_row.linkage)
        first = grid.dendrogram_row.reordered_ind[0]
        branch = root.left if first in root.left.pre_order() else root.right
        members = branch.pre_order()
        row_positions = [grid.dendrogram_row.reordered_ind.index(i) for i in members]
        col_positions = [grid.dendrogram_col.reordered_ind.index(i) for i in members]
        grid.ax_heatmap.add_patch(Rectangle(
            (min(col_positions) + 0.03, min(row_positions) + 0.03),
            max(col_positions) - min(col_positions) + 0.94,
            max(row_positions) - min(row_positions) + 0.94,
            fill=False, edgecolor=OUTLINE_COLOUR, linewidth=3.5, zorder=5,
        ))
        print('Outlined clustermap cluster:', [matrix.columns[i] for i in sorted(members)])
        grid.ax_heatmap.set_xlabel('Mosquito abundance, virus detection and variables')
        grid.ax_heatmap.set_ylabel('Mosquito abundance, virus detection and variables')
        plt.setp(grid.ax_heatmap.get_xticklabels(), rotation=65, ha='right', rotation_mode='anchor')
        grid.savefig(outfile, dpi=200, bbox_inches='tight')
        grid.savefig(str(outfile).rsplit('.', 1)[0] + '.pdf', bbox_inches='tight')
        plt.close(grid.fig)
    else:
        print('Some correlations are undefined; plotting without clustering or a cluster outline.')
        fig, ax = plt.subplots(figsize=(15, 13))
        sns.heatmap(values, cmap=HEATMAP_CMAP, center=0, vmin=vmin, vmax=vmax,
                    linewidths=0.35, linecolor='white', xticklabels=labels, yticklabels=labels, ax=ax)
        fig.tight_layout()
        fig.savefig(outfile, dpi=200, bbox_inches='tight')
        plt.close(fig)


def clean_df(df):
    df = df.with_columns(
        pl.col("season_week").str.replace("_", ".").cast(pl.Float64).alias("season_week_float")
    )
    location_labels = {
        'Urban': 'Urban', 'Met': 'Metropolitan',
        'Met_fringe': 'Met_Fringe', 'Rural': 'Rural',
    }
    df = df.with_columns([
        pl.when(pl.col('location_type').is_in(list(location_labels.values())))
        .then((pl.col('location_type') == category).cast(pl.Int8))
        .otherwise(None)
        .alias(f'location_{label}')
        for label, category in location_labels.items()
    ])
    # also make a column of virus detection (1 where nc2 >= 0.05, 0 otherwise)
    virus_cols = [col for col in df.columns if col.endswith('_nc2')]
    for col in virus_cols:
        df = df.with_columns(
            pl.when(pl.col(col) >= 0.05).then(1).otherwise(0).alias(col.replace('_nc2', '_detected'))
        )
    # make a separate season and week column for plotting
    df = df.with_columns(
        pl.col("season_week").str.split("_").list.get(0).cast(pl.Int32).alias("season"),
        pl.col("season_week").str.split("_").list.get(1).cast(pl.Int16).alias("week"),
    )
    return df

def encode_location(df, encoding='olg_code'):
    """Return encoded data and the location columns to include in correlations."""
    if encoding == 'olg_code':
        return df, ['council_olg_code']
    if encoding == 'binary':
        return df, ['location_Urban', 'location_Met', 'location_Met_fringe', 'location_Rural']
    if encoding == 'ordinal':
        scores = {'Rural': 0.0, 'Urban': 0.33, 'Met_Fringe': 0.66, 'Metropolitan': 1.0}
        return df.with_columns(
            pl.col('location_type').replace_strict(scores, default=None,
                                                  return_dtype=pl.Float64).alias('location_score')
        ), ['location_score']
    raise ValueError("encoding must be 'binary', 'olg_code', or 'ordinal'")

def encode_season(df, encoding='binary'):
    """Return data and season columns; week remains a separate numeric variable."""
    if encoding == 'numeric':
        return df, ['season']
    if encoding == 'binary':
        seasons = sorted(df.get_column('season').drop_nulls().unique().to_list())
        if not seasons:
            raise ValueError('Cannot encode season: no observed seasons')
        cols = [f'season_{season}' for season in seasons]
        return df.with_columns([
            (pl.col('season') == season).cast(pl.Int8).alias(col)
            for season, col in zip(seasons, cols)
        ]), cols
    raise ValueError("encoding must be 'binary' or 'numeric'")


def three_panel_summary_heatmap(correlation_matrix, outfile, desired_mosquito_cols=None,
                                embedding=None, region_markers=None):
    mosquito_cols = [col for col in correlation_matrix.columns if col.startswith('ab_')]
    virus_cols = [col for col in correlation_matrix.columns if col.endswith('_detected')]
    variable_cols = [col for col in correlation_matrix.columns if col not in mosquito_cols + virus_cols]
    if desired_mosquito_cols is not None:
        missing_cols = [col for col in desired_mosquito_cols if col not in mosquito_cols]
        if missing_cols:
            raise ValueError(f'Requested mosquito columns missing from correlation matrix: {missing_cols}')
        mosquito_cols = list(dict.fromkeys(desired_mosquito_cols))
    if not mosquito_cols or not virus_cols or not variable_cols:
        raise ValueError('Three-panel heatmap requires mosquito, virus detection, and extra variable columns')

    def block(rows, columns):
        # Matrix rows follow the original column order, even after select().
        indices = [correlation_matrix.columns.index(row) for row in rows]
        return correlation_matrix.select(columns).to_numpy()[indices, :]

    # Sort only the relevant variables, then rebuild every panel in that order.
    distance_col = 'distance_to_coast_km_log10'
    if distance_col not in correlation_matrix.columns:
        raise ValueError(f'Column ordering requires {distance_col!r}')
    distance_correlations = block([distance_col], mosquito_cols)[0]
    distance_strength = np.where(np.isfinite(distance_correlations),
                                 distance_correlations, -np.inf)
    mosquito_cols = [mosquito_cols[i] for i in np.argsort(-distance_strength, kind='stable')]
    virus_correlations = block(virus_cols, mosquito_cols)
    virus_strength = np.max(np.where(np.isfinite(virus_correlations),
                                     np.abs(virus_correlations), -np.inf), axis=1)
    virus_cols = [virus_cols[i] for i in np.argsort(-virus_strength, kind='stable')]

    mosquito_virus_corr = block(virus_cols, mosquito_cols)
    mosquito_variable_corr = block(variable_cols, mosquito_cols)
    variable_virus_corr = block(virus_cols, variable_cols)
    values = np.concatenate([panel.ravel() for panel in
                             [mosquito_virus_corr, mosquito_variable_corr, variable_virus_corr]])
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        raise ValueError('No finite correlations available for the three-panel heatmap')

    vmin, vmax = finite_values.min(), finite_values.max()
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={
        'width_ratios': [len(mosquito_cols), len(variable_cols)],
        'height_ratios': [len(virus_cols), len(variable_cols)],
    }, sharex='col', sharey='row')
    # Mosquito columns align vertically; virus rows align horizontally.
    cmap = HEATMAP_CMAP
    sns.heatmap(mosquito_virus_corr, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                linewidths=0.35, linecolor='white',
                cbar=False, xticklabels=mosquito_cols, yticklabels=virus_cols, ax=axs[0, 0])
    sns.heatmap(mosquito_variable_corr, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                linewidths=0.35, linecolor='white',
                cbar=False, xticklabels=mosquito_cols, yticklabels=variable_cols, ax=axs[1, 0])
    sns.heatmap(variable_virus_corr, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                linewidths=0.35, linecolor='white',
                cbar=False, xticklabels=variable_cols, yticklabels=virus_cols, ax=axs[0, 1])
    axs[1, 1].axis('off')
    if embedding is not None:
        # The scatter must not share the heatmaps' categorical coordinates.
        spec = axs[1, 1].get_subplotspec()
        axs[1, 1].remove()
        embedding_grid = spec.subgridspec(
            2, 2, width_ratios=[0.12, 0.88], height_ratios=[0.12, 0.88],
            wspace=0, hspace=0,
        )
        axs[1, 1] = fig.add_subplot(embedding_grid[1, 1])
        plot_mosquito_embedding(axs[1, 1], embedding, region_markers=region_markers)
    for ax in (axs[0, 0], axs[0, 1]):
        ax.xaxis.tick_top()
        ax.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False)
        plt.setp(ax.get_xticklabels(), rotation=55, ha='left', rotation_mode='anchor')
    axs[1, 0].tick_params(axis='x', top=False, labeltop=False, bottom=False, labelbottom=False)
    axs[0, 1].tick_params(axis='y', left=False, labelleft=False)
    # add cell annotations if absolute value is > threshold
    threshold = 0.1

    def select_text_color(value, mesh):
        # Use Seaborn's rendered colormap and normalisation (including center=0).
        rgb = np.asarray(mesh.to_rgba(value)[:3])
        linear_rgb = np.where(rgb <= 0.04045, rgb / 12.92,
                              ((rgb + 0.055) / 1.055) ** 2.4)
        luminance = linear_rgb @ np.array([0.2126, 0.7152, 0.0722])
        black_contrast = (luminance + 0.05) / 0.05
        white_contrast = 1.05 / (luminance + 0.05)
        return 'white' if white_contrast > black_contrast else 'black'

    for ax, values in [(axs[0, 0], mosquito_virus_corr),
                       (axs[1, 0], mosquito_variable_corr),
                       (axs[0, 1], variable_virus_corr)]:
        mesh = ax.collections[0]
        for i, j in np.ndindex(values.shape):
            value = values[i, j]
            if np.isfinite(value) and abs(value) > threshold:
                ax.text(j + 0.5, i + 0.5, f'{value:.2f}', ha='center', va='center',
                        color=select_text_color(value, mesh))

    # Shared label formatter keeps this figure and the supplement consistent.
    axs[0, 0].set_yticklabels([heatmap_label(col) for col in virus_cols])
    axs[1, 0].set_yticklabels([heatmap_label(col) for col in variable_cols])
    axs[0, 1].set_xticklabels([heatmap_label(col) for col in variable_cols],
                              rotation=55, ha='left', rotation_mode='anchor')
    axs[0, 0].set_xticklabels([heatmap_label(col) for col in mosquito_cols],
                              rotation=55, ha='left', rotation_mode='anchor')
    axs[0, 0].set_xlabel('Mosquito abundance', labelpad=12)
    axs[0, 0].xaxis.set_label_position('top')
    axs[0, 0].set_ylabel('Virus detection', labelpad=24)
    axs[1, 0].set_ylabel('Environmental and temporal\nvariables', labelpad=24)
    axs[0, 1].set_xlabel('Environmental and temporal\nvariables', labelpad=12)
    axs[0, 1].xaxis.set_label_position('top')

    # Highlight the three associations discussed in the main figure.
    for ax, row, rows in [
        (axs[0, 0], 'MVEV_detected', virus_cols),
        (axs[1, 0], 'distance_to_coast_km_log10', variable_cols),
        (axs[1, 0], 'location_Rural', variable_cols),
    ]:
        if row in rows and 'ab_S-AUS_Cx_annu' in mosquito_cols:
            ax.add_patch(Rectangle(
                (mosquito_cols.index('ab_S-AUS_Cx_annu') + 0.03, rows.index(row) + 0.03),
                0.94, 0.94, fill=False, edgecolor=OUTLINE_COLOUR, linewidth=3.5, zorder=5,
            ))

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.align_xlabels([axs[0, 0], axs[0, 1]])
    fig.align_ylabels([axs[0, 0], axs[1, 0]])
    for ax in (axs[0, 0], axs[0, 1]):
        ax.xaxis.label.set_verticalalignment('top')
    for ax in (axs[0, 0], axs[1, 0]):
        ax.yaxis.label.set_verticalalignment('top')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(str(outfile).rsplit('.', 1)[0] + '.pdf', bbox_inches='tight')
    plt.close(fig)

def compute_PCA_of_mosquito_abundance(df, method='pca', pseudo_count=1e-9,
                                     loadings_outfile=None):
    """Embed species abundances, using total Cx. annulirostris."""
    species = [col for col in df.columns
               if col.startswith('ab_') and col not in {'ab_S-AUS_Cx_annu', 'ab_AUS_Cx_annu'}]
    if not species:
        raise ValueError('No mosquito abundance columns found for embedding')
    complete = df.drop_nulls(species).filter(
        pl.all_horizontal([pl.col(col).is_finite() for col in species])
    )
    data = complete.select(species).to_numpy().astype(float)
    if complete.height < 3 or len(species) < 2:
        raise ValueError('Embedding requires at least three complete samples and two abundance columns')
    if (data < 0).any() or not np.isfinite(pseudo_count) or pseudo_count <= 0:
        raise ValueError('CLR requires nonnegative abundance and a positive finite pseudocount')
    if method == 'pcoa':
        from skbio.stats.ordination import pcoa
        # Bray-Curtis uses the original nonnegative abundance, without CLR.
        keep = data.sum(axis=1) > 0
        complete = complete.filter(pl.Series(keep))
        data = data[keep]
        distances = beta_diversity('braycurtis', counts=data,
                                   ids=complete['sample'].to_list())
        result = pcoa(distances, dimensions=complete.height)
        coordinates = result.samples.loc[complete['sample'].to_list(), ['PC1', 'PC2']].to_numpy()
        labels = [f'PCoA{i + 1} ({result.proportion_explained.iloc[i]:.1%})' for i in range(2)]
    elif method == 'pca':
        transformed = centre_log_ratio_transform(data, pseudo_count=pseudo_count)
        from sklearn.decomposition import PCA
        model = PCA(n_components=2)
        coordinates = model.fit_transform(transformed)
        labels = [f'PC{i + 1} ({ratio:.1%})' for i, ratio in enumerate(model.explained_variance_ratio_)]
        if loadings_outfile is not None:
            loadings = pl.DataFrame({
                'species': species, 'PC1': model.components_[0], 'PC2': model.components_[1],
            })
            plot_pca_loadings(loadings, labels, loadings_outfile)
    elif method == 'umap':
        transformed = centre_log_ratio_transform(data, pseudo_count=pseudo_count)
        from umap import UMAP
        coordinates = UMAP(n_components=2, random_state=42,
                           n_neighbors=min(15, complete.height - 1)).fit_transform(transformed)
        labels = ['UMAP1', 'UMAP2']
    else:
        raise ValueError("Embedding method must be 'pca', 'umap', or 'pcoa'")
    preparation = 'Bray-Curtis on original abundance' if method == 'pcoa' else f'CLR pseudocount={pseudo_count:g}'
    print(f'{method.upper()} embedding: {complete.height}/{df.height} samples, {len(species)} abundance columns (S-AUS/AUS split excluded); {preparation}')
    return complete.with_columns(
        pl.Series('embedding_x', coordinates[:, 0]),
        pl.Series('embedding_y', coordinates[:, 1]),
    ), labels


def plot_pca_loadings(loadings, component_labels, outfile):
    """Plot signed PCA component coefficients for every abundance column."""
    from pathlib import Path

    loadings = loadings.sort('PC1', descending=True)

    def species_label(name):
        lineage_labels = {
            'ab_S-AUS_Cx_annu': 'Cx. annulirostris (S-AUS)',
            'ab_AUS_Cx_annu': 'Cx. annulirostris (AUS)',
        }
        if name in lineage_labels:
            return lineage_labels[name]
        return shorted_genus(name.removeprefix('ab_').replace('_', ' '))

    names = [species_label(name) for name in loadings['species']]
    fig, axes = plt.subplots(1, 2, figsize=(12, max(6, 0.30 * len(names))), sharey=True)
    positions = np.arange(len(names))
    for ax, component, title in zip(axes, ['PC1', 'PC2'], component_labels):
        values = loadings[component].to_numpy()
        ax.barh(positions, values, color=np.where(values >= 0, '#F28E2B', '#0072B2'))
        ax.axvline(0, color='0.3', linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel('PCA coefficient (signed loading)')
        ax.set_yticks(positions, labels=names, fontstyle='italic')
        sns.despine(ax=ax)
    axes[0].invert_yaxis()
    limit = np.max(np.abs(loadings.select('PC1', 'PC2').to_numpy())) * 1.1
    for ax in axes:
        ax.set_xlim(-limit, limit)
    fig.tight_layout()
    out = Path(outfile)
    fig.savefig(out, dpi=200, bbox_inches='tight')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    loadings.write_csv(out.with_suffix('.tsv'), separator='\t')
    plt.close(fig)
    print(f'Saved PCA loadings: {out}')


def plot_mosquito_embedding(ax, embedding, region_markers=None):
    from matplotlib.lines import Line2D
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    region_markers = region_markers or {
        'Rural': 'o', 'Urban': 's', 'Met_Fringe': '^', 'Metropolitan': 'D',
    }
    data, labels = embedding
    # Same blue (AUS) -> orange (S-AUS) gradient as results_figure.py.
    cmap = sns.blend_palette(['#0072B2', '#F28E2B'], as_cmap=True)
    norm = Normalize(0, 1)
    saus = data['mean_S-AUS_subtype_ratio'].to_numpy().astype(float)
    annu = data['ab_Culex_annulirostris'].to_numpy().astype(float)
    has_lineage_colour = np.isfinite(annu) & (annu > 0) & np.isfinite(saus)
    colors = cmap(norm(np.clip(saus, 0, 1)))
    colors[~has_lineage_colour, :3] = 0.7
    colors[:, 3] = 1.0
    detected = data['MVEV_detected'].fill_null(0).to_numpy().astype(bool)
    regions = data['location_type'].fill_null('Unknown').to_numpy()
    x, y = data['embedding_x'].to_numpy(), data['embedding_y'].to_numpy()
    handles = []
    for region in sorted(set(regions)):
        marker = region_markers.get(region, 'X')
        mask = regions == region
        ax.scatter(x[mask], y[mask], c=colors[mask], marker=marker, s=30,
                   edgecolors='#777777', linewidths=0.35)
        positive = mask & detected
        if positive.any():
            # Draw MVEV outlines above the filled points.
            ax.scatter(x[positive], y[positive], facecolors='none', edgecolors=OUTLINE_COLOUR,
                       marker=marker, s=30, linewidths=2, zorder=4)
        handles.append(Line2D([], [], marker=marker, linestyle='none', color='0.4',
                              markersize=5, label=region.replace('_', ' ')))
    handles.append(Line2D([], [], marker='o', linestyle='none', markerfacecolor='none',
                          markeredgecolor=OUTLINE_COLOUR, markeredgewidth=2, label='MVEV detected'))
    ax.set_box_aspect(1)
    ax.set_anchor('W')
    ax.legend(handles=handles, fontsize=6, loc='upper left', bbox_to_anchor=(1.08, 1),
              frameon=False, borderaxespad=0)
    ax.set_xlabel(labels[0], fontsize=8)
    ax.set_ylabel(labels[1], fontsize=8)
    ax.tick_params(labelsize=7)
    cax = ax.inset_axes([1.12, 0.08, 0.045, 0.46])
    bar = ax.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                            orientation='vertical')
    bar.outline.set_visible(False)
    bar.set_label('S-AUS proportion within\nCx. annulirostris', fontsize=7)
    bar.ax.tick_params(labelsize=6)
    sns.despine(ax=ax)

def linear_regression_with_random_effect(df, outcome, continuous_predictors=None, categorical_predictors=None, group_variable=None, suppress_print=False):
    if continuous_predictors is None and categorical_predictors is None:
        raise ValueError("At least one of continuous_predictors or categorical_predictors must be provided")
    if continuous_predictors is None:
        continuous_predictors = []
    if categorical_predictors is None:
        categorical_predictors = []
    # deal with missing values by dropping rows with NaN in any of the relevant columns
    relevant_columns = [outcome] + continuous_predictors + categorical_predictors + [group_variable]
    df = df.dropna(subset=relevant_columns)
    terms = continuous_predictors + [f'C({col})' for col in categorical_predictors]
    formula = f"{outcome} ~ " + ' + '.join(terms)
    if not suppress_print:
        print(f"Fitting mixed linear model with formula: {formula} and random effect for {group_variable}")
    import warnings
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter('always')
        mixed_model = smf.mixedlm(formula=formula, data=df, groups=df[group_variable]).fit()
    fit_notes = '; '.join(dict.fromkeys(str(item.message) for item in fit_warnings))
    if not suppress_print:
        print(mixed_model.summary())
        if fit_notes:
            print(f'Fit notes: {fit_notes}')
    # Fixed effects only: Group Var is not a predictor coefficient.
    fixed_terms = mixed_model.fe_params.index
    intervals = mixed_model.conf_int().loc[fixed_terms]
    return {
        'formula': formula,
        'intercept': mixed_model.fe_params['Intercept'],
        'coefficients': mixed_model.fe_params.drop('Intercept').to_dict(),
        'pvalues': mixed_model.pvalues.loc[fixed_terms].drop('Intercept').to_dict(),
        'confidence_intervals': {term: tuple(intervals.loc[term]) for term in fixed_terms if term != 'Intercept'},
        'n': int(mixed_model.nobs),
        'n_sites': int(df[group_variable].nunique()),
        'converged': bool(mixed_model.converged),
        'fit_notes': fit_notes,
        'reference_locations': {
            col: next(str(category) for category in df[col].dropna().unique()
                      if f'C({col})[T.{category}]' not in fixed_terms)
            for col in categorical_predictors
        },
        'site_variance': float(mixed_model.cov_re.iloc[0, 0]),
        'residual_variance': float(mixed_model.scale),
        'r_squared': None,  # MixedLM does not supply ordinary OLS R-squared.
    }

def logistic_regression_with_random_effect(df, outcome, continuous_predictors=None,
                                           categorical_predictors=None, group_variable=None,
                                           suppress_print=False, effect_increments=None):
    """Site-random-intercept logit; report effects per requested original-scale increment."""
    import pandas as pd
    import warnings
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

    continuous_predictors = continuous_predictors or []
    categorical_predictors = categorical_predictors or []
    effect_increments = effect_increments or {}
    for col, increment in effect_increments.items():
        if col not in continuous_predictors or not np.isfinite(increment) or increment <= 0:
            raise ValueError('Effect increments must be positive and refer to continuous predictors')
    if not continuous_predictors and not categorical_predictors:
        raise ValueError('At least one predictor must be provided')
    if group_variable is None:
        raise ValueError('A grouping column is required for the random intercept')
    relevant = [outcome] + continuous_predictors + categorical_predictors + [group_variable]
    data = df.dropna(subset=relevant).copy()
    numeric = [outcome] + continuous_predictors
    data = data.loc[np.isfinite(data[numeric].to_numpy(dtype=float)).all(axis=1)].copy()
    if set(data[outcome].unique()) != {0, 1}:
        raise ValueError('Logistic regression requires both 0 and 1 outcomes')
    if data[group_variable].nunique() < 2:
        raise ValueError('Random-intercept model requires at least two sites')
    data[outcome] = data[outcome].astype(int)
    scales = {}
    for col in continuous_predictors:
        mean, sd = float(data[col].mean()), float(data[col].std())
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f'Cannot standardise constant predictor {col!r}')
        data[col] = (data[col] - mean) / sd
        scales[col] = {'mean': mean, 'sd': sd}
    terms = continuous_predictors + [f'C({col})' for col in categorical_predictors]
    formula = f'{outcome} ~ ' + ' + '.join(terms)
    # Gaussian N(0, 2^2) fixed-effect priors and N(0, 0.5^2) prior on log site SD.
    model = BinomialBayesMixedGLM.from_formula(
        formula, {'site': f'0 + C({group_variable})'}, data, fe_p=2, vcp_p=0.5,
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        fit = model.fit_vb(mean=np.zeros(model.k_fep + model.k_vcp + model.k_vc),
                           sd=np.ones(model.k_fep + model.k_vcp + model.k_vc),
                           minim_opts={'maxiter': 1000})
    converged = bool(fit.optim_retvals['success'])
    notes = '; '.join(dict.fromkeys(str(item.message) for item in captured))
    if not converged:
        notes = (notes + '; ' if notes else '') + str(fit.optim_retvals['message'])
    # Fit on standardised predictors to retain the same prior; convert the
    # posterior back to the requested original-scale change for reporting.
    factors = np.array([
        effect_increments[term] / scales[term]['sd'] if term in effect_increments else 1.0
        for term in model.exog_names
    ])
    means, sds = fit.fe_mean * factors, fit.fe_sd * factors
    lower, upper = means - 1.96 * sds, means + 1.96 * sds
    table = pd.DataFrame({
        'term': model.exog_names,
        'effect_unit': [(f'per +{effect_increments[term]:g} week' if term == 'week' else
                         f'per +{effect_increments[term]:g} subtype ratio' if 'subtype_ratio' in term else
                         f'per +{effect_increments[term]:g} abundance proportion') if term in effect_increments else
                        ('per 1 SD' if term in continuous_predictors else
                        ('intercept' if term == 'Intercept' else 'vs reference category')
                        ) for term in model.exog_names],
        'log_odds_posterior_mean': means,
        'log_odds_posterior_sd': sds,
        'odds_ratio': np.exp(means),
        'or_credible_lower': np.exp(lower), 'or_credible_upper': np.exp(upper),
    })
    if not suppress_print:
        print(f'Bayesian logistic mixed model: {formula}; site random intercept: {group_variable}')
        print(f'N={len(data)}, detections={int(data[outcome].sum())}, sites={data[group_variable].nunique()}, converged={converged}')
        print('Predictors are standardised internally; reported odds ratios use effect_unit. Intervals are approximate 95% variational posterior credible intervals; no frequentist p-values.')
        print(table.to_string(index=False, float_format=lambda value: f'{value:.4g}'))
        if notes:
            print(f'Fit notes: {notes}')
    return {
        'method': 'Bayesian logistic random-intercept model (variational Bayes)',
        'formula': formula, 'table': table, 'predictor_scales': scales,
        'n': len(data), 'n_events': int(data[outcome].sum()),
        'n_sites': int(data[group_variable].nunique()),
        'converged': converged, 'fit_notes': notes,
        'site_log_sd_posterior_mean': float(fit.vcp_mean[0]),
        'fe_prior_sd': 2.0, 'log_site_sd_prior_sd': 0.5,
    }


def save_regression_tables(results_for_saving, outdir=f'{RESULTS_DIR}/abundance-geospatial/tables'):
    """Save unrounded fixed effects and one readable summary row per outcome."""
    import pandas as pd
    from pathlib import Path

    def outcome_label(name):
        if name == 'simpson_diversity_index':
            return "Simpson's diversity"
        if name == 'species_count_in_trap':
            return 'Mosquito species richness'
        if name == 'ab_S_AUS_Cx_annu':
            return 'Cx. annulirostris (S-AUS)'
        if name == 'ab_AUS_Cx_annu':
            return 'Cx. annulirostris (AUS)'
        return shorted_genus(name.removeprefix('ab_').replace('_', ' '))

    def predictor_label(term, reference):
        if term == 'distance_to_coast_km_log10':
            return 'Log10 distance to coast'
        if term.startswith('C(location_type)[T.'):
            category = term.removeprefix('C(location_type)[T.').removesuffix(']')
            return f'{category.replace("_", " ")} vs {reference.replace("_", " ")}'
        return term

    rows, summary = [], []
    for outcome, result in results_for_saving:
        reference = result['reference_locations'].get('location_type', '')
        readable = {'Outcome': outcome_label(outcome), 'N': result['n'],
                    'Sites': result['n_sites'], 'Reference location': reference.replace('_', ' '),
                    'Converged': result['converged'], 'Fit notes': result['fit_notes']}
        terms = sorted(result['coefficients'], key=lambda term: term != 'distance_to_coast_km_log10')
        for term in terms:
            estimate = result['coefficients'][term]
            low, high = result['confidence_intervals'][term]
            pvalue = result['pvalues'][term]
            label = predictor_label(term, reference)
            rows.append({
                'outcome': outcome, 'outcome_label': outcome_label(outcome),
                'term': term, 'predictor': label, 'estimate': estimate,
                'ci_lower': low, 'ci_upper': high, 'p_value': pvalue,
                'n': result['n'], 'n_sites': result['n_sites'],
                'reference_location': reference, 'converged': result['converged'],
                'fit_notes': result['fit_notes'],
                'site_variance': result['site_variance'],
                'residual_variance': result['residual_variance'], 'formula': result['formula'],
            })
            readable[f'{label}: estimate (95% CI)'] = f'{estimate:.3g} ({low:.3g}, {high:.3g})'
            readable[f'{label}: p'] = f'{pvalue:.3g}' if np.isfinite(pvalue) else 'NA'
        summary.append(readable)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    numeric = pd.DataFrame(rows)
    formatted = pd.DataFrame(summary)
    numeric.to_csv(outdir / 'mixed_model_coefficients.tsv', sep='\t', index=False)
    formatted.to_csv(outdir / 'mixed_model_summary.tsv', sep='\t', index=False)
    note = (
        'Separate linear mixed models with a random intercept for trap site. '
        'Fixed effects: log10 distance to coast and location type. '
        'Estimates are on the original outcome scale; abundance is a proportion, not CLR or z-scored. '
        'Distance estimates describe a one-unit increase in log10(distance + 1e-10). '
        'Location estimates compare with the reference listed in the table, adjusting for distance. '
        'Confidence intervals are 95% Wald intervals; p-values are unadjusted. '
        'Interpret nonconverged fits cautiously; Group Var is excluded from predictor coefficients.'
    )
    html = '<!doctype html><html><head><meta charset="utf-8"><title>Mixed-model results</title>'
    html += '<style>body{font:14px system-ui;margin:24px}table{border-collapse:collapse}th,td{padding:8px;border:1px solid #ddd;text-align:left}th{background:#eee}p{max-width:1100px}</style></head><body>'
    html += '<h1>Mixed-model regression results</h1><p>' + note + '</p>'
    html += formatted.to_html(index=False, border=0, na_rep='NA') + '</body></html>'
    (outdir / 'mixed_model_summary.html').write_text(html)
    print(f'Saved regression tables for {len(summary)} outcomes to {outdir}')
    return numeric, formatted


def main():
    # TODO: check simpsons diversity index

    # Analysis settings
    LOCATION_ENCODING = 'binary'  # 'binary', 'olg_code', or 'ordinal'
    SEASON_ENCODING = 'binary'  # 'binary' indicators or 'numeric' chronological trend
    ABUNDANCE_TRANSFORM = 'zscore'  # 'clr' or 'zscore'
    EMBEDDING_METHOD = 'pcoa'  # Bray-Curtis PCoA; 'pca'/'umap' use CLR
    EMBEDDING_PSEUDOCOUNT = 1e-9
    REGION_MARKERS = {'Rural': 'o', 'Urban': 's', 'Met_Fringe': '^', 'Metropolitan': 'D'}
    HEATMAP_OUTFILE = None  # Optional path; None names the file from the settings above.

    # make data
    df = prepare_community_analysis_data()

    # convert season and week to string - check this is a good idea
    df = clean_df(df) # trasnform season_week to float

    # log transfor distance to coast
    df = log_transform_desired_columns(df, ['distance_to_coast_km']) # added _log10 to the column name

    # kruskal-wallis test on D
    kruskal_wallis_on_D(df, 'location_type', 'simpson_diversity_index')

    # permanova test on mosquito abundance by location type
    run_permanova(df, 'location_type')

    df, location_cols = encode_location(df, LOCATION_ENCODING)
    df, season_cols = encode_season(df, SEASON_ENCODING)
    variable_cols = [
        'distance_to_coast_km_log10', 
        'simpson_diversity_index', 
        'virus_count_in_trap',
        'week'
        ] + season_cols + location_cols
    drivers_df = corr_mosquito_abundance_and_vars(
        df, variable_cols, abundance_transform=ABUNDANCE_TRANSFORM
    )

    desired_mosquito_cols = [
        "ab_S-AUS_Cx_annu", "ab_AUS_Cx_annu", 
        "ab_Aedes_notoscriptus",
        "ab_Culex_australicus",
        "ab_Anopheles_annulipes",
        "ab_Culex_pipiens",
        "ab_Mansonia_uniformis",
        "ab_Aedes_vigilax",
        "ab_Culex_sitiens",
        "ab_Ochlerotatus_procax",
        "ab_Coquillettidia_linealis",
    ]

    # get correlation between all mosquito abundance and virus counts and other variables
    correlation_matrix = corr_all_with_all(
        df, variable_cols, abundance_transform=ABUNDANCE_TRANSFORM
    )
    outfile = HEATMAP_OUTFILE or (
        f'{RESULTS_DIR}/abundance-geospatial/correlation_matrix_{LOCATION_ENCODING}_{ABUNDANCE_TRANSFORM}.png'
    )
    cluster_map(correlation_matrix, outfile)

    # three panel summary heatmap
    outfile = f'{RESULTS_DIR}/abundance-geospatial/three_panel_summary_heatmap_{LOCATION_ENCODING}_{ABUNDANCE_TRANSFORM}.png'
    embedding = compute_PCA_of_mosquito_abundance(
        df, method=EMBEDDING_METHOD, pseudo_count=EMBEDDING_PSEUDOCOUNT,
        loadings_outfile=f'{RESULTS_DIR}/abundance-geospatial/mosquito_pca_loadings.png',
    )
    three_panel_summary_heatmap(
        correlation_matrix, outfile, desired_mosquito_cols=desired_mosquito_cols,
        embedding=embedding, region_markers=REGION_MARKERS,
    )

    # Fit diversity, richness and abundance outcomes, then save a readable table.
    regression_df = df.rename({col: col.replace('-', '_') for col in df.columns}).to_pandas()
    results_for_saving = []
    outcomes = ['simpson_diversity_index', 'species_count_in_trap'] + [
        col.replace('-', '_') for col in desired_mosquito_cols
    ]
    for outcome in outcomes:
        results = linear_regression_with_random_effect(
            df=regression_df,
            outcome=outcome,
            continuous_predictors=['distance_to_coast_km_log10'],
            categorical_predictors=['location_type'],
            group_variable='location_site',
            suppress_print=True,
        )
        results_for_saving.append((outcome, results))
    save_regression_tables(results_for_saving)

    # Separate small models: MVEV ~ one mosquito + categorical season + (1 | site).
    import pandas as pd
    logistic_tables, logistic_summary = [], []
    predictors = {
        'ab_S_AUS_Cx_annu': 'Cx. annulirostris (S-AUS)',
        'ab_AUS_Cx_annu': 'Cx. annulirostris (AUS)',
        'ab_Culex_australicus': 'Cx. australicus',
    }
    for predictor, label in predictors.items():
        result = logistic_regression_with_random_effect(
            df=regression_df, outcome='MVEV_detected',
            continuous_predictors=[predictor], categorical_predictors=['season'],
            group_variable='location_site', suppress_print=True,
            effect_increments={predictor: 0.10},
        )
        table = result['table'].assign(
            mosquito=label, n=result['n'], detections=result['n_events'],
            sites=result['n_sites'], converged=result['converged'],
            fit_notes=result['fit_notes'], formula=result['formula'],
            predictor_mean=result['predictor_scales'][predictor]['mean'],
            predictor_sd=result['predictor_scales'][predictor]['sd'],
        )
        logistic_tables.append(table)
        effect = table.loc[table['term'] == predictor].iloc[0]
        logistic_summary.append({
            'Mosquito': label, 'N': result['n'], 'MVEV detections': result['n_events'],
            'Sites': result['n_sites'], 'Abundance increase': '10 percentage points',
            'Odds ratio': effect['odds_ratio'],
            '95% credible lower': effect['or_credible_lower'],
            '95% credible upper': effect['or_credible_upper'],
            'Converged': result['converged'], 'Fit notes': result['fit_notes'],
        })
    outdir = f'{RESULTS_DIR}/abundance-geospatial/tables'
    pd.concat(logistic_tables, ignore_index=True).to_csv(
        f'{outdir}/mvev_logistic_mixed_model.tsv', sep='\t', index=False,
    )
    summary = pd.DataFrame(logistic_summary)
    summary.to_csv(f'{outdir}/mvev_logistic_summary.tsv', sep='\t', index=False)
    summary.to_html(f'{outdir}/mvev_logistic_summary.html', index=False,
                    float_format=lambda value: f'{value:.3g}')
    print('MVEV models adjusted for categorical season with a site random intercept; '
          'OR per 10-percentage-point abundance increase; approximate Bayesian credible intervals.')
    print(summary.to_string(index=False, float_format=lambda value: f'{value:.3g}'))

    # Within-season models: MVEV ~ one mosquito + week + (1 | location-site).
    season_df = regression_df.loc[regression_df['season'] == 2223].copy()
    print('\n2022-23 only: MVEV detection ~ mosquito abundance + week + (1 | location-site)')
    print('Abundance ORs are per 10-percentage-point increase; week ORs are per one week. '
          'Intervals are approximate 95% Bayesian credible intervals, not frequentist confidence intervals.')
    within_season_predictors = {**predictors, 'ab_Culex_annulirostris': 'Cx. annulirostris (total)'}
    for predictor, label in within_season_predictors.items():
        result = logistic_regression_with_random_effect(
            df=season_df, outcome='MVEV_detected',
            continuous_predictors=[predictor, 'week'],
            group_variable='location_site', suppress_print=True,
            effect_increments={predictor: 0.10, 'week': 1.0},
        )
        effects = result['table'].set_index('term')
        print(f"\n{label}: N={result['n']}, MVEV detections={result['n_events']}, "
              f"sites={result['n_sites']}, converged={result['converged']}")
        for term, description in [(predictor, 'Abundance (+10 percentage points)'),
                                  ('week', 'Week (+1)')]:
            effect = effects.loc[term]
            print(f"  {description}: OR={effect['odds_ratio']:.3g} "
                  f"(95% credible interval {effect['or_credible_lower']:.3g}–{effect['or_credible_upper']:.3g})")
        effect = effects.loc[predictor]
        direction = 'higher' if effect['odds_ratio'] > 1 else 'lower'
        interval_note = ('includes 1' if effect['or_credible_lower'] <= 1 <= effect['or_credible_upper']
                         else 'excludes 1')
        print(f'  Adjusting for week and site, the estimated association is with {direction} odds '
              f'of MVEV detection; the approximate credible interval {interval_note}.')
        if result['fit_notes']:
            print(f"  Fit notes: {result['fit_notes']}")

    # Subtype composition is interpretable only where Cx. annulirostris is present.
    ratio_df = season_df.loc[
        (season_df['ab_Culex_annulirostris'] > 0)
        & np.isfinite(season_df['mean_S_AUS_subtype_ratio'])
    ].copy()
    print('\n2022-23: total Cx. annulirostris + continuous S-AUS subtype ratio + week '
          '+ (1 | location-site); restricted to Cx. annulirostris-positive traps with a measured ratio.')
    ratio_result = logistic_regression_with_random_effect(
        df=ratio_df, outcome='MVEV_detected',
        continuous_predictors=['ab_Culex_annulirostris', 'mean_S_AUS_subtype_ratio', 'week'],
        group_variable='location_site', suppress_print=False,
        effect_increments={'ab_Culex_annulirostris': 0.10,
                           'mean_S_AUS_subtype_ratio': 0.10, 'week': 1.0},
    )


if __name__ == "__main__":
    main()
