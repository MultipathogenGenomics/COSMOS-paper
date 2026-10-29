"""Figure 5 from frozen PCA scores, loadings and association results."""
from pathlib import Path
import argparse
import numpy as np
import polars as pl
import seaborn as sns
from matplotlib import pyplot as plt
from analysis.common.mosquito_ids import shorted_genus
from analysis.common.paths import RESULTS_DIR
from analysis.common.release_data_parsers import RELEASE_DATA_DIR, read_release_pca

def plot_virus_positives_on_PCA_with_colours(df, pca_scores, fig, ax, colour_dict=None):
    # plot PCA scores with points coloured by virus presence/absence
    if colour_dict is None:
        return
    for virus, colour in colour_dict.items():
        virus_col = f"{virus}_detected"
        virus_samples = df.filter(pl.col(virus_col) == 1).select('sample').to_numpy().flatten()
        virus_pca_scores = pca_scores.filter(pl.col('sample').is_in(virus_samples)).to_pandas()
        ax.scatter(virus_pca_scores['PC1'], virus_pca_scores['PC2'], label=virus, color=colour, s=30, alpha=0.8, edgecolor='white', linewidth=0.35, zorder=4)

def plot_loading_barplot(loadings, feature_names, pc_index=0, ax=None, orientation='vertical'):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    order = np.argsort(loadings[:, pc_index], kind='stable')[::-1]
    species_labels = fix_species_labels(feature_names)
    labels = [species_labels[j] for j in order]
    if orientation == 'vertical':
        ax.bar(np.arange(len(order)), loadings[order, pc_index], color='dimgray', width=0.8)
        ax.set_xticks(np.arange(len(order)), labels=labels, rotation=90, fontsize=6.5, fontstyle='italic')
        ax.set_ylabel(f'loading', fontsize=9)
        ax.axhline(0, color='0.7', linewidth=0.6)
    else:
        ax.barh(np.arange(len(order)), loadings[order, pc_index], color='dimgray', height=0.8)
        ax.set_yticks(np.arange(len(order)), labels=labels, fontsize=6.5, fontstyle='italic')
        ax.set_xlabel(f'loading', fontsize=9)
        ax.axvline(0, color='0.7', linewidth=0.6)
    # ax.margins(x=0.01)
    # remove spines, keep ticks
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ax

def fix_species_labels(feature_names):
    lineage_labels = {'ab_S-AUS_Cx_annu': 'Cx. annulirostris (S-AUS)',
                      'ab_AUS_Cx_annu': 'Cx. annulirostris (AUS)'}
    species_labels = [lineage_labels[name] if name in lineage_labels else
                      shorted_genus(name.removeprefix('ab_').replace('_', ' '))
                      for name in feature_names]
    return species_labels

def annotate_scatter_with_virus_centroids(df, pca_scores, ax, colour_dict=None, label_locations=None):
    if colour_dict is None:
        return
    if label_locations is None:
        # (horizontal, vertical) offsets in points; +x is right, +y is up
        # on the page, including when the scatter's x-axis is inverted.
        label_locations = {
            'MVEV': (-23, 4), 'JEV': (-6, 18), 'KUN': (8, 8),
            'RRV': (0, 18), 'BFV': (22, 4), 'SINDV': (0, 18),
            'STRV': (10, 18), 'KOK': (6, -18), 'EHV': (-25, -12), 'LSV': (-10, -18),
        }
    for virus, colour in colour_dict.items():
        virus_col = f"{virus}_detected"
        virus_samples = df.filter(pl.col(virus_col) == 1).select('sample').to_numpy().flatten()
        virus_pca_scores = pca_scores.filter(pl.col('sample').is_in(virus_samples)).to_pandas()
        if not virus_pca_scores.empty:
            centroid_x = virus_pca_scores['PC1'].mean()
            centroid_y = virus_pca_scores['PC2'].mean()
            n = len(virus_pca_scores)
            label = f'{virus}\n(n={n})'
            ax.scatter(centroid_x, centroid_y, s=180, color=colour, alpha=0.65,
                       edgecolor='black', linewidth=1.2, zorder=5)
            dx, dy = label_locations.get(virus, (8, 8))
            ax.annotate(label, xy=(centroid_x, centroid_y),
                        xytext=(dx, dy), textcoords='offset points',
                        ha='center', va='center',
                        fontsize=7, fontweight='bold', color=colour, zorder=6,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

def draw_PCA_on_main_figure(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names, ax):
    # plot PCA scatter
    sns.scatterplot(data=pca_scores.to_pandas(), x='PC1', y='PC2', ax=ax, s=10, color='black', alpha=0.8, edgecolor='white', linewidth=0.35, zorder=3)
    virus_colours = {
        'JEV': 'tab:orange', 'KUN': 'tab:green', 'RRV': 'tab:blue',
        'BFV': 'tab:purple', 'SINDV': 'tab:brown', 'STRV': 'tab:pink',
        'KOK': 'tab:cyan', 'EHV': 'tab:olive', 'LSV': 'black','MVEV': 'tab:red'
    }
    plot_virus_positives_on_PCA_with_colours(df, pca_scores, fig=ax.figure, ax=ax, colour_dict=virus_colours)
    ax.set_xlabel(f'-PC1 ({explained_variance[0]:.1%})' if explained_variance is not None else '-PC1', fontsize=9)
    ax.set_ylabel(f'PC2 ({explained_variance[1]:.1%})' if explained_variance is not None else 'PC2', fontsize=9, labelpad=15)
    # remove spines, ticks and tick labels, add faint grey lines at x=0 and y=0
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis='both', length=0, labelsize=7)
    original_xlim, original_ylim = ax.get_xlim(), ax.get_ylim()
    # Fixed data-coordinate segments keep the grey lines out of the padding.
    ax.plot(original_xlim, [0, 0], color='lightgrey', linewidth=0.5, zorder=1)
    ax.plot([0, 0], original_ylim, color='lightgrey', linewidth=0.5, zorder=1)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    annotate_scatter_with_virus_centroids(df, pca_scores, ax, colour_dict=virus_colours)
    padding = 0.0  # Extra space on each side, relative to the original axis span.
    dx = (original_xlim[1] - original_xlim[0]) * padding
    dy = (original_ylim[1] - original_ylim[0]) * padding
    ax.set_xlim(original_xlim[0] - dx, original_xlim[1] + dx)
    ax.set_ylim(original_ylim[0] - dy, original_ylim[1] + dy)

def draw_beta_forrest_plot(assoc_coefficients, pc_index=0, variables_to_plot=None, ax=None):
    from matplotlib.transforms import blended_transform_factory
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))
    pc = f'PC{pc_index + 1}'
    data = assoc_coefficients.filter(pl.col('PC') == pc)
    if variables_to_plot is not None:
        data = data.filter(pl.col('variable').is_in(variables_to_plot))
    groups = {
        'Environmental': ['distance_to_coast_km_log10', 'E_or_W_of_GDR', 'location_type', 'elevation_m', 'rainfall_mm_3w'],
        'Temporal': ['season', 'week'],
    }
    labels = {'distance_to_coast_km_log10': 'Coast distance (log₁₀)',
              'elevation_m': 'Elevation', 'rainfall_mm_3w': 'Rainfall', 'week': 'Week',
              'E_or_W_of_GDR': 'W vs E of GDR',}
    rows, headings = [], []
    y = 0.0
    for group, variables in groups.items():
        group_rows = []
        for variable in variables:
            for row in data.filter(pl.col('variable') == variable).sort('term').iter_rows(named=True):
                if row['reference_level'] is not None:
                    level = row['term'].split('[T.', 1)[1].removesuffix(']').replace('_', ' ')
                    reference = row['reference_level'].replace('_', ' ')
                    if variable == 'season':
                        level = f'20{level[:2]}-{level[2:]}'
                        reference = f'20{reference[:2]}-{reference[2:]}'
                    label = f'{level} vs {reference}'
                    if variable == 'location_type':
                        label = label.replace('Metropolitan', 'Met.').replace('Met Fringe', 'fringe')
                    elif variable == 'E_or_W_of_GDR':
                        label = f'GDR: {level} vs {reference}'
                else:
                    label = labels.get(variable, variable)
                group_rows.append((row, label))
        if not group_rows:
            continue
        headings.append((group, y))
        y += 0.85
        for row, label in group_rows:
            rows.append((row, label, y))
            y += 1
        y += 0.6  # Small gap before the next group.
    ax.set_axis_off()
    if not rows:
        ax.text(0.5, 0.5, f'No coefficients for {pc}', ha='center', transform=ax.transAxes)
        return ax
    # Keep labels, forest, R² and p inside the supplied panel, rather than
    # letting long tick labels push neighbouring panels apart.
    forest = ax.inset_axes([0.52, 0.09, 0.23, 0.80])
    forest.set_ylim(y - 0.2, -0.9)
    text_transform = blended_transform_factory(ax.transAxes, forest.transData)
    ax.text(0.02, -1.15, '-PC1' if pc == 'PC1' else pc, fontsize=9, fontweight='bold', transform=text_transform)
    ax.text(0.88, -1.15, r'$R^2$', fontsize=8, ha='right', transform=text_transform)
    ax.text(1.04, -1.15, 'Adj. p', fontsize=8, ha='right', transform=text_transform)
    for group, position in headings:
        ax.text(0.02, position - 0.2, group, fontsize=8, fontweight='bold',
                va='center', transform=text_transform)
    for row, label, position in rows:
        beta, low, high = row['beta'], row['ci_lower'], row['ci_upper']
        significant = (row['p_value_fdr'] is not None and row['p_value_fdr'] < 0.05
                       and (low > 0 or high < 0))
        forest.errorbar(beta, position, xerr=[[beta-low], [high-beta]],
                        fmt='o', markersize=6 if significant else 4, color='black', ecolor='0.65', elinewidth=1.4,
                        markerfacecolor='black' if significant else 'white',
                        markeredgecolor='black')
        ax.text(0.02, position, label, fontsize=7, va='center', transform=text_transform)
    # Model-level values appear once, centred across each variable's contrasts.
    for variable in dict.fromkeys(row['variable'] for row, _, _ in rows):
        variable_rows = [(row, position) for row, _, position in rows if row['variable'] == variable]
        row = variable_rows[0][0]
        first, last = variable_rows[0][1], variable_rows[-1][1]
        centre = (first + last) / 2
        r_squared, p_value = row['r_squared'], row['p_value_fdr']
        ax.text(0.88, centre, f'{r_squared:.2f}'.removeprefix('0') if r_squared is not None else 'NA',
                fontsize=7, ha='right', va='center', transform=text_transform)
        p_label = 'NA' if p_value is None else ('<.001' if p_value < 0.001 else f'{p_value:.3f}'.removeprefix('0'))
        ax.text(1.04, centre, p_label, fontsize=7, ha='right', va='center', transform=text_transform)
        if row['reference_level'] is not None and len(variable_rows) > 1:
            # Square bracket to the left of the shared model statistics.
            for x in [0.80]:
                ax.plot([x - 0.012, x, x, x - 0.012], [first, first, last, last],
                        color='0.35', linewidth=0.7, transform=text_transform, clip_on=False)
    forest.axvline(0, color='0.7', linewidth=0.6, linestyle='--')
    forest.set_yticks([])
    forest.tick_params(axis='x', labelsize=7, length=2)
    from matplotlib.ticker import MaxNLocator
    forest.xaxis.set_major_locator(MaxNLocator(2))
    forest.set_xlabel('Beta (95% CI)', fontsize=7, labelpad=2)
    for spine in forest.spines.values():
        spine.set_visible(False)
    return ax

def figure_for_main_text_v2(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names):
    # gridspec 3 cols and 4 rows
    # Fixed panel positions prevent automatic layout from widening the gaps.
    fig = plt.figure(figsize=(11.69, 7.5))
    gs = fig.add_gridspec(4, 3, width_ratios=[0.3, 1.7, 1.4],
                          height_ratios=[1, 1, 1, 1], left=0.27, right=0.98,
                          bottom=0.30, top=0.94, wspace=0.13, hspace=0.35)
    # first col, first three rows
    PC2_loadings_ax = fig.add_subplot(gs[0:3, 0])
    PC1_loadings_ax = fig.add_subplot(gs[3, 1])
    bottom_left_corner_ax = fig.add_subplot(gs[3, 0])
    PCA_ax = fig.add_subplot(gs[0:3, 1])
    PC1_assoc_ax = fig.add_subplot(gs[0:2, 2])
    PC2_assoc_ax = fig.add_subplot(gs[2:4, 2])
    # empty botom left corner
    bottom_left_corner_ax.axis('off')
    bottom_left_corner_ax.text(0.5, 0.5, 'Loading', ha='center', va='center',
                               fontsize=9, transform=bottom_left_corner_ax.transAxes)
    # Bottom strip: species along x, with bars extending vertically.
    plot_loading_barplot(pca.components_.T * np.sqrt(pca.explained_variance_), feature_names, pc_index=0, ax=PC1_loadings_ax, orientation='vertical')
    # Fix the labels
    plt.setp(PC1_loadings_ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    # Left strip: species along y, with bars extending horizontally.
    plot_loading_barplot(pca.components_.T * np.sqrt(pca.explained_variance_), feature_names, pc_index=1, ax=PC2_loadings_ax, orientation='horizontal')
    PC2_loadings_ax.invert_yaxis()
    plt.setp(PC2_loadings_ax.get_yticklabels(), rotation=0, fontsize=7)
    from matplotlib.ticker import MultipleLocator
    loading_values = pca.components_.T * np.sqrt(pca.explained_variance_)
    for i, loading_ax in enumerate([PC1_loadings_ax, PC2_loadings_ax]):
        low, high = loading_values[:, i].min(), loading_values[:, i].max()
        padding = (high - low) * 0.08
        if i == 0:
            loading_ax.set_ylim(low - padding, high + padding)
            loading_ax.yaxis.set_major_locator(MultipleLocator(0.1))
            loading_ax.tick_params(axis='y', labelsize=6)
            loading_ax.set_ylabel('')
        else:
            loading_ax.set_xlim(low - padding, high + padding)
            loading_ax.xaxis.set_major_locator(MultipleLocator(0.1))
            loading_ax.tick_params(axis='x', labelsize=6)
            loading_ax.set_xlabel('')
    # plot PCA scatter
    draw_PCA_on_main_figure(df, pca_scores, explained_variance, pca, assoc_coefficients, variables_to_plot, out_path, feature_names, ax=PCA_ax)


    draw_beta_forrest_plot(assoc_coefficients, pc_index=0,
                           variables_to_plot=variables_to_plot, ax=PC1_assoc_ax)
    draw_beta_forrest_plot(assoc_coefficients, pc_index=1,
                           variables_to_plot=variables_to_plot, ax=PC2_assoc_ax)

    # Scale every label consistently, including text in the inset forest axes.
    from matplotlib.text import Text
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize() * 1.25)

    # savefig
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(out.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=RELEASE_DATA_DIR)
    parser.add_argument('--output', type=Path, default=RESULTS_DIR / 'figure5.png')
    args = parser.parse_args()
    df, scores, variance, pca, coefficients, features = read_release_pca(args.data_dir)
    variables = ['season', 'location_type', 'distance_to_coast_km_log10', 'elevation_m', 'rainfall_mm_3w', 'week', 'E_or_W_of_GDR']
    figure_for_main_text_v2(df, scores, variance, pca, coefficients, variables, args.output, features)


if __name__ == '__main__':
    main()
