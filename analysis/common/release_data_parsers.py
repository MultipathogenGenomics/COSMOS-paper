"""Readers and schema adapters for the frozen data release.

Add figure-specific parsers here as figures are connected to release files.
Keep release inputs read-only and write derived outputs under results/.
"""

from __future__ import annotations

from pathlib import Path

from analysis.common.paths import REPO_ROOT
RELEASE_DATA_DIR = REPO_ROOT / "data"


def release_data_path(relative_path: str | Path) -> Path:
    """Resolve a path relative to the repository data release."""
    return RELEASE_DATA_DIR / relative_path


def read_release_metadata(data_dir=RELEASE_DATA_DIR):
    """Read sample metadata; location-site also supplies the display location."""
    import polars as pl
    df = pl.read_csv(Path(data_dir) / 'samples/sample_metadata.csv')
    if df['sample'].n_unique() != df.height:
        raise ValueError('Release metadata must have one row per sample')
    return df.with_columns(pl.col('location-site').alias('location'))


def read_release_viruses(data_dir=RELEASE_DATA_DIR):
    """Preserve absent metrics as nulls; detection flags remain authoritative."""
    import polars as pl
    return pl.read_csv(Path(data_dir) / 'virus_detections/virus_detections_long_form.csv')


def read_release_abundance(data_dir=RELEASE_DATA_DIR):
    """Read the published species proportions, retaining release precision."""
    import polars as pl
    df = pl.read_csv(Path(data_dir) / 'mosquito_species_calls/tableS2.mosquito_species_calls_summary.csv')
    species = [c for c in df.columns if c not in {
        'sample', 'batch', 'location-site', 'season', 'coastal', 'region_type',
        'simpsons_diversity_index'}]
    return df.select(['sample'] + species).rename({c: c.replace('_', ' ') for c in species})


def read_release_lineages(data_dir=RELEASE_DATA_DIR):
    import polars as pl
    base = Path(data_dir) / 'cryptic_mosquito_subpopulations/Cx_annulirostris'
    summary = pl.read_csv(base / 'tableS6.Cx_annulirostris_lineage_informative_snp_summary.csv')
    counts = pl.read_csv(base / 'Cx_annulirostris_lineage_informative_snp_counts.csv')
    depth = counts.group_by('sample').agg(pl.col('depth_at_position').min().alias('min_count'))
    return summary.join(depth, on='sample', how='left', validate='1:1').rename({
        'mean_S-AUS_subtype_ratio': 'avg_subtype_ratio'})


def prepare_release_results_figure(data_dir=RELEASE_DATA_DIR):
    """Assemble Figure 4 inputs without historical batch directories."""
    import polars as pl
    from analysis.figures.scripts.process import order_results_data_for_plotting
    metadata = read_release_metadata(data_dir)
    abundance = read_release_abundance(data_dir)
    diversity = pl.read_csv(Path(data_dir) / 'mosquito_species_calls/tableS2.mosquito_species_calls_summary.csv').select(
        'sample', pl.col('simpsons_diversity_index').alias('Simpsons_Diversity_Index'))
    viruses = read_release_viruses(data_dir).pivot(on='virus_acronym', index='sample', values='prop_covered_depth2')
    cts = read_release_viruses(data_dir).with_columns(
        pl.col('ct').cast(pl.Float64, strict=False)).pivot(
            on='virus_acronym', index='sample', values='ct')
    cts = cts.rename({c: 'ct_'+c for c in cts.columns if c!='sample'})
    lineages = read_release_lineages(data_dir).with_columns(
        pl.when(pl.col('included_in_figure4')).then(pl.col('avg_subtype_ratio'))
        .otherwise(None).alias('avg_subtype_ratio'))
    df = (metadata.join(diversity, on='sample', validate='1:1')
          .join(abundance, on='sample', validate='1:1')
          .join(viruses, on='sample', validate='1:1')
          .join(cts, on='sample', validate='1:1')
          .join(lineages.select('sample', 'avg_subtype_ratio'), on='sample', validate='1:1')
          .rename({'Culex annulirostris': 'Culex annulirostris proportion'})
          .with_columns((pl.col('distance_to_coast_km') < 20).cast(pl.Int8).alias('Coastal'))
          .with_columns(
              (pl.col('Culex annulirostris proportion') * pl.col('avg_subtype_ratio')).alias('Culex annulirostris S-AUS'),
              (pl.col('Culex annulirostris proportion') * (1-pl.col('avg_subtype_ratio'))).alias('Culex annulirostris AUS')))
    return order_results_data_for_plotting(df)


def read_release_pca(data_dir=RELEASE_DATA_DIR):
    """Use frozen scores/models; recover scaled loadings from score variances.

    Explained-variance ratios are read alongside the frozen PCA tables.
    """
    import numpy as np
    import polars as pl
    from types import SimpleNamespace
    base = Path(data_dir) / 'PCA'
    scores = pl.read_csv(base / 'PCA_scores.csv')
    loadings = pl.read_csv(base / 'PCA_loadings.csv')
    models = pl.read_csv(base / 'PCA_association_with_variables_results.csv')
    pcs = [c for c in scores.columns if c.startswith('PC')]
    eigenvalues = scores.select(pcs).to_numpy().var(axis=0, ddof=1)
    model = SimpleNamespace(components_=loadings.select(pcs).to_numpy().T,
                            explained_variance_=eigenvalues)
    variance = pl.read_csv(base / 'PCA_explained_variance_ratio.csv')
    lookup = dict(zip(variance['PC'], variance['explained_variance_ratio']))
    ratios = np.array([lookup[pc] for pc in pcs], dtype=float)
    if not np.isfinite(ratios).all() or (ratios < 0).any() or ratios.sum() > 1+1e-8:
        raise ValueError('Explained variance ratios must be finite, nonnegative and sum to <=1')
    viruses = read_release_viruses(data_dir).select('sample', 'virus_acronym',
        pl.col('COSMOS_virus_detection').cast(pl.Int8).alias('detected')).pivot(
            on='virus_acronym', index='sample', values='detected')
    viruses = viruses.rename({c: c+'_detected' for c in viruses.columns if c!='sample'})
    return viruses, scores, ratios, model, models, loadings['features'].to_list()


def read_release_virus_plot_data(data_dir=RELEASE_DATA_DIR):
    """Convert release column names to the original virus plotting schema."""
    import polars as pl
    return read_release_viruses(data_dir).rename({
        'ct': 'ct_value', 'prop_covered_depth2': 'nc2', 'num_unique_reads': 'dedup_reads',
    }).with_columns(pl.col('ct_value').cast(pl.Float64, strict=False))


def read_release_virus_coverage(data_dir=RELEASE_DATA_DIR):
    import polars as pl
    return pl.read_csv(Path(data_dir) / 'enrichment/virus_coverage_per_position_per_sample.tsv',
                       separator='\t', schema_overrides={'coverage': pl.Float64})


def read_release_morphology_comparison(data_dir=RELEASE_DATA_DIR):
    import polars as pl
    df = pl.read_csv(Path(data_dir) / 'mosquito_species_calls/COSMOS_vs_ASMMP_two_species.tsv', separator='\t')
    return df.select('sample', 'species', 'estimated_abundance', 'surveillance_proportion').drop_nulls()


def read_mitochondrial_coverage(path):
    """Read a released sample,position,coverage CSV/TSV without changing depth."""
    import polars as pl
    path = Path(path)
    df = pl.read_csv(path, separator=',' if path.suffix=='.csv' else '\t')
    required = {'sample', 'coverage'}
    if not required.issubset(df.columns) or not {'position','aligned_position'}.intersection(df.columns):
        raise ValueError('Mitochondrial coverage needs sample, position (or aligned_position), coverage')
    if 'is_cox1' not in df.columns:
        position = 'position' if 'position' in df.columns else 'aligned_position'
        df = df.with_columns(pl.col(position).is_between(1430,2966).alias('is_cox1'))
    return df


def prepare_release_control_heatmap(data_dir=RELEASE_DATA_DIR):
    """Build the original heatmap structure from released controls and assignments.

    Simulation estimates use both tools, the published count/coverage/identity
    filters, and KMA counts normalised over passing species, as in the original
    abundance workflow. No FASTQ, generated pickle or external annotations needed.
    """
    import re
    import polars as pl
    from analysis.common.mosquito_ids import shorted_genus
    controls = pl.read_csv(Path(data_dir) / 'controls_and_validation/controls_mosquito_abundance_results.tsv', separator='\t')
    columns = [c for c in controls.columns if c.startswith('ab_')]
    names = {c: c.removeprefix('ab_').replace('_',' ').replace('Ochlerotatus ', 'Aedes ') for c in columns}
    # Preserve the genus ordering used in the original validation figure.
    genera = {'Mansonia':0, 'Anopheles':1, 'Coquillettidia':2, 'Aedes':3, 'Culex':4}
    overrides = {'Aedes vigilax':10000, 'Culex annulirostris':-1,
                 'Culex australicus':10000, 'Coquillettidia linealis':10000}
    columns = sorted(columns, key=lambda c:(genera[names[c].split()[0]],overrides.get(names[c],columns.index(c))))
    controls = controls.sort('sample')
    data = {'column_metadata':[{'species':names[c], 'display_species':shorted_genus(names[c]),
                               'genus':names[c].split()[0], 'label':str(i+1)} for i,c in enumerate(columns)]}
    def group(kind):
        return controls.filter(pl.col('type')==kind)
    for key,kind in [('species','single species pool'),('morphology test','morphologically similar taxa'),
                     ('laboratory mixture','laboratory mixture'),('laboratory serial dilution','laboratory serial dilution')]:
        rows = group(kind).to_dicts()
        if key=='species':
            rows.sort(key=lambda r:(genera[r['content'].split()[0]], next((i for i,c in enumerate(columns) if names[c]==r['content']),999)))
        if key=='laboratory serial dilution':
            # Display sitiens:pipiens above vigilax:annulirostris as in the original figure.
            rows.sort(key=lambda r:('vigilax' in r['content'], float(re.findall(r'([\d.]+)%',r['content'])[0])))
        labels=[]
        for i,r in enumerate(rows):
            if key=='species':labels.append(shorted_genus(r['content']))
            elif key=='morphology test':labels.append(f'Mixture {i+1}')
            elif key=='laboratory mixture':
                labels.append('Ae.proc/Ae.mult/Cx.siti/Cx.pipi' if 'procax' in r['content'] else 'Cq.line/Ae.noto/Ae.vigi/Cx.annu')
            else:
                percentages=re.findall(r'([\d.]+)%',r['content'])
                if 'sitiens' in r['content']:percentages=percentages[::-1]
                labels.append(':'.join(percentages))
        data[key]=(pl.DataFrame([{names[c]:r[c] for c in columns} for r in rows]).to_pandas(),labels)
    base=Path(data_dir)/'controls_and_validation/in_silico_mixtures'
    meta=pl.read_csv(base/'in_silico_mixtures_metadata.tsv',separator='\t')
    quant=pl.read_csv(base/'in_silico_mixtures_kallisto_kma_quantification_results_joined.tsv',separator='\t')
    quant=quant.with_columns((pl.col('kma_template_identity') / pl.min_horizontal(
        pl.lit(100.),pl.col('kma_template_coverage')).clip(lower_bound=1e-9)*100).alias('identity_covered'),
        ((pl.col('kma_read_count')>0)&(pl.col('kallisto_estimated_counts')>0)).alias('both_tools'))
    quant=quant.group_by('mixture_id','species').agg(
        pl.col('kma_read_count','kallisto_estimated_counts').sum(),
        pl.col('kma_template_coverage','identity_covered','both_tools').max())
    quant=quant.filter((pl.col('kma_read_count')>=5)&(pl.col('kallisto_estimated_counts')>=5)&
                       (pl.col('kma_template_coverage')>=90)&(pl.col('identity_covered')>=97.5)&pl.col('both_tools'))
    quant=quant.with_columns((pl.col('kma_read_count')/pl.col('kma_read_count').sum().over('mixture_id')).alias('abundance'))
    mixtures=[{'Ochlerotatus_procax','Ochlerotatus_multiplex','Culex_sitiens','Culex_pipiens'},
              {'Coquillettidia_linealis','Aedes_notoscriptus','Aedes_vigilax','Culex_annulirostris'}]
    simulated=[]
    for members in mixtures:
        matched=[r['mixture_id'] for r in meta.to_dicts() if {r[f'sample_{i}'] for i in range(1,5)}==members]
        if len(matched)!=1:raise ValueError('Expected exactly one simulation for the selected species')
        estimates={r['species'].replace('_',' ').replace('Ochlerotatus ','Aedes '):r['abundance']
                   for r in quant.filter(pl.col('mixture_id')==matched[0]).to_dicts()}
        simulated.append({names[c]:estimates.get(names[c],0.) for c in columns})
    data['simulated mixture']=(pl.DataFrame(simulated).to_pandas(),
        ['Ae.proc/Ae.mult/Cx.siti/Cx.pipi','Cq.line/Ae.noto/Ae.vigi/Cx.annu'])
    return data
