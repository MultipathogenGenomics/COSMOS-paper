"""Figure 2: mitochondrial enrichment, viral depth and Ct/coverage curve."""
import argparse
from pathlib import Path
import polars as pl
from analysis.common.paths import RESULTS_DIR
from analysis.common.release_data_parsers import (
    RELEASE_DATA_DIR, read_release_virus_plot_data, read_release_virus_coverage,
    read_mitochondrial_coverage,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=RELEASE_DATA_DIR)
    parser.add_argument('--mitochondrial-coverage-file', type=Path,
                        help='Sample,position,coverage CSV/TSV; defaults to data/enrichment/Cx_annu_mitochondria_coverage_plotting_data.csv')
    parser.add_argument('--output', type=Path, default=RESULTS_DIR / 'figure2.png')
    args = parser.parse_args()
    mitochondrial_file = args.mitochondrial_coverage_file or (
        args.data_dir / 'enrichment/Cx_annu_mitochondria_coverage_plotting_data.csv')
    if not mitochondrial_file.is_file():
        parser.error(f'Mitochondrial coverage is required: {mitochondrial_file}. Supply --mitochondrial-coverage-file.')
    from analysis.figures.scripts.enrichment_and_estimation_figure import plot_enrichment_and_nc2
    data = read_release_virus_plot_data(args.data_dir)
    # Match the original curve's seven PCR-tested viruses and zero handling.
    ct_data = data.filter(pl.col('virus_acronym').is_in(
        ['BFV','EHV','JEV','KOK','MVEV','RRV','STRV']) &
        (pl.col('ct_value')>0) & pl.col('ct_value').is_finite()).with_columns(pl.col('nc2').fill_null(0))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_enrichment_and_nc2(
        None, None, None, 'released virus detections', str(args.output),
        mitochondrial_data=read_mitochondrial_coverage(mitochondrial_file),
        virus_coverage_data=read_release_virus_coverage(args.data_dir), positive_ct_data=ct_data)


if __name__ == '__main__':
    main()
