"""Figure 3: virus Ct recovery, ASMMP comparison and mosquito controls."""
import argparse
from pathlib import Path
from analysis.common.paths import RESULTS_DIR
from analysis.common.release_data_parsers import (
    RELEASE_DATA_DIR, read_release_virus_plot_data,
    read_release_morphology_comparison, prepare_release_control_heatmap,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=RELEASE_DATA_DIR)
    parser.add_argument('--output', type=Path, default=RESULTS_DIR / 'figure3.png')
    args = parser.parse_args()
    from analysis.figures.scripts.virus_scatters_and_heatmap import plot_virus_scatters_and_heatmap
    from analysis.validation.scripts.surveillance_validation.process import compute_surveillance_validation_correlations
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mosquito = read_release_morphology_comparison(args.data_dir)
    compute_surveillance_validation_correlations(mosquito).write_csv(
        args.output.with_suffix('.correlations.tsv'), separator='\t')
    plot_virus_scatters_and_heatmap(
        None, None, None, None, None, str(args.output),
        virus_data=read_release_virus_plot_data(args.data_dir), mosquito_data=mosquito,
        control_heatmap=prepare_release_control_heatmap(args.data_dir),
        stats_output_file=str(args.output.with_suffix('.virus_regression.tsv')))


if __name__ == '__main__':
    main()
