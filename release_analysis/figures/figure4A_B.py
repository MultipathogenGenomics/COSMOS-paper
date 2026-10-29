"""Figure 4A–B: released mosquito/virus results and lineage map."""
import argparse
from pathlib import Path
from analysis.common.paths import RESULTS_DIR
from analysis.common.release_data_parsers import (
    RELEASE_DATA_DIR, prepare_release_results_figure, read_release_lineages,
    read_release_metadata,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=RELEASE_DATA_DIR)
    parser.add_argument('--map-file', type=Path, required=True,
                        help='ABS GeoPackage with STATE_NAME_2021 and NSW/ACT polygons')
    parser.add_argument('--sydney-inset-bounds', type=float, nargs=4,
                        default=(150.5, 151.43, -34.1, -33.5),
                        metavar=('XMIN', 'XMAX', 'YMIN', 'YMAX'),
                        help='Inset extent in longitude/latitude; release lacks historical council labels')
    parser.add_argument('--site-averages-file', type=Path, help='Optional TSV export of map groups')
    parser.add_argument('--output', type=Path, default=RESULTS_DIR / 'figure4A-B.png')
    args = parser.parse_args()
    if not args.map_file.is_file():
        parser.error(f'Map file not found: {args.map_file}')
    from analysis.figures.scripts.results_figure import plot_results_figure
    df = prepare_release_results_figure(args.data_dir)
    lineages = read_release_lineages(args.data_dir)
    points = lineages.filter(lineages['included_in_figure4']).select(
        'sample', 'avg_subtype_ratio').join(read_release_metadata(args.data_dir), on='sample', validate='1:1')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_results_figure(None, None, None, None, str(args.map_file), str(args.output),
                        prepared_data=df, release_lineage_points=points,
                        sydney_inset_bounds=args.sydney_inset_bounds,
                        site_averages_file=args.site_averages_file)


if __name__ == '__main__':
    main()
