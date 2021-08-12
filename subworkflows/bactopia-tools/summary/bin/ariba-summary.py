#! /usr/bin/env python3
"""
usage: bactopia tools summary ariba [-h] [--exclude STR] [--include_all]
                                    [--version]
                                    BACTOPIA_DIRECTORY

bactopia tools summary ariba (v1.2.4) - Aggregate Ariba summaries

positional arguments:
  BACTOPIA_DIRECTORY  Directory containing Bactopia output.

optional arguments:
  -h, --help          show this help message and exit
  --exclude STR       A list of samples to exclude from summary report.
  --include_all       Include all hits, not just the matches.
  --version           show program's version number and exit
"""
import logging
from collections import Counter, OrderedDict
PROGRAM = "bactopia tools summary ariba"
VERSION = "1.7.1"
IGNORE_LIST = ['.nextflow', '.nextflow.log', 'bactopia-info', 'work', 'bactopia-tools']
EXCLUDE_LIST = []

def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO

def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())

def parse_summary(summary_file, include_all=False):
    """Find all Ariba summary reports in a directory."""
    hits = OrderedDict()
    first_line = True
    with open(summary_file, 'rt') as fh:
        cols = None
        for line in fh:
            line = line.rstrip()
            if first_line:
                cols = line.split(',')
                first_line = False
            else:
                row = dict(zip(cols, line.split(',')))
                for key, val in row.items():
                    if key != 'name':
                        cluster, field = key.split('.')
                        if cluster not in hits:
                            hits[cluster] = {
                                'cluster': cluster,
                            }
                        hits[cluster][field] = val

    summary = {}
    for key, val in hits.items():
        if val['match'] == 'yes' or include_all:
            summary[key] = val
    return [summary, list(summary.keys())]

if __name__ == '__main__':
    import argparse as ap
    import os
    import sys
    from pathlib import Path
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Aggregate Ariba summaries'
    )
    parser.add_argument(
        'bactopia', metavar="BACTOPIA_DIRECTORY", type=str,
        help='Directory containing Bactopia output.'
    )

    parser.add_argument('--exclude', metavar="STR", type=str,
                        help='A list of samples to exclude from summary report.')

    parser.add_argument('--include_all', action='store_true',
                        help='Include all hits, not just the matches.')
    parser.add_argument('--verbose', action='store_true',
                        help='Increase the verbosity of output.')
    parser.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.exclude:
        with open(args.exclude, 'rt') as exclude_fh:
            for line in exclude_fh:
                EXCLUDE_LIST.append(line.rstrip().split('\t')[0])
    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    datasets = {}
    with os.scandir(args.bactopia) as dirs:
        for directory in dirs:
            if directory.is_dir():
                if directory.name not in IGNORE_LIST and directory.name not in EXCLUDE_LIST:
                    logging.debug(f"{directory.name} found")
                    ariba_directory = f"{args.bactopia}/{directory.name}/ariba"
                    if os.path.exists(ariba_directory):
                        for dataset in [ f.name for f in os.scandir(ariba_directory) if f.is_dir() ]:
                            logging.debug(f"\t{dataset} found")
                            if dataset not in datasets:
                                datasets[dataset] = []
                            datasets[dataset].append(directory.name)
                    else:
                        logging.debug(f"{directory.name} missing Ariba directory (possibly single-end reads, or ariba not run)")

    for dataset, samples in sorted(datasets.items()):
        cluster_set = set()
        cluster_counts = Counter()
        results = {}
        logging.debug(f"Processing {dataset} ({len(samples)} total samples)")
        for sample in sorted(samples):
            summary_file = f"{args.bactopia}/{sample}/ariba/{dataset}/summary.csv"
            if os.path.exists(summary_file):
                summary_results, clusters = parse_summary(summary_file, include_all=args.include_all)
                results[sample] = summary_results
                for c in clusters:
                    cluster_counts[c] += 1
                cluster_set.update(clusters)

        summary_report = []
        row = ['sample_name']
        for cluster in sorted(cluster_set):
            row.append(cluster)
        summary_report.append(row)
        cluster_hits = Counter()
        for sample, result in sorted(results.items()):
            row = [sample]
            for cluster in sorted(cluster_set):
                if cluster in result:
                    row.append("True")
                    cluster_hits[cluster] += 1
                else:
                    row.append("False")
            summary_report.append(row)

        with open(f"ariba-{dataset}-summary.txt", 'w') as fh_out:
            for row in summary_report:
                row_string = '\t'.join(row)
                fh_out.write(f'{row_string}\n')

        cols = ['cluster', 'assembled', 'match', 'pct_id', 'ctg_cov', 'known_var', 'novel_var']
        with open(f"ariba-{dataset}-detailed-summary.txt", 'w') as fh_out:
            row = ['sample']
            for col in cols:
                row.append(col)
            fh_out.write('\t'.join(row))
            fh_out.write('\n')
            for sample, result in sorted(results.items()):
                for cluster, vals in sorted(result.items()):
                    row = [sample]
                    for col in cols:
                        row.append(vals[col])
                    fh_out.write('\t'.join(row))
                    fh_out.write('\n')
