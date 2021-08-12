#! /usr/bin/env python3
"""

"""
import logging
from collections import Counter, OrderedDict
PROGRAM = "bactopia tools summary amrfinder"
VERSION = "1.7.1"
IGNORE_LIST = ['.nextflow', '.nextflow.log', 'bactopia-info', 'work', 'bactopia-tools']
EXCLUDE_LIST = []


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO

def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())

def parse_amr_report(report_file, subclass=False):
    """Parse through AMRFinder results."""
    hits = []
    classes = {}
    first_line = True
    cols = None
    with open(report_file, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if first_line:
                cols = line.split('\t')
                first_line = False
            else:
                row = dict(zip(cols, line.split('\t')))
                hits.append(row)
                if subclass:
                    classes[row['Subclass']] = True
                else:
                    classes[row['Class']] = True

    return [hits, list(classes.keys()), cols]

if __name__ == '__main__':
    import argparse as ap
    import os
    import sys
    from pathlib import Path
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Aggregate AMRFinderPlus results'
    )
    parser.add_argument(
        'bactopia', metavar="BACTOPIA_DIRECTORY", type=str,
        help='Directory containing Bactopia output.'
    )

    parser.add_argument('--exclude', metavar="STR", type=str,
                        help='A list of samples to exclude from summary report.')
    parser.add_argument(
        '--subclass', action='store_true',
        help='Group the report by subclass (ex. streptomycin). Default: group by class (ex. aminoglycoside)'
    )
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

    datasets = {'gene': [], 'protein': []}
    with os.scandir(args.bactopia) as dirs:
        for directory in dirs:
            if directory.is_dir():
                if directory.name not in IGNORE_LIST and directory.name not in EXCLUDE_LIST:
                    logging.debug(f"{directory.name} found")
                    amr_directory = f"{args.bactopia}/{directory.name}/antimicrobial-resistance"
                    if os.path.exists(amr_directory):
                        datasets['gene'].append(directory.name)
                        datasets['protein'].append(directory.name)
                    else:
                        logging.debug(f"{directory.name} missing AMRFinderPlus directory")

    
    for dataset, samples in sorted(datasets.items()):
        class_set = set()
        class_counts = Counter()
        results = {}
        sample_class = {}
        columns = None
        logging.debug(f"Processing {dataset} ({len(samples)} total samples)")
        for sample in sorted(samples):
            sample_class[sample] = {}
            amr_report = f"{args.bactopia}/{sample}/antimicrobial-resistance/{sample}-{dataset}-report.txt"
            if os.path.exists(amr_report):
                summary_results, classes, cols = parse_amr_report(amr_report, subclass=args.subclass)
                if not columns:
                    columns = cols
                results[sample] = summary_results
                for c in classes:
                    sample_class[sample][c] = True
                    class_counts[c] += 1
                class_set.update(classes)

        summary_report = []
        row = ['sample_name']
        for c in sorted(class_set):
            row.append(c)
        summary_report.append(row)
        class_hits = Counter()
        for sample, result in sorted(results.items()):
            row = [sample]
            for c in sorted(class_set):
                if c in sample_class[sample]:
                    row.append("True")
                    class_hits[c] += 1
                else:
                    row.append("False")
            summary_report.append(row)

        with open(f"amrfinder-{dataset}-summary.txt", 'w') as fh_out:
            for row in summary_report:
                row_string = '\t'.join(row)
                fh_out.write(f'{row_string}\n')

        with open(f"amrfinder-{dataset}-detailed-summary.txt", 'w') as fh_out:
            row = ['sample']
            for col in columns:
                row.append(col)
            fh_out.write('\t'.join(row))
            fh_out.write('\n')
            for sample, result in sorted(results.items()):
                for hit in result:
                    row = [sample]
                    for col in cols:
                        row.append(hit[col])
                    fh_out.write('\t'.join(row))
                    fh_out.write('\n')
