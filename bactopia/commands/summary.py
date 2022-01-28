import logging
import os
from collections import defaultdict
import bactopia
from bactopia.summary import get_rank, gather_results, print_failed, print_cutoffs
from bactopia.parse import parse_bactopia_files
from bactopia.const import IGNORE_LIST

PROGRAM = 'bactopia summary'
VERSION = bactopia.__version__
COUNTS = defaultdict(int)
FAILED = defaultdict(list)
CATEGORIES = defaultdict(list)
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


def process_errors(name: str, errors: dict):
    """
    Process a set of errors.

    Args:
        name (str): the sample name
        errors (dict): Errors encountered during processing (keys: 'error_type', 'description')
    """
    error_msg = []
    for error in errors:
        error_msg.append(error['description'])
        COUNTS[error['error_type']] += 1
        FAILED[error['error_type']].append(name)
    COUNTS['total-excluded'] += 1
    COUNTS['qc-failure'] += 1
    CATEGORIES['failed'].append([name, f"Not processed, reason: {';'.join(error_msg)}"])
    logging.debug(f"\t{name}: Not processed, reason: {';'.join(error_msg)}")
    return None


def increment_and_append(key: str, name: str) -> None:
    """
    Increment COUNTS and append to CATEGORIES.

    Args:
        key (str): The key value to use
        name (str): The value to append
    """
    COUNTS[key] += 1
    CATEGORIES[key].append(name)


def process_sample(sample: dict, rank_cutoff: dict) -> dict:
    """
    Process the results of a sample.

    Args:
        sample (dict): all the parsed results associated with a sample
        rank_cutoff (dict): the set of cutoffs for each rank

    Returns:
        list: 0: the sample rank, [description]
    """
    qc = sample['results']['quality-control']['final']['qc_stats']
    assembly = sample['results']['assembly']['stats']
    rank, reason = get_rank(
        rank_cutoff, qc['coverage'], qc['qual_mean'], qc['read_mean'],
        assembly['total_contig'], sample['genome_size'], sample['is_paired']
    )
    increment_and_append('processed', sample['sample'])
    increment_and_append(rank, sample['sample'])

    if rank == 'exclude':
        COUNTS['total-excluded'] += 1
        FAILED['failed-cutoff'].append(sample['sample'])
        CATEGORIES['failed'].append([sample['sample'], f'Failed to pass minimum cutoffs, reason: {reason}'])
    else:
        COUNTS['pass'] += 1

    return gather_results(sample, rank, reason)


def main():
    import argparse as ap
    import csv
    import json
    import sys
    import textwrap

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Create a summary of Bactopia outputs',
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} BACTOPIA_DIRECTORY
        ''')
    )

    parser.add_argument(
        'bactopia', metavar="BACTOPIA_DIRECTORY", type=str,
        help='Directory containing Bactopia output.'
    )

    group1 = parser.add_argument_group('Gold Cutoffs')
    group1.add_argument(
        '--gold_coverage', metavar="FLOAT", type=float, default=100,
        help='Minimum amount of coverage required for Gold status (Default: 100x)'
    )
    group1.add_argument(
        '--gold_quality', metavar="INT", type=int, default=30,
        help='Minimum per-read mean quality score required for Gold status (Default: Q30)'
    )

    group1.add_argument(
        '--gold_read_length', metavar="INT", type=int, default=95,
        help='Minimum mean read length required for Gold status (Default: 95bp)'
    )

    group1.add_argument(
        '--gold_contigs', metavar="INT", type=int, default=100,
        help='Maximum contig count required for Gold status (Default: 100)'
    )

    group2 = parser.add_argument_group('Silver Cutoffs')
    group2.add_argument(
        '--silver_coverage', metavar="FLOAT", type=float, default=50,
        help='Minimum amount of coverage required for Silver status (Default: 50x)'
    )
    group2.add_argument(
        '--silver_quality', metavar="INT", type=int, default=20,
        help='Minimum per-read mean quality score required for Silver status (Default: Q20)'
    )

    group2.add_argument(
        '--silver_read_length', metavar="INT", type=int, default=75,
        help='Minimum mean read length required for Silver status (Default: 75bp)'
    )

    group2.add_argument(
        '--silver_contigs', metavar="INT", type=int, default=200,
        help='Maximum contig count required for Silver status (Default: 200)'
    )

    group3 = parser.add_argument_group('Fail Cutoffs')
    group3.add_argument(
        '--min_coverage', metavar="FLOAT", type=float, default=20,
        help='Minimum amount of coverage required to pass (Default: 20x)'
    )
    group3.add_argument(
        '--min_quality', metavar="INT", type=int, default=12,
        help='Minimum per-read mean quality score required to pass (Default: Q12)'
    )

    group3.add_argument(
        '--min_read_length', metavar="INT", type=int, default=49,
        help='Minimum mean read length required to pass (Default: 49bp)'
    )

    group3.add_argument(
        '--max_contigs', metavar="INT", type=int, default=500,
        help='Maximum contig count required to pass (Default: 500)'
    )

    group3.add_argument(
        '--min_assembled_size', metavar="FLOAT", type=float,
        help='Minimum assembled genome size.'
    )

    group3.add_argument(
        '--max_assembled_size', metavar="FLOAT", type=float,
        help='Maximum assembled genome size.'
    )

    group5 = parser.add_argument_group('Helpers')
    group5.add_argument(
        '--outdir', metavar="OUTPUT_DIRECTORY", type=str, default="./",
        help='Directory to write output. (Default: ./)'
    )

    group5.add_argument(
        '--prefix', metavar="STR", type=str, default="bactopia",
        help='Prefix to use for output files. (Default: bactopia)'
    )
    group5.add_argument('--force', action='store_true',
                        help='Overwrite existing reports.')
    group5.add_argument('--depends', action='store_true',
                        help='Verify dependencies are installed.')
    group5.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    group5.add_argument('--verbose', action='store_true',
                        help='Increase the verbosity of output.')
    group5.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if os.path.exists(f'{args.outdir}/{args.prefix}-exclude.txt') and not args.force:
        print(f"Existing reports found in {args.outdir}. Will not overwirte unless --force is used. Exiting.",
              file=sys.stderr)
        sys.exit(1)

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(logging.ERROR if args.silent else logging.DEBUG if args.verbose else logging.INFO)
    RANK_CUTOFF = {
        'gold': {
            'coverage': args.gold_coverage,
            'quality': args.gold_quality,
            'length': args.gold_read_length,
            'contigs': args.gold_contigs
        },
        'silver': {
            'coverage': args.silver_coverage,
            'quality': args.silver_quality,
            'length': args.silver_read_length,
            'contigs': args.silver_contigs
        },
        'bronze': {
            'coverage': args.min_coverage,
            'quality': args.min_quality,
            'length': args.min_read_length,
            'contigs': args.max_contigs
        },
        'min-assembled-size': args.min_assembled_size,
        'max-assembled-size': args.max_assembled_size
    }

    processed_samples = {}
    fields = []
    results = []
    logging.debug(f"Working on {args.bactopia}...")
    with os.scandir(args.bactopia) as dirs:
        for i, directory in enumerate(dirs):
            if directory.name not in IGNORE_LIST:
                logging.debug(f"Working on {directory.name} ({i+1})")
                sample = parse_bactopia_files(args.bactopia, directory.name)
                if sample['ignored']:
                    logging.debug(f"\t{sample['sample']} is not a Bactopia directory, ignoring...")
                    increment_and_append('ignore-unknown', sample['sample'])
                else:
                    COUNTS['total'] += 1
                    if sample['has_errors']:
                        process_errors(sample['sample'], sample['errors'])
                    else:
                        processed = process_sample(sample, RANK_CUTOFF)
                        fields = list(dict.fromkeys(fields + list(processed.keys())))
                        results.append(processed)
                        processed_samples[sample['sample']] = True

    # Write outputs
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Tab-delimited report
    txt_report = f'{outdir}/{args.prefix}-report.txt'
    with open(txt_report, 'w') as txt_fh:
        outputs = []
        for result in results:
            output = {}
            for field in fields:
                if field in result:
                    if isinstance(result[field], float):
                        output[field] = f"{result[field]:.3f}"
                    else:
                        output[field] = result[field]
                else:
                    output[field] = ""
            outputs.append(output)
        if outputs:
            writer = csv.DictWriter(txt_fh, fieldnames=outputs[0].keys(), delimiter='\t')
            writer.writeheader()
            for row in outputs:
                writer.writerow(row)

    # Exclusion report
    exclusion_report = f'{outdir}/{args.prefix}-exclude.txt'
    cutoff_counts = defaultdict(int)
    with open(exclusion_report, 'w') as exclude_fh:
        exclude_fh.write('sample\tstatus\treason\n')
        for name, reason in CATEGORIES['failed']:
            if name in processed_samples:
                reasons = reason.split(':')[1].split(';')
                cutoffs = []
                for r in reasons:
                    cutoffs.append(r.split('(')[0].strip().title())
                cutoff_counts[';'.join(sorted(cutoffs))] += 1
                exclude_fh.write(f'{name}\texclude\t{reason}\n')
            else:
                exclude_fh.write(f'{name}\tqc-fail\t{reason}\n')

    # Screen report
    summary_report = f'{outdir}/{args.prefix}-summary.txt'
    with open(summary_report, 'w') as summary_fh:
        summary_fh.write("Bactopia Summary Report\n")
        summary_fh.write(textwrap.dedent(f'''
            Total Samples: {COUNTS['total']}
            
            Passed: {COUNTS["pass"]}
                Gold: {COUNTS["gold"]}
                Silver: {COUNTS["silver"]}
                Bronze: {COUNTS["bronze"]}

            Excluded: {COUNTS["total-excluded"]}
                Failed Cutoff: {COUNTS["exclude"]}\n'''))
        summary_fh.write(print_cutoffs(cutoff_counts))
        summary_fh.write(f'    QC Failure: {COUNTS["qc-failure"]}\n')
        summary_fh.write(print_failed(FAILED))
        summary_fh.write(textwrap.dedent(f'''
            Reports:
                Full Report (txt): {txt_report}
                Exclusion: {exclusion_report}
                Summary: {summary_report}

            Rank Cutoffs:
                Gold:
                    Coverage >= {RANK_CUTOFF['gold']['coverage']}x
                    Quality >= Q{RANK_CUTOFF['gold']['quality']}
                    Read Length >= {RANK_CUTOFF['gold']['length']}bp
                    Total Contigs < {RANK_CUTOFF['gold']['contigs']}
                Silver:
                    Coverage >= {RANK_CUTOFF['silver']['coverage']}x
                    Quality >= Q{RANK_CUTOFF['silver']['quality']}
                    Read Length >= {RANK_CUTOFF['silver']['length']}bp
                    Total Contigs < {RANK_CUTOFF['silver']['contigs']}
                Bronze:
                    Coverage >= {RANK_CUTOFF['bronze']['coverage']}x
                    Quality >= Q{RANK_CUTOFF['bronze']['quality']}
                    Read Length >= {RANK_CUTOFF['bronze']['length']}bp
                    Total Contigs < {RANK_CUTOFF['bronze']['contigs']}
            
            Assembly Length Exclusions:
                Minimum: {RANK_CUTOFF['min-assembled-size']}
                Maximum: {RANK_CUTOFF['min-assembled-size']}
        '''))


if __name__ == '__main__':
    main()
