#! /usr/bin/env python3
"""

"""
import glob
import jinja2
import json
import logging
import os
import sys

from collections import defaultdict
from statistics import mean

PROGRAM = "bactopia tools summary"
VERSION = "1.2.4"
TEMPLATE_DIR = os.path.dirname(os.path.realpath(__file__)).replace('bin', 'templates')

STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")

IGNORE_LIST = ['.nextflow', '.nextflow.log', 'bactopia-info', 'work']
COUNTS = defaultdict(int)
FAILED = defaultdict(list)
CATEGORIES = defaultdict(list)
COUNTS_BY_RANK = {
    'total': defaultdict(list),
    'gold': defaultdict(list),
    'silver': defaultdict(list),
    'bronze': defaultdict(list),
    'fail': defaultdict(list)
}

def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO

def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())

def check_bactopia(path, name):
    """Check if input appears to be a Bactopia directory."""
    found_bactopia = False
    has_error = []
    files = [
        f"{path}/{name}/{name}-genome-size.txt",
        f"{path}/{name}/{name}-genome-size-error.txt",
        f"{path}/{name}/{name}-low-read-count-error.txt",
        f"{path}/{name}/{name}-low-sequence-depth-error.txt",
        f"{path}/{name}/{name}-paired-end-error.txt"
    ]

    for filename in files:
        if os.path.exists(filename):
            found_bactopia = True
            if filename.endswith("error.txt"):
                has_error.append(filename)

    return [found_bactopia, has_error]

def parse_error(error):
    """Return the error name."""
    if error.endswith("genome-size-error.txt"):
        return "genome-size-error"
    elif error.endswith("low-read-count-error.txt"):
        return "low-read-count-error"
    elif error.endswith("low-sequence-depth-error.txt"):
        return "low-sequence-depth-error"
    elif  error.endswith("paired-end-error.txt"):
        return "paired-end-error"
    return "unknown-error"

def parse_read_stats(json_file, ):
    """ """
    with open(json_file) as json_fh:
        return json.load(json_fh)

def parse_json(json_file):
    """Return a dict of the input json_file"""
    with open(json_file) as json_fh:
        return json.load(json_fh)

def parse_annotation(txt_file):
    """Parse Prokka summary text file."""
    results = {}
    with open(txt_file) as txt_fh:
        for line in txt_fh:
            line = line.rstrip()
            key, val = line.split(":")
            results[key] = val.lstrip()
    return results

def get_files(path, sample):
    "Return a list of files to read."
    files = None
    assemblies = []
    missing = []
    end_type = 'paired-end'

    if os.path.exists(f"{path}/{sample}/quality-control/{sample}.fastq.gz"):
        # Single End
        end_type = 'single-end'
        files = {
            'annotation': f"{path}/{sample}/annotation/{sample}.txt",
            'assembly': f"{path}/{sample}/assembly/{sample}.fna.json",
            'original': f"{path}/{sample}/quality-control/summary-original/{sample}-original.json",
            'final': f"{path}/{sample}/quality-control/summary-final/{sample}-final.json",
        }
    else:
        # Paired End
        files = {
            'annotation': f"{path}/{sample}/annotation/{sample}.txt",
            'assembly': f"{path}/{sample}/assembly/{sample}.fna.json",
            'original-r1': f"{path}/{sample}/quality-control/summary-original/{sample}_R1-original.json",
            'original-r2': f"{path}/{sample}/quality-control/summary-original/{sample}_R2-original.json",
            'final-r1': f"{path}/{sample}/quality-control/summary-final/{sample}_R1-final.json",
            'final-r2': f"{path}/{sample}/quality-control/summary-final/{sample}_R2-final.json"
        }

    for key, val in files.items():
        if not os.path.exists(val):
            missing.append(val)

    return {'files': files, 'missing': missing, 'end_type': end_type, 'has_error': False}

def get_rank(cutoff, coverage, quality, length, contigs, is_paired):
    """Return the rank (gold, silver, bronze, fail) of the sample."""
    gold = cutoff['gold']
    silver = cutoff['silver']
    bronze = cutoff['bronze']
    if coverage >= gold['coverage'] and quality >= gold['quality'] and length >= gold['length'] and contigs <= gold['contigs'] and is_paired:
        return 'gold'
    elif coverage >= silver['coverage'] and quality >= silver['quality'] and length >= silver['length'] and contigs <= silver['contigs'] and is_paired:
        return 'silver'
    elif coverage >= bronze['coverage'] and quality >= bronze['quality'] and length >= bronze['length'] and contigs <= bronze['contigs']:
        return 'bronze'
    else:
        return 'fail'

def merge_stats(r1, r2):
    if r2:
        merged = {}
        for key, val in r1['qc_stats'].items():
            if key in ['total_bp', 'coverage', 'read_total']:
                merged[key] = r1['qc_stats'][key] + r1['qc_stats'][key]
            else:
                merged[key] = mean([r1['qc_stats'][key], r1['qc_stats'][key]])
        return merged
    return r1

def add_to_counts(dictionary, rank):
    """Append values to rank counts."""
    for key, val in dictionary.items():
        COUNTS_BY_RANK[rank][key].append(val)
        COUNTS_BY_RANK['total'][key].append(val)

def fastani(samples):
    for key, val in samples.items():
        if val['has_error'] or len(val['missing']):
            #
        else:

if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    from statistics import mean
    import textwrap
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Create a summary of Bacopia outputs'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} bactopiadir outdir
        ''')
    )

    parser.add_argument(
        'bactopia', metavar="BACTOPIA_DIRECTORY", type=str,
        help='Directory containing Bactopia output.'
    )

    parser.add_argument(
        'outdir', metavar="OUTPUT_DIRECTORY", type=str,
        help='Directory to write output.'
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
        '--silver_read_length', metavar="INT", type=int, default=45,
        help='Minimum mean read length required for Silver status (Default: 45bp)'
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
        '--min_read_length', metavar="INT", type=int, default=30,
        help='Minimum mean read length required to pass (Default: 30bp)'
    )

    group3.add_argument(
        '--max_contigs', metavar="INT", type=int, default=500,
        help='Maximum contig count required to pass (Default: 500)'
    )

    group4 = parser.add_argument_group('Taxon Check')
    group4.add_argument(
        '--fastani', action='store_true',
        help='Determine the pairwise ANI of each sample using FastANI.')
    )

    group4.add_argument(
        '--ani_threshold', metavar="FLOAT", type=float, default=0.8
        help='Exclude samples with a mean ANI less than the threshold. (Default: 0.80).'
    )

    group4.add_argument(
        '--min_sourmash', metavar="FLOAT", type=float, default=0.5
        help='Minimum percentage of k-mers that must have matched (Default: 50%).'
    )
    group5.add_argument('--cpus', metavar="INT", type=int, default=1,
                        help='Number of CPUs available for FastANI. (Default: 1)')

    group5 = parser.add_argument_group('Helpers')
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

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))
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
        }
    }

    samples = {}
    with os.scandir(args.bactopia) as dirs:
        for directory in dirs:
            if directory.name in IGNORE_LIST:
                logging.debug(f"Ignoring {directory.name}")
                COUNTS['ignore-list'] += 1
                CATEGORIES['ignore-list'].append(directory.name)
            else:
                is_bactopia, has_error = check_bactopia(args.bactopia, directory.name)
                if not is_bactopia:
                    logging.debug(f"{directory.name} is not a Bactopia directory, ignoring...")
                    COUNTS['ignore-unknown'] += 1
                    CATEGORIES['ignore-unknown'].append(directory.name)
                else:
                    COUNTS['total'] += 1
                    if has_error:
                        logging.debug(f"{directory.name} {has_error} ")
                        for error in has_error:
                            error_msg = parse_error(error)
                            COUNTS[error_msg] += 1
                            COUNTS['error'] += 1
                            COUNTS['exclude'] += 1
                            FAILED[error_msg].append(directory.name)
                            CATEGORIES['exclude'].append([directory.name, f"Not processed due to '{error_msg}' error"])
                            samples[directory.name] = {'has_error': True, 'missing': []}
                    else:
                        logging.debug(f"{directory.name} found ")
                        COUNTS['processed'] += 1
                        CATEGORIES['processed'].append(directory.name)
                        samples[direcotry.name] = get_files(args.bactopia, directory.name)

    if args.fastani:
        # Use FastANI to flag samples that may be another organism
        pairwise_ani = fastani(samples)

    for key, val in samples.items():
        if not val['has_error']:
            end_type = val['end_type']
            files = val['files']
            missing = val['missing']

            COUNTS[end_type] += 1
            if missing:
                COUNTS['missing'] += 1
                CATEGORIES['exclude'].append([directory.name, 'Missing expected files'])
                logging.debug(f"{directory.name} missing files ")
                logging.debug(f"{missing}")
            else:
                COUNTS['found'] += 1
                original_r1 = None
                original_r2 = None
                final_r1 = None
                final_r2 = None
                rank = None
                assembly = parse_json(files['assembly'])
                annotation = parse_annotation(files['annotation'])
                if end_type == 'single-end':
                    original_r1 = parse_json(files['original'])
                    final_r1 = parse_json(files['final'])
                    coverage = final_r1['qc_stats']['coverage']
                else:
                    original_r1 = parse_json(files['original-r1'])
                    original_r2 = parse_json(files['original-r2'])
                    final_r1 = parse_json(files['final-r1'])
                    final_r2 = parse_json(files['final-r2'])
                    coverage = final_r1['qc_stats']['coverage'] * 2
                read_length = round(final_r1['qc_stats']['read_mean'])
                quality = final_r1['qc_stats']['qual_mean']
                contigs = assembly['total_contig']
                is_paired = True if end_type == 'paired-end' else False
                rank = get_rank(RANK_CUTOFF, coverage, quality, read_length, contigs, is_paired)

                print('\t'.join([str(a) for a in [rank, coverage, quality, read_length, contigs, is_paired]]))
                COUNTS[rank] += 1
                CATEGORIES[rank].append(directory.name)
                add_to_counts(assembly, rank)
                add_to_counts(annotation, rank)
                add_to_counts(merge_stats(final_r1, final_r2), rank)
                if rank == 'fail':
                    FAILED['failed-cutoff'].append(directory.name)
                    CATEGORIES['exclude'].append([directory.name, 'Failed to pass minimum cutoffs'])
                    COUNTS['exclude'] += 1
                else:
                    COUNTS['pass'] += 1


    #print("rank\tcount\tcoverage\tlength\tquality\tcontigs\tn50")
    for key, val in COUNTS_BY_RANK.items():
        coverage = int(mean(val['coverage'])) if val['coverage'] else 0
        read_length = int(mean(val['read_mean'])) if val['read_mean'] else 0
        quality = int(mean(val['qual_mean'])) if val['qual_mean'] else 0
        contigs = int(mean(val['total_contig'])) if val['total_contig'] else 0
        n50 = int(mean(val['n50_contig_length'])) if val['n50_contig_length'] else 0
        #print(f'{key}\t{COUNTS[key]}\t{coverage}\t{read_length}\t{quality}\t{contigs}\t{n50}')

    templateLoader = jinja2.FileSystemLoader(searchpath=TEMPLATE_DIR)
    templateEnv = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=TEMPLATE_DIR),
        autoescape=True
    )
    template = templateEnv.get_template('summary.j2')
    outputText = template.render(
        data={
            'counts': COUNTS,
            'failed': FAILED,
            'rank_counts': COUNTS_BY_RANK
        }
    )
    #print(outputText)
