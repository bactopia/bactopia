#! /usr/bin/env python3
"""
usage: bactopia tools summary [-h] [--gold_coverage FLOAT]
                              [--gold_quality INT] [--gold_read_length INT]
                              [--gold_contigs INT] [--silver_coverage FLOAT]
                              [--silver_quality INT]
                              [--silver_read_length INT]
                              [--silver_contigs INT] [--min_coverage FLOAT]
                              [--min_quality INT] [--min_read_length INT]
                              [--max_contigs INT] [--min_genome_size FLOAT]
                              [--max_genome_size FLOAT]
                              [--outdir OUTPUT_DIRECTORY] [--prefix STR]
                              [--depends] [--version] [--verbose] [--silent]
                              BACTOPIA_DIRECTORY

bactopia tools summary - Create a summary of Bactopia outputs

positional arguments:
  BACTOPIA_DIRECTORY    Directory containing Bactopia output.

optional arguments:
  -h, --help            show this help message and exit

Gold Cutoffs:
  --gold_coverage FLOAT
                        Minimum amount of coverage required for Gold status
                        (Default: 100x)
  --gold_quality INT    Minimum per-read mean quality score required for Gold
                        status (Default: Q30)
  --gold_read_length INT
                        Minimum mean read length required for Gold status
                        (Default: 95bp)
  --gold_contigs INT    Maximum contig count required for Gold status
                        (Default: 100)

Silver Cutoffs:
  --silver_coverage FLOAT
                        Minimum amount of coverage required for Silver status
                        (Default: 50x)
  --silver_quality INT  Minimum per-read mean quality score required for
                        Silver status (Default: Q20)
  --silver_read_length INT
                        Minimum mean read length required for Silver status
                        (Default: 75bp)
  --silver_contigs INT  Maximum contig count required for Silver status
                        (Default: 200)

Fail Cutoffs:
  --min_coverage FLOAT  Minimum amount of coverage required to pass (Default:
                        20x)
  --min_quality INT     Minimum per-read mean quality score required to pass
                        (Default: Q12)
  --min_read_length INT
                        Minimum mean read length required to pass (Default:
                        49bp)
  --max_contigs INT     Maximum contig count required to pass (Default: 500)
  --min_genome_size FLOAT
                        Minimum assembled genome size.
  --max_genome_size FLOAT
                        Maximum assembled genome size.

Helpers:
  --outdir OUTPUT_DIRECTORY
                        Directory to write output. (Default:
                        BACTOPIA_DIR/bactopia-tools)
  --prefix STR          Prefix to use for output files. (Default: bactopia)
  --depends             Verify dependencies are installed.
  --version             show program's version number and exit
  --verbose             Increase the verbosity of output.
  --silent              Only critical errors will be printed.

example usage:
  bactopia tools summary bactopiadir outdir
"""
import json
import logging
import os

from collections import defaultdict
from statistics import mean

PROGRAM = "bactopia tools summary"
VERSION = "1.7.1"
TEMPLATE_DIR = os.path.dirname(os.path.realpath(__file__)).replace('bin', 'templates')

STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")

IGNORE_LIST = ['.nextflow', '.nextflow.log', 'bactopia-info', 'work', 'bactopia-tools']
COUNTS = defaultdict(int)
FAILED = defaultdict(list)
CATEGORIES = defaultdict(list)
COUNTS_BY_RANK = {
    'total': defaultdict(list),
    'gold': defaultdict(list),
    'silver': defaultdict(list),
    'bronze': defaultdict(list),
    'exclude': defaultdict(list),
    'qc-failure': defaultdict(list)
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
        f"{path}/{name}/{name}-low-basepair-proportion-error.txt",
        f"{path}/{name}/{name}-different-read-count-error.txt",
        f"{path}/{name}/{name}-paired-end-error.txt",
        f"{path}/{name}/{name}-assembly-error.txt",
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
        return ["genome-size-error", "Poor estimate of genome size"]
    elif error.endswith("low-read-count-error.txt"):
        return ["low-read-count-error", "Low number of reads"]
    elif error.endswith("low-sequence-depth-error.txt"):
        return ["low-sequence-depth-error", "Low depth of sequencing"]
    elif  error.endswith("paired-end-error.txt"):
        return ["paired-end-error", "Paired-end reads were not in acceptable format"]
    elif  error.endswith("different-read-count-error.txt"):
        return ["different-read-count-error", "Paired-end read count mismatch"]
    elif  error.endswith("low-basepair-proportion-error.txt"):
        return ["low-basepair-proportion-error", "Paired-end basepair counts are out of accesptable proportions"]
    elif  error.endswith("assembly-error.txt"):
        return ["assembly-error", "Assembled size was not withing an acceptable range"]

    return ["unknown-error", "Unknown Error"]

def parse_genome_size(txt_file):
    """Pull out the genome size used for analysis."""
    with open(txt_file, 'rt') as txt_fh:
        return txt_fh.readline().rstrip()

def parse_json(json_file):
    """Return a dict of the input json_file"""
    with open(json_file) as json_fh:
        return json.load(json_fh)

def parse_annotation(txt_file, prefix='annotation'):
    """Parse Prokka summary text file."""
    results = {}
    with open(txt_file, 'rt') as txt_fh:
        for line in txt_fh:
            line = line.rstrip()
            key, val = line.split(":")
            results[f'{prefix}_{key}'] = val.lstrip()
    return results

def gather_stats(files, rank_cutoff):
    """Return a dictionary of combined stats."""
    stats = {}
    stats['estimated_genome_size'] = parse_genome_size(files['genome_size'])
    if 'original' in files:
        stats['is_paired'] = False
        stats.update(merge_qc_stats(parse_json(files['original']), None))
        stats.update(merge_qc_stats(parse_json(files['final']), None, prefix='final'))
    else:
        stats['is_paired'] = True
        stats.update(merge_qc_stats(parse_json(files['original-r1']), parse_json(files['original-r2'])))
        stats.update(merge_qc_stats(parse_json(files['final-r1']), parse_json(files['final-r2']), prefix='final'))
    stats.update(parse_json(files['assembly']))
    stats.update(parse_annotation(files['annotation']))
    rank, reason = get_rank(rank_cutoff, stats['final_coverage'], stats['final_qual_mean'], stats['final_read_mean'],
                            stats['total_contig'], stats['total_contig_length'], stats['is_paired'])
    stats['rank'] = rank
    stats['reason'] = reason
    
    return stats

def get_files(path, sample):
    "Return a list of files to read."
    files = None
    missing = []
    end_type = 'paired-end'

    if os.path.exists(f"{path}/{sample}/quality-control/{sample}.fastq.gz"):
        # Single End
        end_type = 'single-end'
        files = {
            'genome_size': f"{path}/{sample}/{sample}-genome-size.txt",
            'annotation': f"{path}/{sample}/annotation/{sample}.txt",
            'assembly': f"{path}/{sample}/assembly/{sample}.fna.json",
            'original': f"{path}/{sample}/quality-control/summary-original/{sample}-original.json",
            'final': f"{path}/{sample}/quality-control/summary-final/{sample}-final.json",
        }
    else:
        # Paired End
        files = {
            'genome_size': f"{path}/{sample}/{sample}-genome-size.txt",
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

def get_rank(cutoff, coverage, quality, length, contigs, genome_size, is_paired):
    """Return the rank (gold, silver, bronze, fail) of the sample."""
    rank = None
    reason = []
    coverage = float(f'{coverage:.2f}')
    gold = cutoff['gold']
    silver = cutoff['silver']
    bronze = cutoff['bronze']

    if coverage >= gold['coverage'] and quality >= gold['quality'] and length >= gold['length'] and contigs <= gold['contigs'] and is_paired:
        reason.append('passed all cutoffs')
        rank = 'gold'
    elif coverage >= silver['coverage'] and quality >= silver['quality'] and length >= silver['length'] and contigs <= silver['contigs'] and is_paired:
        if coverage < gold['coverage']:
            reason.append(f"Low coverage ({coverage:.2f}x, expect >= {gold['coverage']}x)")
        if quality < gold['quality']:
            reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{gold['quality']})")
        if length < gold['length']:
            reason.append(f"Short read length ({length:.2f}bp, expect >= {gold['length']} bp)")
        if contigs > gold['contigs']:
            reason.append(f"Too many contigs ({contigs}, expect <= {gold['contigs']})")
        rank = 'silver'
    elif coverage >= bronze['coverage'] and quality >= bronze['quality'] and length >= bronze['length'] and contigs <= bronze['contigs']:
        if coverage < silver['coverage']:
            reason.append(f"Low coverage ({coverage:.2f}x, expect >= {silver['coverage']}x)")
        if quality < silver['quality']:
            reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{silver['quality']})")
        if length < silver['length']:
            reason.append(f"Short read length ({length:.2f}bp, expect >= {silver['length']} bp)")
        if contigs > silver['contigs']:
            reason.append(f"Too many contigs ({contigs}, expect <= {silver['contigs']})")
        if not is_paired:
            reason.append(f"Single-end reads")
        rank = 'bronze'


    if not rank:
        rank = 'exclude'

    if coverage < bronze['coverage']:
        reason.append(f"Low coverage ({coverage:.2f}x, expect >= {bronze['coverage']}x)")
    if quality < bronze['quality']:
        reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{bronze['quality']})")
    if length < bronze['length']:
        reason.append(f"Short read length ({length:.2f}bp, expect >= {bronze['length']} bp)")
    if contigs > bronze['contigs']:
        reason.append(f"Too many contigs ({contigs}, expect <= {bronze['contigs']})")

    if cutoff['min-genome-size']:
        if genome_size < cutoff['min-genome-size']:
            reason.append(f"Assembled genome size is too small ({genome_size} bp, expect <= {cutoff['min-genome-size']})")

    if cutoff['max-genome-size']:
        if genome_size < cutoff['max-genome-size']:
            reason.append(f"Assembled genome size is too large ({genome_size} bp, expect <= {cutoff['max-genome-size']})")

    reason = ";".join(sorted(reason))
    return [rank, reason]

def merge_qc_stats(r1, r2, prefix='original'):
    merged = {}
    for key in r1['qc_stats']:
        prefixed_key = f'{prefix}_{key}'
        if key in ['total_bp', 'coverage', 'read_total']:
            merged[prefixed_key] = r1['qc_stats'][key] + r2['qc_stats'][key] if r2 else r1['qc_stats'][key]
        else:
            merged[prefixed_key] = mean([r1['qc_stats'][key], r2['qc_stats'][key]]) if r2 else r1['qc_stats'][key]
    return merged

def add_to_counts(dictionary):
    """Append values to rank counts."""
    rank = dictionary['rank']
    for key, val in dictionary.items():
        if key != 'rank':
            COUNTS_BY_RANK[rank][key].append(val)
            COUNTS_BY_RANK['total'][key].append(val)

def generate_html_report():
    """Generate a HTML report using Jinja2."""
    import jinja2
    templateEnv = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=TEMPLATE_DIR),
        autoescape=True
    )
    template = templateEnv.get_template('summary.j2')
    return template.render(
        data={
            'counts': COUNTS,
            'failed': FAILED,
            'rank_counts': COUNTS_BY_RANK
        }
    )

def generate_txt_report(output_file):
    return None

def print_failed(failed):
    """Print failed counts."""
    lines = []
    for key, val in sorted(failed.items()):
        if key != 'failed-cutoff':
            pretty_key = key.replace('-', ' ').title()
            lines.append(f'        {pretty_key}: {len(val)}')
    return '\n'.join(lines)


def print_cutoffs(cutoffs):
    """Print cutoff counts."""
    lines = []
    for key, val in sorted(cutoffs.items()):
        lines.append(f'        {key}: {val}')
    return '\n'.join(lines)

if __name__ == '__main__':
    import argparse as ap
    import csv
    import sys
    import textwrap
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Create a summary of Bactopia outputs',
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
        '--min_genome_size', metavar="FLOAT", type=float,
        help='Minimum assembled genome size.'
    )

    group3.add_argument(
        '--max_genome_size', metavar="FLOAT", type=float,
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
        },
        'min-genome-size': args.min_genome_size,
        'max-genome-size': args.max_genome_size
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
                        error_msg = []
                        for error in has_error:
                            error_type, error_name = parse_error(error)
                            error_msg.append(error_name)
                            COUNTS[error_type] += 1
                            FAILED[error_type].append(directory.name)
                        COUNTS['total-excluded'] += 1
                        COUNTS['qc-failure'] += 1
                        CATEGORIES['failed'].append([directory.name, f"Not processed, reason: {';'.join(error_msg)}"])
                        samples[directory.name] = {'has_error': True, 'missing': []}
                    else:
                        logging.debug(f"{directory.name} found ")
                        COUNTS['processed'] += 1
                        CATEGORIES['processed'].append(directory.name)
                        samples[directory.name] = get_files(args.bactopia, directory.name)

    sample_stats = {}
    stat_fields = []
    for sample, val in samples.items():
        if not val['has_error']:
            end_type = val['end_type']
            files = val['files']
            missing = val['missing']
            COUNTS[end_type] += 1
            if missing:
                COUNTS['missing'] += 1
                CATEGORIES['failed'].append([sample, 'Missing expected files'])
                logging.debug(f"{sample} missing files ")
                logging.debug(f"{missing}")
            else:
                COUNTS['found'] += 1
                stats = gather_stats(files, RANK_CUTOFF)
                add_to_counts(stats)
                COUNTS[stats['rank']] += 1
                CATEGORIES[stats['rank']].append(sample)
                if stats['rank'] == 'exclude':
                    FAILED['failed-cutoff'].append(sample)
                    CATEGORIES['failed'].append([sample, f'Failed to pass minimum cutoffs, reason: {stats["reason"]}'])
                    COUNTS['total-excluded'] += 1
                else:
                    COUNTS['pass'] += 1
                sample_stats[sample] = stats
                for key in stats:
                    if key not in stat_fields:
                        stat_fields.append(key)

    for key, val in COUNTS_BY_RANK.items():
        coverage = int(mean(val['coverage'])) if val['coverage'] else 0
        read_length = int(mean(val['read_mean'])) if val['read_mean'] else 0
        quality = int(mean(val['qual_mean'])) if val['qual_mean'] else 0
        contigs = int(mean(val['total_contig'])) if val['total_contig'] else 0
        n50 = int(mean(val['n50_contig_length'])) if val['n50_contig_length'] else 0

    # Write outputs
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # HTML report
    """
    html_report = f'{outdir}/{args.prefix}-report.html'
    with open(html_report, 'w') as html_fh:
        html_fh.write(generate_html_report())
    """

    # Tab-delimited report
    txt_report = f'{outdir}/{args.prefix}-report.txt'
    with open(txt_report, 'w') as txt_fh:
        outputs = []
        for sample, stats in sorted(sample_stats.items()):
            output = {'sample': sample}
            for field in stat_fields:
                if field in stats:
                    if isinstance(stats[field], float):
                        output[field] = f"{stats[field]:.3f}"
                    else:
                        output[field] = stats[field]
                else:
                    output[field] = ""
            outputs.append(output)

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
            if name in sample_stats:
                reasons = reason.split(':')[1].split(';')
                cutoffs = []
                for r in reasons:
                    cutoffs.append(r.split('(')[0].strip().title())
                cutoff_counts[';'.join(sorted(cutoffs))] += 1
                exclude_fh.write(f'{name}\texclude\t{reason}\n')
            else:
                exclude_fh.write(f'{name}\tqc-fail\t{reason}\n')

    # Screen report
    with open(f'{outdir}/{args.prefix}-summary.txt', 'w') as summary_fh:
        summary_fh.write("Bactopia Summary Report\n")
        summary_fh.write(textwrap.dedent(f'''
            Total Samples: {len(samples)}
            
            Passed: {COUNTS["pass"]}
                Gold: {COUNTS["gold"]}
                Silver: {COUNTS["silver"]}
                Bronze: {COUNTS["bronze"]}

            Excluded: {COUNTS["total-excluded"]}
                Failed Cutoff: {COUNTS["exclude"]}\n'''))
        summary_fh.write(print_cutoffs(cutoff_counts))
        summary_fh.write('\n')
        summary_fh.write(f'    QC Failure: {COUNTS["qc-failure"]}\n')
        summary_fh.write(print_failed(FAILED))
        summary_fh.write(textwrap.dedent(f'''\n
            Reports:
                Full Report (txt): {txt_report}
                Exclusion: {exclusion_report}

            Rank Cutoffs:\n'''))
        summary_fh.write(json.dumps(RANK_CUTOFF, indent=4))
        summary_fh.write('\n')
