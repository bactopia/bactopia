#! /usr/bin/env python3
"""
usage: cleanup-coverage [-h] [--mincov INT] [--version] COVERAGE

cleanup-coverage - Reduce redundancy in per-base coverage.

positional arguments:
  COVERAGE      Output from genomeBedCoverage

optional arguments:
  -h, --help    show this help message and exit
  --version     show program's version number and exit
"""
PROGRAM = "cleanup-coverage"
VERSION = "1.7.1"
import sys

def read_coverage(coverage):
    """Read the per-base coverage input."""
    import re
    accession = None
    length = None
    first_line = True
    coverages = {}
    with open(coverage, 'rt') as coverage_fh:
        for line in coverage_fh:
            line = line.rstrip()
            if line.startswith('##'):
                # ##contig=<ID=NZ_CP020108,length=5407749>
                contig = re.search(r'contig=<ID=(.*),length=([0-9]+)>', line)
                if contig:
                    accession = contig.group(1)
                    length = contig.group(2)
                    coverages[accession] = {'length':int(length), 'positions': []}
                else:
                    print(f'{line} is an unexpected format.', file=sys.stderr)
                    sys.exit(1)
            else:
                accession, position, coverage = line.split('\t')
                coverages[accession]['positions'].append(int(coverage))

    for accession, vals in coverages.items():
        if len(vals['positions']) != vals['length']:
            print(f'Observed bases ({len(vals["positions"])} in {accession} not expected length ({vals["length"]}).', file=sys.stderr)
            sys.exit(1)

    return coverages

if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Snippy consensus (subs) with coverage masking.'
        )
    )
    parser.add_argument('coverage', metavar="COVERAGE", type=str,
                        help='Directory where BLAST databases are stored')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    coverages = read_coverage(args.coverage)
    for accession, vals in coverages.items():
        print(f'##contig=<ID={accession},length={vals["length"]}>')
        for cov in vals['positions']:
            print(cov)
