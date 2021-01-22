#! /usr/bin/env python3
"""
"""
PROGRAM = "split-coverages"
VERSION = "1.6.0"

if __name__ == '__main__':
    import argparse as ap
    import os
    import sys
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Split a genomeCoverageBed output into separate files based on FASTA entry'
        )
    )

    parser.add_argument(
        'mapping', metavar="FILE", type=str,
        help='Tab-delimited file used to map entry names to original fasta file.'
    )
    parser.add_argument(
        'coverage', metavar="FILE", type=str,
        help='genomeCoverageBed output file'
    )
    parser.add_argument(
        '--outdir', metavar="STR", type=str, default='coverages',
        help='Directory to output split coverages into. (Default: coverages)'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    mappings = {}
    with open(args.mapping, 'rt') as mapping_fh:
        for line in mapping_fh:
            fasta, entry = line.rstrip().split('\t')
            mappings[entry] = fasta

    coverages = {}
    with open(args.coverage, 'rt') as coverage_fh:
        for line in coverage_fh:
            entry, position, depth = line.rstrip().split('\t')
            if mappings[entry] not in coverages:
                coverages[mappings[entry]] = {}

            if entry not in coverages[mappings[entry]]:
                coverages[mappings[entry]][entry] = []

            coverages[mappings[entry]][entry].append(depth)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    for fasta in coverages:
        with open(f'{args.outdir}/{fasta}-coverage.txt', 'wt') as coverage_out:
            total_entries = len(coverages[fasta])
            coverage_out.write(f'##total={total_entries}\n')
            for entry, depths in coverages[fasta].items():
                coverage_out.write(f'##contig=<ID={entry},length={len(depths)}>\n')
                for depth in depths:
                    coverage_out.write(f'{depth}\n')
 