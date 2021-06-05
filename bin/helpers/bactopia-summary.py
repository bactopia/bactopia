#! /usr/bin/env python3
"""


"""
import os
import sys

VERSION = "1.7.1"
PROGRAM = "bactopia summary"
DESCRIPTION = 'Create a summary report for samples'

def get_output_files():
    """Return a dictionary opf output files to include in the summary."""
    """
    ${SAMPLE_NAME}/
    ├── annotation
    ├── antimicrobial_resistance
    ├── ariba
    ├── assembly
    ├── blast
    ├── kmers
    ├── logs
    ├── mapping
    ├── minmers
    ├── mlst
    ├── quality-control
    ├── variants
    └── ${SAMPLE_NAME}-genome-size.txt
    """


if __name__ == '__main__':
    import argparse as ap
    import textwrap

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - {DESCRIPTION}'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument('--bactopia', metavar="STR", type=str,
                        help='Directory where Bactopia repository is stored.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    citations = validate_args(args.bactopia)

    for ref_type, entries in sorted(citations.items()):
        print(f'# {ref_type} potentially used by Bactopia')
        print('# ----------')
        for entry in entries:
            print(f'## {entry["name"]}')
            print(textwrap.fill(entry['citation'], width=100))
            print()
