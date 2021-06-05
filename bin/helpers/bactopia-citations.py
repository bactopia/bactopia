#! /usr/bin/env python3
"""
usage: bactopia citations [-h] [--bactopia STR] [--version] STR

bactopia citations - Prints the citations of datasets and tools used by Bactopia

optional arguments:
  -h, --help      show this help message and exit
  --bactopia STR  Directory where Bactopia repository is stored.
  --version       show program's version number and exit
"""
import os
import sys

VERSION = "1.7.1"
PROGRAM = "bactopia citations"
DESCRIPTION = 'Prints the citations of datasets and tools used by Bactopia'


def validate_args(bactopia_repo):
    bactopia_citations = f'{bactopia_repo}/data/bactopia-datasets-software.txt'
    if not os.path.exists(bactopia_citations):
        print(f"cannot access '{bactopia_citations}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)
    else:
        citations = {}
        with open(bactopia_citations, 'rt') as citation_fh:
            for line in citation_fh:
                line.rstrip()
                if not line.startswith('name'):
                    name, ref_type, citation = line.split('\t')
                    if ref_type not in citations:
                        citations[ref_type] = []
                    citations[ref_type].append({'name':name, 'citation': citation})
        return citations

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
