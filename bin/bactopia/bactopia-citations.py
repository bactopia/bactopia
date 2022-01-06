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

VERSION = "2.0.1"
PROGRAM = "bactopia citations"
DESCRIPTION = 'Prints the citations of datasets and tools used by Bactopia'


def validate_args(bactopia_repo):
    import yaml
    bactopia_citations = f'{bactopia_repo}/citations.yml'
    if not os.path.exists(bactopia_citations):
        print(f"Cannot access '{bactopia_citations}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)
    else:
        citations = {}
        module_citations = {}
        counts = {}
        with open(bactopia_citations, "rt") as citations_fh:
            citations = yaml.safe_load(citations_fh)
            for group, refs in citations.items():
                for ref, vals in refs.items():
                    module_citations[ref] = vals
        return [citations, module_citations]

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
    parser.add_argument('--name', metavar="STR", type=str,
                        help='Only print citation matching a given name.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    citations, module_citations = validate_args(args.bactopia)

    if args.name:
        if args.name in module_citations:
            print(module_citations[args.name]['cite'].rstrip())
        else:
            print(f'"{args.name}" does not match available citations', file=sys.stderr)
    else:
        for group, refs in citations.items():
            print(f'# {group} potentially used by Bactopia')
            print('# ----------')
            for ref, vals in refs.items():
                print(f'## {vals["name"]}')
                print(textwrap.fill(vals['cite'], width=100))
