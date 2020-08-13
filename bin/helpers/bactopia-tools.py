#! /usr/bin/env python3
"""
usage: bactopia tools [-h] [--bactopia STR] [--version] STR

bactopia tools - A suite of comparative analyses for Bactopia outputs

positional arguments:
  STR             Name of the Bactopia tool to execute.

optional arguments:
  -h, --help      show this help message and exit
  --bactopia STR  Directory where Bactopia repository is stored.
  --version       show program's version number and exit
"""
import sys

VERSION = "1.4.5"
PROGRAM = "bactopia tools"
DESCRIPTION = 'A suite of comparative analyses for Bactopia outputs'
AVAILABLE_TOOLS = {
    'eggnog': 'Functional annotation using orthologous groups',
    'fastani': 'Pairwise average nucleotide identity',
    'gtdb': 'Identify marker genes and assign taxonomic classifications',
    'ismapper': 'Identify positions of insertion sites',
    'mashtree': 'Trees based on Mash distances',
    'pirate': 'Pan-genome with optional core-genome tree',
    'phyloflash': '16s assembly, alignment and tree',
    'roary': 'Pan-genome with optional core-genome tree',
    'summary': 'A report summarizing Bactopia project',
}


def print_available_tools():
    """Print the available Bactopia Tools."""
    print(f"{PROGRAM} (v{VERSION}) - {DESCRIPTION}")
    print("")
    print(available_tools())

def available_tools():
    """Return a string of available tools."""
    usage = ['Available Tools:']
    for k,v in sorted(AVAILABLE_TOOLS.items()):
        usage.append(f'  {k: <12}{v}')
    return '\n'.join(usage)

def validate_args(tool, bactopia_repo):
    import os

    if args.tool not in AVAILABLE_TOOLS:
        print(f'"{args.tool}" is not available.\n', file=sys.stderr)
        print(available_tools(), file=sys.stderr)
        sys.exit(1)

    tool_nf = f'{bactopia_repo}/tools/{tool}/main.nf'
    if not os.path.exists(tool_nf):
        print(f"cannot access '{tool_nf}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)

    return tool_nf

if __name__ == '__main__':
    import argparse as ap
    import textwrap

    parser = ap.ArgumentParser(
        prog='bactopia tools',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - {DESCRIPTION}'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=available_tools()
    )
    parser.add_argument('tool', metavar="STR", type=str,
                        help='Name of the Bactopia tool to execute.')
    parser.add_argument('--bactopia', metavar="STR", type=str,
                        help='Directory where Bactopia repository is stored.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        print_available_tools()
        sys.exit(0)

    args = parser.parse_args()
    tool_nf = validate_args(args.tool, args.bactopia)
    print(tool_nf)
