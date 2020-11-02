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

VERSION = "1.4.11"
PROGRAM = "bactopia tools"
DESCRIPTION = 'A suite of comparative analyses for Bactopia outputs'
AVAILABLE_TOOLS = {
    'eggnog': {'info': 'Functional annotation using orthologous groups', 'mac': True},
    'fastani': {'info': 'Pairwise average nucleotide identity', 'mac': True},
    'gtdb': {'info': 'Identify marker genes and assign taxonomic classifications', 'mac': False},
    'ismapper': {'info': 'Identify positions of insertion sites', 'mac': True},
    'mashtree': {'info': 'Trees based on Mash distances', 'mac': True},
    'pirate': {'info': 'Pan-genome with optional core-genome tree', 'mac': True},
    'phyloflash': {'info': '16s assembly, alignment and tree', 'mac': True},
    'roary': {'info': 'Pan-genome with optional core-genome tree', 'mac': True},
    'summary': {'info': 'A report summarizing Bactopia project', 'mac': True},
}

def get_platform():
    from sys import platform
    if platform == "darwin":
        return 'mac'
    elif platform == "win32":
        # Windows is not supported
        print("Windows is not supported.", file=sys.stderr)
        sys.exit(1)
    return 'linux'

def print_available_tools():
    """Print the available Bactopia Tools."""
    print(f"{PROGRAM} (v{VERSION}) - {DESCRIPTION}")
    print("")
    print(available_tools())

def available_tools():
    """Return a string of available tools."""
    usage = ['Available Tools:']
    for k,v in sorted(AVAILABLE_TOOLS.items()):
        usage.append(f'  {k: <12}{v["info"]}')
    return '\n'.join(usage)

def validate_args(tool, bactopia_repo):
    import os
    platform = get_platform()

    if tool not in AVAILABLE_TOOLS:
        print(f'"{tool}" is not available.\n', file=sys.stderr)
        print(available_tools(), file=sys.stderr)
        sys.exit(1)
    elif platform == 'mac' and not AVAILABLE_TOOLS[tool]['mac']:
        print(f'"{tool}" is not available on Mac OSX.\n', file=sys.stderr)
        sys.exit(1)


    tool_nf = f'{bactopia_repo}/tools/{tool}/main.nf'
    if not os.path.exists(tool_nf):
        print(f"cannot access '{tool_nf}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)

    tool_env = f'{bactopia_repo}/tools/{tool}/environment-linux.yml'
    if platform == 'mac':
        tool_env = f'{bactopia_repo}/tools/{tool}/environment-osx.yml'

    return f"{tool_nf} --conda_yaml {tool_env}" 

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
    print(validate_args(args.tool, args.bactopia))
