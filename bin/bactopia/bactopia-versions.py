#! /usr/bin/env python3
"""
usage: bactopia versions [-h] [--bactopia STR] [--version] STR

bactopia versions - Prints the version of tools used by Bactopia

optional arguments:
  -h, --help      show this help message and exit
  --bactopia STR  Directory where Bactopia repository is stored.
  --version       show program's version number and exit
"""

import os
import sys

VERSION = "1.7.1"
PROGRAM = "bactopia versions"
DESCRIPTION = 'Prints the version of tools used by Bactopia'


def get_platform():
    from sys import platform
    if platform == "darwin":
        return 'mac'
    elif platform == "win32":
        # Windows is not supported
        print("Windows is not supported.", file=sys.stderr)
        sys.exit(1)
    return 'linux'


def validate_args(bactopia_repo):
    import json 

    bactopia_json = f'{bactopia_repo}/conda/bactopia-programs.json'
    if not os.path.exists(bactopia_json):
        print(f"cannot access '{bactopia_json}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)
    else:
        with open(bactopia_json, 'rt') as json_fh:
            return json.load(json_fh)


def read_yaml(yaml):
    versions = {}
    with open(yaml, 'rt') as yaml_fh:
        for line in yaml_fh:
            line = line.strip()
            if '=' in line:
                program, version = line.replace('- ', '').split('=')[0:2]
                versions[program] = version
    return versions


if __name__ == '__main__':
    import argparse as ap

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
    ostype = get_platform()
    tools = validate_args(args.bactopia)

    conda_dir = f'{args.bactopia}/conda/{ostype}'
    yamls = [f'{f.name}' for f in os.scandir(conda_dir) if f.name.endswith('.yml')]
    versions = {}
    for yaml in yamls:
        versions[yaml] = read_yaml(f'{conda_dir}/{yaml}')

    final_versions = {}
    for tool, info in sorted(tools.items()):
        yaml = info['conda']['yaml']
        if yaml not in versions:
            if yaml.startswith("tools"):
                versions[yaml] = read_yaml(f'{args.bactopia}/{yaml}')
            else:
                versions[yaml] = read_yaml(f'{conda_dir}/{yaml}')

        if ostype == 'mac' and info['conda']['name'] == 'checkm-genome':
            continue
        else:
            final_versions[tool.lower()] = {
                'name': tool,
                'version': versions[yaml][info['conda']['name']],
                'description': info['description'],
                'link': info['link']
            }

    print(f'name\tversion\tdescription\tlink')
    for tool, cols in sorted(final_versions.items()):
        print(f'{cols["name"]}\t{cols["version"]}\t{cols["description"]}\t{cols["link"]}')
