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
import logging
import os
import sys

STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")

VERSION = "1.7.1"
PROGRAM = "bactopia tools"
DESCRIPTION = 'A suite of comparative analyses for Bactopia outputs'
AVAILABLE_TOOLS = {
    'eggnog': {'info': 'Functional annotation using orthologous groups', 'mac': True},
    'fastani': {'info': 'Pairwise average nucleotide identity', 'mac': True},
    'gtdb': {'info': 'Identify marker genes and assign taxonomic classifications', 'mac': False},
    'hicap': {'info': 'in silico typing of the H. influenzae cap locus', 'mac': True},
    'ismapper': {'info': 'Identify positions of insertion sites', 'mac': True},
    'mashtree': {'info': 'Trees based on Mash distances', 'mac': True},
    'pirate': {'info': 'Pan-genome with optional core-genome tree', 'mac': True},
    'phyloflash': {'info': '16s assembly, alignment and tree', 'mac': True},
    'roary': {'info': 'Pan-genome with optional core-genome tree', 'mac': True},
    'staph-typer': {'info': 'Tools for typing Staphylococcus aureus.', 'mac': True},
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


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


def check_md5sum(expected_md5, current_md5):
    """Compare the two md5 files to see if a rebuild is needed."""
    expected = None
    current = None
    with open(expected_md5, 'r') as f:
        expected = f.readline().rstrip()

    with open(current_md5, 'r') as f:
        current = f.readline().rstrip()

    return expected == current


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture=False, stdout_file=None,
            stderr_file=None):
    """A simple wrapper around executor."""
    from executor import ExternalCommand
    command = ExternalCommand(
        cmd, directory=directory, capture=True, capture_stderr=True,
        stdout_file=stdout_file, stderr_file=stderr_file
    )

    command.start()
    if get_log_level() == 'DEBUG':
        logging.log(STDOUT, command.decoded_stdout)
        logging.log(STDERR, command.decoded_stderr)

    if capture:
        return command.decoded_stdout


def validate_args(tool, bactopia_repo, skip_conda=False, force_rebuild=False):
    import os
    platform = get_platform()

    if tool not in AVAILABLE_TOOLS:
        print(f'"{tool}" is not available.\n', file=sys.stderr)
        print(available_tools(), file=sys.stderr)
        sys.exit(1)
    elif platform == 'mac' and not AVAILABLE_TOOLS[tool]['mac']:
        print(f'"{tool}" is not available on Mac OSX.\n', file=sys.stderr)
        sys.exit()
    tool_nf = f'{bactopia_repo}/tools/{tool}/main.nf'
    if not os.path.exists(tool_nf):
        print(f"cannot access '{tool_nf}': No such file or directory\n",
              file=sys.stderr)
        print("Please make sure the correct path to Bactopia's repo is given.",
              file=sys.stderr)
        sys.exit(1)

    conda_prefix = f'{bactopia_repo}/tools/{tool}/environment-linux'
    if platform == 'mac':
        conda_prefix = f'{bactopia_repo}/tools/{tool}/environment-osx'

    if skip_conda:
        return f"{tool_nf}"
    else:
        # Check if conda env exists
        major, minor, patch = VERSION.split('.')
        CONTAINER_VERSION = f'{major}.{minor}.x'
        needs_build = False
        condadir = f'{bactopia_repo}/conda/envs/tools-{tool}-{CONTAINER_VERSION}'
        envbuilt_file = f'{condadir}/env-built.txt'
        if os.path.exists(envbuilt_file) and not force_rebuild:
            build_is_current = check_md5sum(f'{conda_prefix}.md5', envbuilt_file)
            if build_is_current:
                logging.info(f'Existing env ({condadir}) found, skipping unless --force_rebuild is used')
            else:
                needs_build = True
                force_rebuild = True
                logging.info(f'Existing env ({condadir}) is out of sync, it will be updated')
        else:
            needs_build = True

        if needs_build:
            logging.info(f'Found {conda_prefix}.yml, begin build to {condadir}')
            force = '--force' if force_rebuild else ''
            execute(f'conda env create -f {conda_prefix}.yml --prefix {condadir} {force}')
            execute(f'cp {conda_prefix}.md5 {envbuilt_file}')

        return f"{tool_nf} --condadir {condadir}"


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
    parser.add_argument('--force_rebuild', action='store_true',
                        help='Force overwrite of existing Conda environments.')
    parser.add_argument('--skip_conda', action='store_true',
                        help='Skip all things conda related.')
    parser.add_argument('--verbose', action='store_true',
                        help='Print debug related text.')
    parser.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        print_available_tools()
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))
    print(validate_args(
        args.tool, args.bactopia,
        skip_conda=args.skip_conda, 
        force_rebuild=args.force_rebuild
    ))
