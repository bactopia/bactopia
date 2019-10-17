#! /usr/bin/env python3
"""
usage: bactopia build [-h] [-e STR] [--force] [--verbose] [--silent]
                      [--version]
                      STR STR

bactopia build - Build Conda environments for use by Bactopia

positional arguments:
  STR                Directory containing Conda environment files to build.
  STR                Directory to install Conda environments to.

optional arguments:
  -h, --help         show this help message and exit
  -e STR, --ext STR  Extension of the Conda environment files. Default: .yml
  --force            Force overwrite of existing Conda environments.
  --verbose          Print debug related text.
  --silent           Only critical errors will be printed.
  --version          show program's version number and exit
"""
import logging
import os

VERSION = "1.2.0"
PROGRAM = "bactopia build"
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


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

if __name__ == '__main__':
    import argparse as ap
    import glob
    import sys

    parser = ap.ArgumentParser(
        prog='bactopia build',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Build Conda environments for use by Bactopia'
        )
    )

    parser.add_argument('conda_envs', metavar="STR", type=str,
                        help='Directory containing Conda environment files to build.')

    parser.add_argument('install_path', metavar="STR", type=str,
                        help='Directory to install Conda environments to.')
    parser.add_argument(
        '-e', '--ext', metavar='STR', type=str,
        default="yml",
        help='Extension of the Conda environment files. Default: .yml'
    )
    parser.add_argument('--force', action='store_true',
                        help='Force overwrite of existing Conda environments.')
    parser.add_argument('--verbose', action='store_true',
                        help='Print debug related text.')
    parser.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    # https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    env_path = os.path.abspath(args.conda_envs)
    install_path = os.path.abspath(args.install_path)
    finish_file = f'{install_path}/envs-build-{VERSION}.txt'
    if not args.force and os.path.exists(finish_file):
        logging.error(
            f'Conda envs are already built in {install_path}, will not rebuild without --force'
        )
        sys.exit(1)

    env_files = sorted(glob.glob(f'{env_path}/*.{args.ext}'))
    if env_files:
        for i, env_file in enumerate(env_files):
            envname = os.path.splitext(os.path.basename(env_file))[0]
            prefix = f'{install_path}/{envname}-{VERSION}'
            logging.info(f'Found {env_file} ({i+1} or {len(env_files)}), begin build to {prefix}')
            force = '--force' if args.force else ''
            execute(f'conda env create -f {env_file} --prefix {prefix} {force}')
        execute(f'touch {install_path}/envs-build-{VERSION}.txt')
    else:
        logging.error(
            f'Unable to find Conda *.{args.ext} files in {env_path}, please verify'
        )
        sys.exit(1)
