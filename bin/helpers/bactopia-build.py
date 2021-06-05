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
import sys

VERSION = "1.7.1"
PROGRAM = "bactopia build"
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


def get_platform():
    from sys import platform
    if platform == "darwin":
        return 'mac'
    elif platform == "win32":
        # Windows is not supported
        print("Windows is not supported.", file=sys.stderr)
        sys.exit(1)
    return 'linux'


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


def check_needs_build(observed_md5, expected_md5, prefix, force=False, is_bactopia=False):
    """Check if a new environment needs to be built."""
    needs_build = False
    if os.path.exists(observed_md5) and not force:
        if check_md5sum(expected_md5, observed_md5):
            if not is_bactopia:
                logging.info(f'Existing env ({prefix}) found, skipping unless --force is used')
        else:
            needs_build = True
            logging.info(f'Existing env ({prefix}) is out of sync, it will be updated')                    
    else:
        needs_build = True
    return needs_build


def build_conda_env(env_file, prefix, max_retry=5, force=False, is_bactopia=False):
    """Build Conda env, with chance to retry."""
    force = '--force' if force else ''
    if is_bactopia:
        force = '--force'
    retry = 0
    allow_fail = False
    success = False
    while not success:
        result = execute(f'conda env create -f {env_file} --prefix {prefix} {force}', allow_fail=allow_fail)
        if not result:
            if retry > max_retry:
                allow_fail = True
            retry += 1
            logging.log(STDERR, "Error creating Conda environment, retrying after short sleep.")
            time.sleep(30 * retry)
        else:
            success = True
    return success


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture=False, stdout_file=None,
            stderr_file=None, allow_fail=False):
    """A simple wrapper around executor."""
    from executor import ExternalCommand, ExternalCommandFailed
    try:
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
        return True
    except ExternalCommandFailed as e:
        if allow_fail:
            print(e, file=sys.stderr)
            sys.exit(e.returncode)
        else:
            return None


if __name__ == '__main__':
    import argparse as ap
    import glob
    import sys
    import time

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
    parser.add_argument('--envname', metavar='STR', type=str,
                        help='Build Conda environment with the given name')
    parser.add_argument('--default', action='store_true',
                        help='Builds Conda environments to the default Bactopia location.')
    parser.add_argument('--max_retry', metavar='INT', type=int, default=5,
                        help='Maximum times to attempt creating Conda environment. (Default: 5)')           
    parser.add_argument('--force', action='store_true',
                        help='Force overwrite of existing Conda environments.')
    parser.add_argument('--is_bactopia', action='store_true',
                        help='This is an automated call by bactopia not a user')
    parser.add_argument('--include_tools', action='store_true',
                        help='Builds Conda environments for Bactopia tools as well.')
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
    ostype = get_platform()
    major, minor, patch = VERSION.split('.')
    CONTAINER_VERSION = f'{major}.{minor}.x'

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    # https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    env_path = f'{os.path.abspath(os.path.expanduser(args.conda_envs))}/{ostype}'
    install_path = os.path.abspath(os.path.expanduser(args.install_path))
    finish_file = f'{install_path}/envs-built-{CONTAINER_VERSION}.txt'
    if os.path.exists(finish_file):
        print(f'Found Conda environments in {install_path}, if a complete rebuild is needed please use --force')

    env_files = sorted(glob.glob(f'{env_path}/*.{args.ext}'))
    if env_files:
        for i, env_file in enumerate(env_files):
            envname = os.path.splitext(os.path.basename(env_file))[0]
            md5_file = env_file.replace('.yml', '.md5')
            prefix = f'{install_path}/{envname}-{CONTAINER_VERSION}'
            envbuilt_file = f'{install_path}/{envname}-{CONTAINER_VERSION}/env-built.txt'
            build = True
            if args.envname:
                if not args.envname == envname:
                    build = False
            
            if build:
                if check_needs_build(envbuilt_file, md5_file, prefix, force=args.force, is_bactopia=args.is_bactopia):
                    logging.info(f'Found {env_file} ({i+1} of {len(env_files)}), begin build to {prefix}')

                    built = build_conda_env(env_file, prefix, max_retry=args.max_retry, force=args.force, is_bactopia=args.is_bactopia)
                    if built:
                        execute(f'cp {md5_file} {envbuilt_file}')
        execute(f'touch {install_path}/envs-built-{CONTAINER_VERSION}.txt')
    else:
        logging.error(f'Unable to find Conda *.{args.ext} files in {env_path}, please verify')
        sys.exit(1)

    if args.include_tools:
        tool_path = os.path.abspath(args.conda_envs).replace('conda', 'tools')
        tools = sorted(glob.glob(f'{tool_path}/*/'))
        for i, tool in enumerate(tools):
            tool = os.path.basename(os.path.dirname(tool))
            if not tool.startswith('.'):
                env_file = f'{tool_path}/{tool}/environment-{ostype}.yml'
                if os.path.exists(env_file):
                    md5_file = f'{tool_path}/{tool}/environment-{ostype}.md5'
                    prefix = f'{install_path}/tools-{tool}-{CONTAINER_VERSION}'
                    envbuilt_file = f'{prefix}/env-built.txt'
                    force = '--force' if args.force else ''
                    build = True
                    if args.envname:
                        if not args.envname == tool:
                            build = False

                    if build:
                        if check_needs_build(envbuilt_file, md5_file, prefix, force=args.force, is_bactopia=args.is_bactopia):
                            logging.info(f'Found {env_file} ({i+1} of {len(tools)}), begin build to {prefix}')
                            built = build_conda_env(env_file, prefix, max_retry=args.max_retry, force=args.force, is_bactopia=args.is_bactopia)
                            if built:
                                execute(f'cp {md5_file} {envbuilt_file}')
