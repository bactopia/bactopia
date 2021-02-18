#! /usr/bin/env python3
"""
usage: bactopia pull [-h] [--envname STR] [--singularity_cache STR]
                     [--registry STR] [--max_retry INT] [--include_tools]
                     [--default] [--is_bactopia] [--force] [--verbose]
                     [--silent] [--version]
                     STR

bactopia pull - Build Singularity images used by Bactopia

positional arguments:
  STR                   Directory containing Dockerfiles.

optional arguments:
  -h, --help            show this help message and exit
  --envname STR         Build Singularity images with the given name
  --singularity_cache STR
                        Directory where Singularity images will be stored.
  --registry STR        Docker registry to pull containers from
  --max_retry INT       Maximum times to attempt creating Conda environment.
                        (Default: 5)
  --include_tools       Singularity images for Bactopia Tools will also be
                        built.
  --default             Builds Singularity images to the default Bactopia
                        location.
  --is_bactopia         This is an automated call by bactopia not a user
  --force               Force overwrite of existing Conda environments.
  --verbose             Print debug related text.
  --silent              Only critical errors will be printed.
  --version             show program's version number and exit
"""
import logging
import os
import sys

VERSION = "1.6.0"
PROGRAM = "bactopia pull"
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
            logging.log(STDERR, e)
            sys.exit(e.returncode)
        else:
            return None


def get_docker_prefix(registry):
    """Return the proper prefix based on registry."""
    if registry == "quay":
        return 'quay.io'
    elif registry == "github":
        return 'ghcr.io'
    else:
        return ''


def check_needs_build(image, force=False, is_bactopia=False):
    """Check if a new image needs to be built."""
    if os.path.exists(image) and not force:
        if not is_bactopia:
            logging.info(f'Existing image ({image}) found, skipping unless --force is used')
        return False
    return True


def build_singularity_image(image, docker, max_retry=5, force=False, is_bactopia=False):
    """Build Conda env, with chance to retry."""
    force = '--force' if force else ''
    if is_bactopia:
        force = '--force'
    retry = 0
    allow_fail = False
    success = False
    while not success:
        result = execute(f'singularity build {force} {image} {docker}', allow_fail=allow_fail)
        if not result:
            if retry > max_retry:
                allow_fail = True
            retry += 1
            logging.log(STDERR, "Error creating image, retrying after short sleep.")
            time.sleep(30 * retry)
        else:
            success = True
    return success


if __name__ == '__main__':
    import argparse as ap
    import glob
    import sys
    import time
    from pathlib import Path

    parser = ap.ArgumentParser(
        prog='bactopia pull',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Build Singularity images used by Bactopia'
        )
    )

    parser.add_argument('docker', metavar="STR", type=str,
                        help='Directory containing Dockerfiles.')
    parser.add_argument('--envname', metavar='STR', type=str,
                        help='Build Singularity images with the given name')
    parser.add_argument('--singularity_cache', metavar='STR', type=str, default="~/.bactopia/singularity",
                        help='Directory where Singularity images will be stored.')
    parser.add_argument('--registry', metavar='STR', type=str, default="dockerhub",
                        help='Docker registry to pull containers from')
    parser.add_argument('--max_retry', metavar='INT', type=int, default=5,
                        help='Maximum times to attempt creating Conda environment. (Default: 5)')
    parser.add_argument('--include_tools', action='store_true',
                        help='Singularity images for Bactopia Tools will also be built.')
    parser.add_argument('--default', action='store_true',
                        help='Builds Singularity images to the default Bactopia location.')
    parser.add_argument('--is_bactopia', action='store_true',
                        help='This is an automated call by bactopia not a user')
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

    args, unknown = parser.parse_known_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    env_path = os.path.abspath(os.path.expanduser(args.docker))
    install_path = os.path.abspath(os.path.expanduser(args.singularity_cache))
    finish_file = f'{install_path}/{args.registry}-images-built-{VERSION}.txt'
    if os.path.exists(finish_file):
        print(f'Found Singularity images in {install_path}, if a complete rebuild is needed please use --force')
    
    if not os.path.exists(install_path):
        logging.info(f'Creating {install_path} to save images to')
        execute(f'mkdir -p {install_path}')

    registry = get_docker_prefix(args.registry)
    docker_prefix = f'docker://{registry}/bactopia' if registry else f'docker://bactopia'
    env_files = sorted(glob.glob(f'{env_path}/linux/*.yml'))
    if env_files:
        for i, env_file in enumerate(env_files):
            envname = os.path.basename(env_file).replace(".yml", "")
            img_name = f"{install_path}/{registry}-bactopia-{envname}-{VERSION}.img" if registry else f"{install_path}/bactopia-{envname}-{VERSION}.img"
            pull_name = f"{docker_prefix}/{envname}:{VERSION}"
            build = True
            if args.envname:
                if not args.envname == envname:
                    build = False
                    
            if build:
                if check_needs_build(img_name, force=args.force, is_bactopia=args.is_bactopia):
                    logging.info(f'Found {envname} ({i+1} of {len(env_files)}), begin build to {img_name}')

                    build_singularity_image(img_name, pull_name, max_retry=args.max_retry, force=args.force,
                                            is_bactopia=args.is_bactopia)
        execute(f'touch {finish_file}')
    else:
        logging.error(f'Unable to find *.Dockerfiles in {env_path}, please verify')
        sys.exit(1)

    if args.include_tools:
        tool_path = os.path.abspath(args.conda_envs).replace('conda', 'tools')
        tools = sorted(glob.glob(f'{tool_path}/*/'))
        for i, tool in enumerate(tools):
            tool = os.path.basename(os.path.dirname(tool))
            if not tool.startswith('.'):
                img_name = f"{install_path}/{registry}-bactopia-tools-{tool}-{VERSION}.img" if registry else f"{install_path}/bactopia-tools-{tool}-{VERSION}.img"
                pull_name = f"{docker_prefix}/tools-{tool}:{VERSION}"
                build = True
                if args.envname:
                    if not args.envname == tool:
                        build = False

                if build:
                    if check_needs_build(img_name, force=args.force, is_bactopia=args.is_bactopia):
                        logging.info(f'Found {tool} ({i+1} of {len(env_files)}), begin build to {img_name}')

                        build_singularity_image(img_name, pull_name, max_retry=args.max_retry, force=args.force,
                                                is_bactopia=args.is_bactopia)
