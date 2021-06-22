#! /usr/bin/env python3
"""
usage: bactopia build [-h] [--github] [--quay] [--max_retry INT] [--force]
                      [--verbose] [--silent] [--version]
                      STR STR

setup-docker-builds.py - Build Docker containers for use by Bactopia

positional arguments:
  STR              Directory containing Bactopia repository
  STR              JSON file with latest releases

optional arguments:
  -h, --help       show this help message and exit
  --github         Push to GitHub container registry.
  --quay           Push to Quay.io container registry.
  --max_retry INT  Maximum times to attemp creating Conda environment.
                   (Default: 5)
  --force          Force rebuild of Docker containers.
  --verbose        Print debug related text.
  --silent         Only critical errors will be printed.
  --version        show program's version number and exit
"""
import glob
import json
import logging
import os
import sys

from executor import ExternalCommand, ExternalCommandFailed

PROGRAM = "setup-docker-builds.py"
VERSION = "1.6.0"
REPO = "bactopia"
MAX_RETRY = 5
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
    from executor import ExternalCommand
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
            logging.log(STDERR, e)
            return None


def get_previous_version(json_file):
    """Get the previous version of Bactopia."""
    json_data = None
    with open(json_file, 'rt') as json_fh:
        json_data = json.load(json_fh)

    for node in json_data['repository']['releases']['nodes']:
        this_version = node['name'].lstrip('v')
        if this_version != VERSION:
            return this_version


def check_md5sum(current_md5, image):
    """Compare the two md5 files to see if a rebuild is needed."""
    current = None
    with open(current_md5, 'r') as f:
        current = f.readline().rstrip()

    previous = None
    data = json.loads(execute(f'skopeo inspect docker://docker.io/{image}', capture=True))
    if data:
        if 'conda.md5' in data['Labels']:
            previous = data['Labels']['conda.md5']
            logging.info(f'Found {previous} from {image}')

    logging.info(f'Testing {current} == {previous}')
    return previous == current


def docker_push(image):
    """Push Docker image, with multiple attempts incase of failure."""
    import time
    retry = 0
    allow_fail = False
    success = False
    logging.info(f'Push on {image}')
    while not success:
        result = execute(f'docker push {image}')
        if not result:
            if retry > MAX_RETRY:
                allow_fail = True
            retry += 1
            logging.log(STDERR, "Retrying after short sleep.")
            time.sleep(30 * retry)
        else:
            success = True
    return True


def docker_retag(previous, current, github=False, quay=False):
    """Pull previous version's container, apply current versions to tag."""
    execute(f'docker pull {previous}')
    execute(f'docker tag {previous} {current}')
    docker_push(current)

    if github:
        execute(f'docker tag {previous} ghcr.io/{current}')
        docker_push(f'ghcr.io/{current}')
    if quay:
        execute(f'docker tag {previous} quay.io/{current}')
        docker_push(f'quay.io/{current}')


def docker_tag(image, tag):
    """Tag and push Docker container."""
    logging.info(f'Tagging {tag} to {image}')
    execute(f'docker tag {image} {tag}')
    docker_push(f'{tag}')


def docker_build(recipe, image, latest=None, github=False, quay=False):
    """Build and push latest Docker container."""
    logging.info(f'Building on {image}')
    execute(f'docker build --rm -t {image} -f {recipe} .')
    docker_push(f'{image}')

    if latest:
        docker_tag(image, latest)

    if github:
        docker_tag(image, f'ghcr.io/{image}')
        if latest:
            docker_tag(image, f'ghcr.io/{latest}')

    if quay:
        docker_tag(image, f'quay.io/{image}')
        if latest:
            docker_tag(image, f'quay.io/{latest}')


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='bactopia build',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Build Docker containers for use by Bactopia'
        )
    )

    parser.add_argument('bactopia', metavar="STR", type=str,
                        help='Directory containing Bactopia repository')
    parser.add_argument('releases', metavar="STR", type=str,
                        help='JSON file with latest releases')
    parser.add_argument('--github', action='store_true',
                        help='Push to GitHub container registry.')
    parser.add_argument('--quay', action='store_true',
                        help='Push to Quay.io container registry.')      
    parser.add_argument('--force', action='store_true',
                        help='Force rebuild of Docker containers.')
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
    major, minor, patch = VERSION.split('.')
    previous_version = get_previous_version(args.releases)

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))
    bactopia_path = args.bactopia.rstrip("/")

    # Bactopia Dockerfile
    logging.info(f'Working on Bactopia Dockerfile')
    docker_build(f'{bactopia_path}/Dockerfile', f'{REPO}/bactopia:{VERSION}', latest=f'{REPO}/bactopia:latest',
                 github=args.github, quay=args.quay)

    # Bactopia Process Dockerfiles
    process_files = sorted(glob.glob(f'{bactopia_path}/containers/*.Dockerfile'))
    for i, dockerfile in enumerate(process_files):
        logging.info(f'Working on {dockerfile} ({i+1} of {len(process_files)})')
        process_name = os.path.splitext(os.path.basename(dockerfile))[0]
        latest_image = f'{REPO}/{process_name}:{VERSION}'
        previous_image = f'{REPO}/{process_name}:{previous_version}'
        if check_md5sum(f"{bactopia_path}/conda/linux/{process_name}.md5", previous_image) and not args.force:
            # MD5s match, just need to retag
            logging.info(f'Conda environment did not change, adding tag to previous version')
            docker_retag(previous_image, latest_image, github=args.github, quay=args.quay)
        else:
            # Need to rebuild
            logging.info(f'Conda environment changed, will need to rebuild container')
            docker_build(dockerfile, latest_image, github=args.github, quay=args.quay)

    # Bactopia Tools Dockerfiles
    tools = sorted(glob.glob(f'{bactopia_path}/tools/*/'))
    for i, tool in enumerate(tools):
        tool = os.path.basename(os.path.dirname(tool))
        if not tool.startswith('.'):
            tool_path = f"{bactopia_path}/tools/{tool}"
            dockerfile = f'{tool_path}/Dockerfile'
            latest_image = f'{REPO}/tools-{tool}:{VERSION}'
            previous_image = f'{REPO}/tools-{tool}:{previous_version}'
            logging.info(f'Working on {dockerfile} ({i+1} of {len(tools)})')
            if check_md5sum(f"{tool_path}/environment-linux.md5", previous_image) and not args.force:
                # MD5s match, just need to retag
                logging.info(f'Conda environment did not change, adding tag to previous version')
                docker_retag(previous_image, latest_image, github=args.github, quay=args.quay)
            else:
                # Need to rebuild
                logging.info(f'Conda environment changed, will need to rebuild container')
                docker_build(dockerfile, latest_image, github=args.github, quay=args.quay)
