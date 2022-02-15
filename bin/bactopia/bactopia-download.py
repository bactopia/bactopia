#! /usr/bin/env python3
"""
usage: bactopia download [-h] [--envtype STR] [--wf STR] [--use_defaults] [--build_all] [--registry STR]
                         [--singularity_cache STR] [--singularity_pull_docker_container] [--condadir STR]
                         [--force_rebuild] [--max_retry INT] [--verbose] [--silent] [--version]
                         STR

bactopia download - Build environments for use by Bactopia

positional arguments:
  STR                   Directory containing the Bactopia repo.

optional arguments:
  -h, --help            show this help message and exit
  --condadir STR        Directory to install Conda environments to
  --max_retry INT       Maximum times to attempt creating Conda environment. (Default: 3)

Useful Options:
  --envtype STR         The type of environment to build. (Default: conda)
  --wf STR              Build a environment for a the given workflow
  --use_defaults        Builds environments to the default Bactopia location.
  --build_all           Builds all environments for Bactopia workflows

Container Related Options:
  --registry STR        Docker registry to pull containers from. (Default: quay)
  --singularity_cache STR
                        Location to download Singularity images (Default: ~/.bactopia/singularity)
  --singularity_pull_docker_container
                        Force conversion of Docker containers, instead downloading Singularity images directly

Custom Options:
  --force_rebuild       Force overwrite of existing pre-built environments.
  --verbose             Print debug related text.
  --silent              Only critical errors will be printed.
  --version             show program's version number and exit
"""
import logging
import os
import sys
import time

VERSION = "2.0.1"
PROGRAM = "bactopia download"
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")

BACTOPIA_MODULES = [
    'annotate_genome', 'antimicrobial_resistance', 'ariba_analysis', 'assemble_genome', 'assembly_qc',
    'blast', 'call_variants', 'gather_samples', 'mapping_query', 'minmer_query', 'minmer_sketch', 'qc_reads',
    'sequence_type'
]


def get_platform():
    from sys import platform
    if platform == "darwin":
        return 'mac'
    elif platform == "win32":
        # Windows is not supported
        logging.error("Windows is not supported.", file=sys.stderr)
        sys.exit(1)
    return 'linux'


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
            return [command.decoded_stdout, command.decoded_stderr]
        return True
    except ExternalCommandFailed as e:
        if allow_fail:
            print(e, file=sys.stderr)
            sys.exit(e.returncode)
        else:
            return None


def cleanup_value(value):
    """Remove some characters Nextflow param values"""
    if value.startswith("["):
        # return a list
        return value.lstrip("[").rstrip("]").replace("'", "").replace(",","").split()
    elif value == "true" or value == "false":
        return bool(value)
    else:
        return value.lstrip("'").rstrip("'").replace("\\", "")


def parse_module(main_nf):
    """Pull out the Conda, Docker and singularity info"""
    envs = {}
    with open(main_nf, 'rt') as main_fh:
        read_container = 0
        for line in main_fh:
            line = line.strip()
            if line.startswith("conda_tools"):
                # Conda
                envs['conda'] = line.replace("conda_tools = ", "").replace('"', '')
            elif line.startswith('container') and "containerEngine" in line:
                # Next two lines are singularity and docker
                read_container = 1
            elif line.startswith('container'):
                # There is not singularity image
                envs['singularity'] = False
                envs['docker'] = line.replace("container ", "").replace("'", "")
            elif read_container == 1:
                # Galaxy Project image
                envs['singularity'] = line.replace("'", "").split()[0].strip()
                read_container = 2
            elif read_container == 2:
                # Galaxy Project image
                envs['docker'] = line.replace("'", "").replace('}"', '').strip()
                read_container = 0
    return envs


def parse_workflows(bactopia_path, include_merlin=False, build_all=False):
    """Parse Bactopia's workflows.conf to get modules per-workflow"""
    workflows = {}
    available_workflows = []
    nf_config, stderr = execute(f'nextflow config -flat {bactopia_path}/main.nf', capture=True)
    for line in nf_config.split("\n"):
        if line.startswith('params.available_workflows.') or line.startswith('params.workflows.'):
            param, val = line.split(' = ')
            if line.startswith('params.available_workflows.'):
                # Available workflow definitions
                for wf in cleanup_value(val):
                    available_workflows.append(wf)
            elif line.startswith('params.workflows.'):
                # Workflow definitions
                wf, key = param.replace('params.workflows.', '').split('.')
                if wf not in workflows:
                    workflows[wf] = {}
                workflows[wf][key] = cleanup_value(val)

    # Merged the two
    final_workflows = {}
    for wf in available_workflows:
        final_workflows[wf] = {}
        modules = {}
        if "includes" in workflows[wf]:
            for include in workflows[wf]["includes"]:
                if "modules" in workflows[include]:
                    for module in workflows[include]["modules"]:
                        modules[module] = True
        if "modules" in workflows[wf]:
            for module in workflows[wf]["modules"]:
                modules[module] = True
        if "path" in workflows[wf]:
            modules[wf] = True

        if include_merlin and wf == "bactopia":
            for module in workflows["merlin"]["modules"]:
                modules[module] = True

        for module in modules:
            if module in BACTOPIA_MODULES:
                final_workflows[wf]['bactopia'] = True
            else:
                final_workflows[wf][module] = parse_module(f'{bactopia_path}/{workflows[module]["path"]}/main.nf')

        final_workflows[wf]["custom_dumpsoftwareversions"] = parse_module(f'{bactopia_path}/{workflows["custom_dumpsoftwareversions"]["path"]}/main.nf')
        final_workflows[wf]["csvtk_concat"] = parse_module(f'{bactopia_path}/{workflows["csvtk_concat"]["path"]}/main.nf')

    return final_workflows


def get_docker_prefix(registry):
    """Return the proper prefix based on registry."""
    if registry == "quay":
        return 'quay.io'
    elif registry == "github":
        return 'ghcr.io'
    else:
        return ''


def build_bactopia_envs(bactopia_path, conda_path, singularity_path, env_type, registry_name="quay", force=False, max_retry=5):
    import glob

    # Determine which environment types to build
    build_conda = True if env_type == "conda" or env_type == "all" else False
    build_docker = True if env_type == "docker" or env_type == "all" else False
    build_singularity = True if env_type == "singularity" or env_type == "all" else False

    # setup config variables
    ostype = get_platform()
    major, minor, patch = VERSION.split('.')
    CONTAINER_VERSION = f'{major}.{minor}.x'
    registry = get_docker_prefix(registry_name)
    docker_prefix = f'docker://{registry}/bactopia' if registry else f'docker://bactopia'
    env_files = sorted(glob.glob(f'{bactopia_path}/conda/{ostype}/*.yml'))

    # Check for completion files
    conda_complete = f'{conda_path}/envs-built-{CONTAINER_VERSION}.txt'
    docker_complete = f'{conda_path}/{registry}-containers-pulled-{VERSION}.txt'
    singularity_complete = f'{singularity_path}/{registry}-images-built-{VERSION}.txt'

    if build_conda and os.path.exists(conda_complete):
        if force:
            logging.info(f'--force_rebuild used, overwriting existing Conda environments in {conda_path}')
        else:
            logging.info(f'Found Conda environments in {conda_path}, if a complete rebuild is needed please use --force_rebuild')
    if build_docker and os.path.exists(docker_complete):
        logging.info(f'Found Docker containers, if a complete rebuild is needed please manually remove the containers')
    if build_singularity and os.path.exists(singularity_complete):
        if force:
            logging.info(f'--force_rebuild used, overwriting existing Singularity images in {singularity_path}')
        else:
            logging.info(f'Found Singularity images in {singularity_path}, if a complete rebuild is needed please use --force_rebuild')

    if env_files:
        for i, yml_file in enumerate(env_files):
            # Conda
            envname = os.path.splitext(os.path.basename(yml_file))[0]
            md5_file = yml_file.replace(".yml", ".md5")
            conda_prefix = f'{conda_path}/{envname}-{CONTAINER_VERSION}'
            envbuilt_file = f'{conda_path}/{envname}-{CONTAINER_VERSION}/env-built.txt'

            # Docker/Singularity
            pull_name = f"{docker_prefix}/{envname}:{VERSION}"
            img_name = f"{singularity_path}/{registry}-bactopia-{envname}-{VERSION}.img" if registry else f"{singularity_path}/bactopia-{envname}-{VERSION}.img"
            
            if build_conda:
                if needs_conda_create(envbuilt_file, md5_file, conda_prefix, force=force):
                    logging.info(f'Found {yml_file} ({i+1} of {len(env_files)}), begin build to {conda_prefix}')
                    built = build_conda_env(yml_file, conda_prefix, max_retry=max_retry, force=force)
                    if built:
                        execute(f'cp {md5_file} {envbuilt_file}')
            if build_docker:
                if needs_docker_pull(pull_name):
                    logging.info(f'Found {pull_name} ({i+1} of {len(env_files)}), begin docker pull')
                    docker_pull(pull_name, max_retry=max_retry)
            if build_singularity:
                if needs_singularity_build(img_name, force=force):
                    execute(f'mkdir -p {singularity_path}')
                    logging.info(f'Found {envname} ({i+1} of {len(env_files)}), begin build to {img_name}')
                    build_singularity_image(img_name, pull_name, max_retry=max_retry, force=force, use_build=True)

        # Create completion files
        if build_conda:
            execute(f'date > {conda_complete}')
        if build_docker:
            execute(f'date > {docker_complete}')
        if build_singularity:
            execute(f'date > {singularity_complete}')
    else:
        logging.error(f'Unable to find Bactopia environment files in {env_path}, please verify')
        sys.exit(1)


def build_nfcore_env(envname, envinfo, conda_path, singularity_path, env_type, force=False, max_retry=5, use_build=False):
    # Determine which environment types to build
    build_conda = True if env_type == "conda" or env_type == "all" else False
    build_docker = True if env_type == "docker" or env_type == "all" else False
    build_singularity = True if env_type == "singularity" or env_type == "all" else False

    # Conda
    conda_envname = envinfo['conda'].replace("=", "-").replace(":", "-").replace(" ", "-")
    conda_prefix = f'{conda_path}/{conda_envname}'
    singularity_name = None
    if use_build:
        singularity_name = envinfo['docker'].replace(":","-").replace("/", "-")
    elif not envinfo['singularity']:
        singularity_name = envinfo['docker'].replace(":","-").replace("/", "-")
        use_build = True
    else:
        singularity_name = envinfo['singularity'].replace("https://", "").replace(":","-").replace("/", "-")
    singularity_img = f"{singularity_path}/{singularity_name}.img"

    # Check for completion files
    conda_complete = f'{conda_path}/{conda_envname}/env-built.txt'

    if build_conda and os.path.exists(conda_complete):
        if force:
            logging.debug(f'Overwriting existing Conda environment in {conda_prefix}')
        else:
            logging.debug(f'Found Conda environment in {conda_path}, if a complete rebuild is needed please use --force_rebuild')
            build_conda = False
    if build_docker and not needs_docker_pull(envinfo['docker']):
        if not force:
            logging.debug(f"Found Docker container for {envinfo['docker']}, if a complete rebuild is needed please manually remove the containers")
            build_docker = False
    if build_singularity and os.path.exists(singularity_img):
        if force:
            logging.debug(f'Overwriting existing Singularity image {singularity_img}')
        else:
            logging.debug(f'Found Singularity image {singularity_img}, if a complete rebuild is needed please use --force_rebuild')
            build_singularity = False

    if build_conda:
        logging.info(f'Begin {envname} create to {conda_prefix}')
        build_conda_env(envinfo['conda'], conda_prefix, max_retry=max_retry, force=force)
    if build_docker:
        if needs_docker_pull(envinfo['docker']):
            logging.info(f"Begin docker pull of {envinfo['docker']}")
            docker_pull(envinfo['docker'], max_retry=max_retry)
    if build_singularity:
        if needs_singularity_build(singularity_img, force=force):
            execute(f'mkdir -p {singularity_path}')
            if use_build:
                logging.info(f'Begin {envname} build to {singularity_img}')
                build_singularity_image(singularity_img, f"docker://{envinfo['docker']}", max_retry=max_retry, force=force, use_build=use_build)
            else:
                logging.info(f'Begin {envname} download to {singularity_img}')
                build_singularity_image(singularity_img, envinfo['singularity'], max_retry=max_retry, force=force, use_build=use_build)

    # Create completion files
    if build_conda:
        execute(f'date > {conda_complete}')

"""
Build checks related
"""
def check_md5sum(expected_md5, current_md5):
    """Compare the two md5 files to see if a rebuild is needed."""
    expected = None
    current = None
    with open(expected_md5, 'r') as f:
        expected = f.readline().rstrip()

    with open(current_md5, 'r') as f:
        current = f.readline().rstrip()

    return expected == current


def needs_conda_create(observed_md5, expected_md5, prefix, force=False):
    """Check if a new Conda environment needs to be built."""
    needs_build = False
    if os.path.exists(observed_md5) and not force:
        if check_md5sum(expected_md5, observed_md5):
            logging.debug(f'Existing env ({prefix}) found, skipping unless --force is used')
        else:
            logging.debug(f'Existing env ({prefix}) is out of sync, it will be updated')
            needs_build = True
    else:
        needs_build = True
    return needs_build


def needs_docker_pull(pull_name):
    """Check if a new container needs to be pulled."""
    output = execute(f'docker inspect {pull_name} || true', capture=True)
    if output[1].startswith("Error: No such object"):
        return True

    logging.debug(f'Existing container ({pull_name}) found, skipping unless manually removed')
    return False


def needs_singularity_build(image, force=False):
    """Check if a new image needs to be built."""
    if os.path.exists(image) and not force:
        logging.debug(f'Existing image ({image}) found, skipping unless --force is used')
        return False
    return True


"""
Envrionment creation related
"""
def build_conda_env(conda_env, conda_path, max_retry=5, force=False):
    """Build Conda env, with chance to retry."""
    force = '--force' if force else ''
    retry = 0
    allow_fail = False
    success = False
    while not success:
        result = None
        if conda_env.endswith(".yml"):
            result = execute(f'mamba env create -f {conda_env} --prefix {conda_path} {force}', allow_fail=allow_fail)
        else:
            result = execute(f'mamba create -p {conda_path} -c conda-forge -c bioconda {force} {conda_env}', allow_fail=allow_fail)
        if not result:
            if retry > max_retry:
                allow_fail = True
            retry += 1
            logging.log(STDERR, "Error creating Conda environment, retrying after short sleep.")
            time.sleep(30 * retry)
        else:
            success = True
    return success


def docker_pull(container, max_retry=5):
    """Pull docker container, with chance to retry."""
    retry = 0
    allow_fail = False
    success = False
    while not success:
        result = execute(f'docker pull {container}', allow_fail=allow_fail)
        if not result:
            if retry > max_retry:
                allow_fail = True
            retry += 1
            logging.log(STDERR, "Error pulling container, retrying after short sleep.")
            time.sleep(30 * retry)
        else:
            success = True
    return success


def build_singularity_image(image, pull, max_retry=5, force=False, use_build=False):
    """Build Conda env, with chance to retry."""
    force = '--force' if force else ''
    retry = 0
    allow_fail = False
    success = False
    while not success:
        result = None
        if use_build:
            result = execute(f'singularity build {force} {image} {pull}', allow_fail=allow_fail)
        else:
            # Download from Galaxy Project
            result = execute(f'wget --quiet -O {image} {pull}', allow_fail=allow_fail)
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
    parser = ap.ArgumentParser(
        prog='bactopia download',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Build environments for use by Bactopia'
        )
    )

    # These should be set by Bactopia
    parser.add_argument('bactopia', metavar="STR", type=str,
                        help='Directory containing the Bactopia repo.')

    group1 = parser.add_argument_group('Useful Options')
    group1.add_argument('--envtype',  metavar='STR', type=str, default="conda", choices=['conda', 'docker', 'singularity', 'all'],
                        help='The type of environment to build. (Default: conda)')
    group1.add_argument('--wf', metavar='STR', type=str, default="bactopia",
                        help='Build a environment for a the given workflow')
    group1.add_argument('--use_defaults', action='store_true',
                        help='Builds environments to the default Bactopia location.')
    group1.add_argument('--build_all', action='store_true',
                        help='Builds all environments for Bactopia workflows')

    group2 = parser.add_argument_group('Container Related Options')
    group2.add_argument('--registry', metavar='STR', type=str, default="quay", choices=['dockerhub', 'quay', 'github'],
                        help='Docker registry to pull containers from. (Default: quay)')
    group2.add_argument('--singularity_cache', metavar='STR', type=str, default="~/.bactopia/singularity",
                        help='Location to download Singularity images (Default: ~/.bactopia/singularity)')
    group2.add_argument('--singularity_pull_docker_container', action='store_true', 
                         help='Force conversion of Docker containers, instead downloading Singularity images directly')

    group3 = parser.add_argument_group('Conda Related Options')
    parser.add_argument('--condadir', metavar="STR", type=str, help='Directory to install Conda environments to')

    group8 = parser.add_argument_group('Custom Options')
    group8.add_argument('--force_rebuild', action='store_true', help='Force overwrite of existing pre-built environments.')
    parser.add_argument('--max_retry', metavar='INT', type=int, default=3, 
                        help='Maximum times to attempt creating Conda environment. (Default: 3)')
    group8.add_argument('--verbose', action='store_true', help='Print debug related text.')
    group8.add_argument('--silent', action='store_true', help='Only critical errors will be printed.')
    group8.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args, unknown = parser.parse_known_args()
    include_merlin = True if "--ask_merlin" in unknown else False

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))

    # Install paths
    bactopia_path = os.path.abspath(os.path.expanduser(args.bactopia))
    conda_path = f'{bactopia_path}/conda/envs'
    if args.condadir:
        conda_path = os.path.abspath(os.path.expanduser(args.condadir))
    singularity_path = os.path.abspath(os.path.expanduser(args.singularity_cache))

    # Current Bactopia workflows
    workflow_modules = parse_workflows(bactopia_path, include_merlin=include_merlin, build_all=args.build_all)
    if args.wf not in workflow_modules:
        # Let nextflow handle unknown workflows
        logging.debug(f"{args.wf} is not a known workflow, skipping")
        sys.exit()

    if args.verbose:
        logging.info("Checking if environment pre-builds are needed")
    else:
        logging.info("Checking if environment pre-builds are needed, use --verbose to see full details.")

    for workflow, modules in workflow_modules.items():
        if workflow == args.wf or args.build_all:
            for module, info in modules.items():
                logging.debug(f"Working on {workflow}")
                if module == "bactopia":
                    # Build all (7) bactopia envs
                    build_bactopia_envs(bactopia_path, conda_path, singularity_path, args.envtype, registry_name=args.registry,
                                        force=args.force_rebuild, max_retry=args.max_retry)
                else:
                    # Build nf-core env, one at a time
                    build_nfcore_env(module, info, conda_path, singularity_path, args.envtype, force=args.force_rebuild, 
                                     max_retry=args.max_retry, use_build=args.singularity_pull_docker_container)
