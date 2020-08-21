#! /usr/bin/env python3
"""
usage: bactopia datasets [-h] [--ariba STR] [--species STR]
                              [--skip_prokka] [--include_genus]
                              [--identity FLOAT] [--overlap FLOAT]
                              [--max_memory INT] [--fast_cluster]
                              [--skip_minmer] [--skip_plsdb] [--cpus INT]
                              [--clear_cache] [--force] [--force_ariba]
                              [--force_mlst] [--force_prokka]
                              [--force_minmer] [--force_plsdb]
                              [--keep_files] [--list_datasets] [--depends]
                              [--version] [--verbose] [--silent]
                              OUTPUT_DIRECTORY

bactopia datasets - Setup datasets for Bactopia

positional arguments:
  OUTPUT_DIRECTORY  Directory to write output.

optional arguments:
  -h, --help        show this help message and exit

Ariba Reference Datasets:
  --ariba STR       Setup Ariba datasets for a given reference or a list of
                    references in a text file. (Default: card,vfdb_core)

Bacterial Species:
  --species STR     Download available MLST schemas and completed genomes
                    for a given species or a list of species in a text file.

Custom Prokka Protein FASTA & Minmer Sketch of Completed Genomes:
  --skip_prokka     Skip creation of a Prokka formatted fasta for each species
  --include_genus   Include all genus members in the Prokka proteins FASTA
  --identity FLOAT  CD-HIT (-c) sequence identity threshold. (Default: 0.9)
  --overlap FLOAT   CD-HIT (-s) length difference cutoff. (Default: 0.8)
  --max_memory INT  CD-HIT (-M) memory limit (in MB). (Default: unlimited
  --fast_cluster    Use CD-HIT's (-g 0) fast clustering algorithm, instead of
                    the accurate but slow algorithm.

Minmer Datasets:
  --skip_minmer     Skip download of pre-computed minmer datasets (mash,
                    sourmash)

PLSDB (Plasmid) BLAST/Sketch:
  --skip_plsdb      Skip download of pre-computed PLSDB datbases (blast, mash)

Helpful Options:
  --cpus INT        Number of cpus to use. (Default: 1)
  --clear_cache     Remove any existing cache.
  --force           Forcibly overwrite existing datasets.
  --force_ariba     Forcibly overwrite existing Ariba datasets.
  --force_mlst      Forcibly overwrite existing MLST datasets.
  --force_prokka    Forcibly overwrite existing Prokka datasets.
  --force_minmer    Forcibly overwrite existing minmer datasets.
  --force_plsdb     Forcibly overwrite existing PLSDB datasets.
  --keep_files      Keep all downloaded and intermediate files.
  --list_datasets   List Ariba reference datasets and MLST schemas
                    available for setup.
  --depends         Verify dependencies are installed.

Adjust Verbosity:
  --version         show program's version number and exit
  --verbose         Print debug related text.
  --silent          Only critical errors will be printed.

example usage:
  bactopia datasets outdir
  bactopia datasets outdir --ariba 'card'
  bactopia datasets outdir --species 'Staphylococcus aureus' --include_genus
"""
import glob
import json
import logging
import os
import sys

from Bio import SeqIO
from executor import ExternalCommand

PROGRAM = "bactopia datasets"
VERSION = "1.4.8"
STDOUT = 11
STDERR = 12
CACHE_DIR = f'{os.path.expanduser("~")}/.bactopia'
CACHE_JSON = f'{CACHE_DIR}/datasets.json'
EXPIRATION = 15 # Refresh db info if cache is older than 15 days
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


def check_cache(clear_cache=False):
    """Check if a local cache exists to avoid re-downloads."""
    import time

    logging.debug(f'Checking for existing cache')
    if not os.path.exists(CACHE_DIR):
        logging.debug(f'Creating cache directory ({CACHE_DIR})')
        execute(f'mkdir -p {CACHE_DIR}')

    cache_data = {}
    if os.path.exists(CACHE_JSON):
        logging.debug(f'Found existing dataset cache ({CACHE_JSON})')
        days_old = (time.time() - os.path.getctime(CACHE_JSON)) // (24 * 3600)
        if days_old >= EXPIRATION or clear_cache:
            logging.debug((f'Deleting {CACHE_JSON}, Reason: older than '
                           f'{EXPIRATION} days or "--clear_cache" used'))
            execute(f'rm {CACHE_JSON}')
        else:
            with open(CACHE_JSON, 'r') as cache_fh:
                cache_data = json.load(cache_fh)

    return cache_data


def get_available_datasets(pubmlst_file, clear_cache):
    """Get a list of available datasets to be set up."""
    data = check_cache(clear_cache=clear_cache)
    expected = ['ariba', 'pubmlst']
    if sum([k in data for k in expected]) != len(expected):
        logging.debug((f'Existing dataset cache ({CACHE_JSON}) is missing '
                       'expected fields, refreshing.'))
        data = {
            'ariba': sorted(ariba_datasets()),
            'pubmlst': pubmlst_schemas(pubmlst_file)
        }

        with open(CACHE_JSON, 'w') as cache_fh:
            logging.debug(f'Created dataset cache ({CACHE_JSON})')
            json.dump(data, cache_fh, indent=4, sort_keys=True)

    return [data['ariba'], data['pubmlst']]


def validate_requirements():
    """Validate the required programs are available, if not exit (1)."""
    from shutil import which
    programs = {
        'ariba': which('ariba'), 'makeblastdb': which('makeblastdb'),
        'cd-hit': which('cd-hit'), 'wget': which('wget'),
        'unzip': which('unzip'), 'gzip': which('gzip')
        # 'mentalist': which('mentalist')
    }

    missing = False
    for prog, path in programs.items():
        if path:
            logging.debug(f'{prog}: command found.')
        else:
            logging.error(f'{prog}: command not found.')
            missing = True

    if missing:
        logging.error("Requirement missing, exiting")
        sys.exit(1)


def validate_species(species):
    """Query input species against ENA to determine if it exists."""
    import requests
    ENDPOINT = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name'
    checks = []
    if "," in species:
        checks = species.split(',')
    else:
        checks.append(species)
    
    for species in checks:
        r = requests.get(f'{ENDPOINT}/{species}?limit=1')
        if r.status_code == requests.codes.ok:
            try:
                json_data = r.json()
                if json_data[0]['scientificName'].lower() != species.lower():
                    # Error! Species/Organism found, but doesn't match input. This shouldn't
                    # (query is case-insensitive exact match) happen, but my grandma could "
                    # probably trigger it, so here it is!
                    logging.error((f'Input species ({species}) does not match return result '
                                f'({json_data[0]["scientificName"]}), please check spelling.'))
                    sys.exit(1)
                logging.info(f'{species} verified in ENA Taxonomy database')
            except json.decoder.JSONDecodeError:
                if r.text == "No results.":
                    logging.error(f'Input species ({species}) not found, please check spelling.')
                    sys.exit(1)
        else:
            # Error! Species/Organism not found. Check spelling?
            # TODO: Implement"Did you mean?" function
            logging.error(f'Input species ({species}) not found, please check spelling.')
            sys.exit(1)

    return len(checks)


def ariba_datasets():
    """Print a list of datasets available with 'ariba getref'."""
    getref_usage = ' '.join([
        line.strip() for line in
        execute('ariba getref --help', capture=True).strip().split('\n')
    ])
    datasets = getref_usage.split('of: ')[1].split(' outprefix')[0]
    return datasets.split()


def pubmlst_schemas(pubmlst_file):
    """Read the PubMLST mappings and return a dict."""
    pubmlst = {}
    with open(pubmlst_file, 'rt') as pubmlst_fh:
        for line in pubmlst_fh:
            if line and not line.startswith('ariba'):
                ariba, species, schema = line.rstrip().split('\t')
                if species not in pubmlst:
                    pubmlst[species] = {}
                pubmlst[species][schema] = ariba
    return pubmlst


def list_datasets(ariba, pubmlst, missing=False):
    """Print available Ariba references, MLST schemas, and exit."""
    print_to = sys.stderr if missing else sys.stdout
    print("Ariba reference datasets available:", file=print_to)
    print("\n".join(sorted(ariba)), file=print_to)

    print("\nMLST schemas available from pubMLST.org:", file=print_to)
    for k,v in sorted(pubmlst.items()):
        if len(v) > 1:
            print(f'{k} ({len(v)} shemas)', file=print_to)
        else:
            print(f'{k}', file=print_to)
    sys.exit(1 if missing else 0)


def setup_requests(request, available_datasets, title, skip_check=False):
    """Return a list of setup requests."""
    datasets = []
    if os.path.exists(request):
        with open(request, 'r') as handle:
            for line in handle:
                dataset = line.rstrip()
                if dataset in available_datasets or skip_check:
                    datasets.append(dataset)
                else:
                    logging.error(f'{dataset} is not available from {title}')
    elif "," in request:
        for dataset in request.split(','):
            dataset = dataset.strip()
            if dataset in available_datasets or skip_check:
                datasets.append(dataset)
            else:
                logging.error(f'{dataset} is not available from {title}')
    elif request in available_datasets or skip_check:
        datasets.append(request)
    else:
        logging.error(f'{request} is not available from {title}')

    return datasets


def setup_ariba(request, available_datasets, outdir, force=False,
                keep_files=False):
    """Setup each of the requested datasets using Ariba."""
    requests = setup_requests(request, available_datasets, 'ariba')
    if requests:
        ariba_dir = f'{outdir}/ariba'
        for request in requests:
            prefix = f'{ariba_dir}/{request}'
            if os.path.exists(f'{prefix}-updated.txt'):
                if force:
                    logging.info(f'--force, removing existing {request} setup')
                    execute(f'rm -rf {prefix}*')
                else:
                    logging.info(f'{request} ({prefix}) exists, skipping')
                    continue
            elif force:
                logging.info(f'--force, removing existing {request} setup')
                execute(f'rm -rf {prefix}*')

            # Setup Ariba dataset
            logging.info(f'Setting up {request} Ariba dataset')
            fa = f'{prefix}.fa'
            tsv = f'{prefix}.tsv'
            execute(f'mkdir -p {ariba_dir}')
            with open(f'{prefix}-log.txt', 'w') as ariba_log:
                execute(
                    f'ariba getref {request} {request}',
                    stdout_file=ariba_log, stderr_file=ariba_log,
                    directory=ariba_dir
                )
            execute(f'ariba prepareref -f {fa} -m {tsv} {prefix}')

            # Clean up
            if not keep_files:
                execute(f'rm {fa} {tsv}')
            execute(f'mv {request}*.* {request}/', directory=ariba_dir)
            execute(f'tar -zcvf {request}.tar.gz {request}/',
                    directory=ariba_dir)
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > {request}-updated.txt',
                    directory=ariba_dir)
            execute(f'rm -rf {request}', directory=ariba_dir)
    else:
        logging.info("No valid Ariba datasets to setup, skipping")


def setup_mlst_request(request, available_schemas):
    """Return a list of mlst schemas to build."""
    requests = []
    if os.path.exists(request):
        with open(request, 'r') as handle:
            for line in handle:
                requests.append(line.rstrip())
    elif "," in request:
        for dataset in request.split(','):
            requests.append(dataset.capitalize().strip())
    else:
        requests.append(request.capitalize())

    schemas = []
    for species in requests:
        genus = species.split()[0]
        if species in available_schemas:
            for schema, ariba_name in available_schemas[species].items():
                schemas.append({'ariba': ariba_name, 'schema': schema, 'species': species})
        elif genus in available_schemas:
            # MLST schema is for a genus not just species
            for schema, ariba_name in available_schemas[genus].items():
                schemas.append({'ariba': ariba_name, 'schema': schema, 'species': species})
        else:
            logging.error(f'{species} is not available from pubMLST.org, skipping')

    return schemas

def setup_mlst(request, available_datasets, outdir, force=False):
    """Setup MLST datasets for each requested schema."""
    import re
    requests = setup_mlst_request(request, available_datasets)
    if requests:
        for request in requests:
            schema = request['schema']
            species = request['species']

            species = re.sub(r'[ /()]', "-", species.lower())
            species = species.replace('--', '-').strip('-')
            mlst_dir = f'{outdir}/{species}/mlst/{schema}'
            if os.path.exists(f'{mlst_dir}/mlst-updated.txt'):
                if force:
                    logging.info(f'--force, removing existing {request["species"]} setup')
                    execute(f'rm -rf {mlst_dir}')
                else:
                    logging.info((f'{request["species"]}MLST Schema ({mlst_dir}) exists'
                                  ', skipping'))
                    continue
            elif force:
                logging.info(f'--force, removing existing {request["species"]}setup')
                execute(f'rm -rf {mlst_dir}')

            # Setup MLST dataset
            logging.info(f'Setting up {schema} MLST schema for {request["species"]}')
            execute(f'mkdir -p {mlst_dir}')

            # Ariba
            species_request = request['ariba']
            logging.info(f'Creating Ariba MLST dataset')
            ariba_dir = f'{mlst_dir}/ariba'
            execute(f'ariba pubmlstget "{species_request}" {ariba_dir}')

            # BLAST
            logging.info(f'Creating BLAST MLST dataset')
            blast_dir = f'{mlst_dir}/blastdb'
            for fasta in glob.glob(f'{ariba_dir}/pubmlst_download/*.tfa'):
                output = os.path.splitext(fasta)[0]
                execute(f'makeblastdb -in {fasta} -dbtype nucl -out {output}')
            execute(f'mv {ariba_dir}/pubmlst_download {blast_dir}')

            # Tarball directories
            execute(f'tar -zcvf {schema}-ariba.tar.gz ariba/', directory=mlst_dir)
            execute(f'rm -rf {ariba_dir}')
            execute(f'tar -zcvf {schema}-blastdb.tar.gz blastdb/', directory=mlst_dir)
            execute(f'rm -rf {blast_dir}')

            # Finish up
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > mlst-updated.txt',
                    directory=mlst_dir)
    else:
        logging.info("No valid MLST schemas to setup, skipping")


def process_cds(cds):
    """Look over the CDS attributes and return passing CDS."""
    header = None
    seq = None
    qualifiers = cds.keys()
    ec_number = ''
    gene = ''
    product = ''
    is_pseudo = ('pseudo' in qualifiers or 'pseudogene' in qualifiers)
    is_hypothetical = (product.lower() == "hypothetical protein")
    if not is_pseudo and not is_hypothetical:
        if 'ec_number' in qualifiers:
            ec_number = cds['ec_number'][0]
        if 'gene' in qualifiers:
            gene = cds['gene'][0]
        if 'product' in qualifiers:
            product = cds['product'][0]
        if 'protein_id' in qualifiers:
            protein_id = cds['protein_id'][0]
        elif 'locus_tag' in qualifiers:
            protein_id = cds['locus_tag'][0]

        header = f'>{protein_id} {ec_number}~~~{gene}~~~{product}'
        seq = cds['translation'][0]


    return [header, seq]


def setup_prokka(request, available_datasets, outdir, force=False,
                 include_genus=False, limit=None, user_accessions=None, identity=0.9, 
                 overlap=0.8, max_memory=0, fast_cluster=False, keep_files=False, 
                 cpus=1):
    """
    Setup a Prokka compatible protein fasta file based on completed genomes.

    Implemented similar approach as Thanh Lê's "make_prokka_db". Check out
    his version for a standalone implementation!
    Github Repo: https://github.com/thanhleviet/make_prokka_db
    """
    import gzip
    import re
    import random
    from statistics import median, mean
    requests = setup_requests(request.capitalize(), available_datasets, 'Prokka Proteins',
                              skip_check=True)
    if requests:
        for request in requests:
            species = re.sub(r'[ /()]', "-", request.lower())
            species = species.replace('--', '-').strip('-')
            prokka_dir = f'{outdir}/{species}/annotation'
            minmer_dir = f'{outdir}/{species}/minmer'
            clean_up = False
            genome_sizes = []

            if os.path.exists(f'{prokka_dir}/proteins.faa'):
                if force:
                    logging.info(f'--force, delete existing {prokka_dir}')
                    clean_up = True
                else:
                    logging.info((f'{prokka_dir} exists, skipping'))
                    continue
            elif os.path.exists(f'{prokka_dir}/'):
                logging.info(f'Incomplete setup, deleting {prokka_dir} to start over')
                clean_up = True
            elif force:
                logging.info(f'--force, delete existing {prokka_dir}')
                clean_up = True

            if clean_up:
                execute(f'rm -rf {prokka_dir}')
                execute(f'rm -rf {minmer_dir}')

            # Setup Prokka proteins file
            logging.info(f'Setting up custom Prokka proteins for {request}')
            execute(f'mkdir -p {prokka_dir}')
            execute(f'mkdir -p {minmer_dir}')

            # Download completed genomes
            logging.info(f'Downloading completed genomes')
            genome_dir = f'{prokka_dir}/genomes'
            genus = ' '.join(request.split()[0:2])
            execute(f'mkdir {genome_dir}')
            accessions = []
            accession_file = f'{genome_dir}/accessions.txt'
            if user_accessions:
                execute(f'cp {user_accessions} {accession_file}')
                if include_genus:
                    logging.info(f'Ignoring `--include_genus` since a file of accessions was given.')
                if limit:
                    logging.info(f'Ignoring `--limit {limit}` since a file of accessions was given.')
            else:
                if include_genus:
                    genus = genus.split()[0]

                results = execute((f'ncbi-genome-download bacteria --genera "{genus}" '
                                f'-l complete -F genbank -r 80 --dry-run'), capture=True)

                for line in results.split('\n'):
                    if line and not line.startswith('Considering'):
                        accessions.append(line.split()[0])

                accessions = accessions
                if limit:
                    if len(accessions) > limit:
                        logging.info(f'Downloading {limit} genomes from a random subset of {len(accessions)} genomes.')
                        accessions = random.sample(accessions, limit)
                    else:
                        logging.info(f'There are less available genomes than the given limit ({limit}), downloading all.')
            
                with open(accession_file, 'w') as accession_fh:
                    for accession in accessions:
                        accession_fh.write(f'{accession}\n')

            execute((f'ncbi-genome-download bacteria -A {accession_file} '
                    f'-l complete -o {prokka_dir}/genomes -F genbank -r 80 '
                    f'-m {prokka_dir}/ncbi-metadata.txt'))

            # Extract information from Genbank files
            genbank_files = execute(
                'find -name "*.gbff.gz"', directory=prokka_dir, capture=True
            ).split('\n')
            count = 0
            passing_cds = f'{prokka_dir}/passing-cds.faa'
            minmer = f'{minmer_dir}/minmer.ffn'
            logging.info(f'Processing {len(genbank_files)-1} Genbank files')
            with open(passing_cds, 'w') as cds_fh, open(minmer, 'w') as ffn_fh:
                for genbank in genbank_files:
                    if genbank:
                        sizes = []
                        genbank = genbank.replace('./', f'{prokka_dir}/')
                        seq_name = None
                        seqs = []
                        gap = "N" * 102
                        with gzip.open(genbank, 'rt') as genbank_fh:
                            for record in SeqIO.parse(genbank_fh, 'genbank'):
                                # Aggregate chromosome and plasmids
                                sizes.append(len(record.seq))
                                for dbxref in record.dbxrefs:
                                    if dbxref.startswith('Assembly'):
                                        seq_name = dbxref.split(':')[1]
                                        seqs.append(str(record.seq))
                                        seqs.append(gap)

                                for feature in record.features:
                                    if feature.type == 'CDS':
                                        header, seq = process_cds(
                                            feature.qualifiers
                                        )

                                        if header and seq:
                                            count += 1
                                            cds_fh.write(f'{header}\n')
                                            cds_fh.write(f'{seq}\n')
                            # Write sequence
                            ffn_fh.write(f'>{seq_name}\n')
                            gap = "N" * 102
                            sequence = "".join(seqs)
                            ffn_fh.write(f'{sequence}\n')

                        # Only add genome sizes for the species, incase the
                        # option '--inlude_genus' was used.
                        if record.annotations["organism"].startswith(request):
                            logging.debug(
                                f'Added {record.annotations["organism"]} '
                                f'({sum(sizes)}) to median genome size '
                                'calculation.'
                            )
                            genome_sizes.append(sum(sizes))
                        else:
                            logging.debug(
                                f'Skip adding {record.annotations["organism"]} '
                                f'({sum(sizes)}) to median genome size '
                                f'calculation (not {request}).'
                            )

            total_genome = len(genome_sizes)
            median_genome = int(median(genome_sizes))
            logging.info(
                f'Median genome size: {median_genome} (n={total_genome})'
            )
            cdhit_cds = f'{prokka_dir}/proteins.faa'
            logging.info(f'Running CD-HIT on {count} proteins')
            g = 0 if fast_cluster else 1
            execute((f'cd-hit -i {passing_cds} -o {cdhit_cds} -s {overlap} '
                     f'-g {g} -c {identity} -T {cpus} -M {max_memory}'))

            # Make sketch/signatures
            execute(
                f'mash sketch -i -k 31 -s 10000 -o refseq-genomes minmer.ffn',
                directory=minmer_dir
            )

            # Finish up
            with open(f'{prokka_dir}/genome_size.json', 'w') as genome_size_fh:
                gs_dict = {
                    'min': min(genome_sizes),
                    'median': int(median(genome_sizes)),
                    'mean': int(median(genome_sizes)),
                    'max': max(genome_sizes),
                    'total': total_genome,
                    'description': (
                        f'Genome size values are based on {total_genome} '
                        'completed genomes (RefSeq).'
                    )
                }
                json.dump(gs_dict, genome_size_fh, indent=4)
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > proteins-updated.txt',
                    directory=prokka_dir)
            execute(f'grep -H -c "^>" *.faa > cdhit-stats.txt',
                    directory=prokka_dir)
            execute(f'sed -i "s=passing-cds.faa:=original\t=" cdhit-stats.txt',
                    directory=prokka_dir)
            execute(
                f'sed -i "s=proteins.faa:=after_cd-hit\t=" cdhit-stats.txt',
                directory=prokka_dir
            )
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > minmer-updated.txt',
                    directory=minmer_dir)

            # Clean up
            if not keep_files:
                execute(f'rm -rf {minmer} {passing_cds} {genome_dir}/')

    else:
        logging.info("No valid species to setup, skipping")


def setup_minmer(outdir, force=False):
    """Download precomputed Refseq (Mash) and Genbank (Sourmash) datasets."""
    datasets = {
        # Last updated: 2019-03-04
        'genbank-k21.json.gz': 'https://osf.io/d7rv8/download',
        'genbank-k31.json.gz': 'https://osf.io/4f8n3/download',
        'genbank-k51.json.gz': 'https://osf.io/nemkw/download',
        'refseq-k21-s1000.msh': (
            'https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh'
        )
    }

    minmer_dir = f'{outdir}/minmer'
    update_timestamp = False
    if force:
        logging.info(f'--force, removing existing {minmer_dir} setup')
        execute(f'rm -rf {minmer_dir}')

    execute(f'mkdir -p {minmer_dir}')
    for filename, url in datasets.items():
        filepath = f'{minmer_dir}/{filename}'
        if os.path.exists(filepath):
            if force:
                logging.info(f'--force, removing existing {filepath} setup')
                execute(f'rm -rf {filepath}')
                update_timestamp = True
            else:
                logging.info(f'{filepath} exists, skipping')
                continue

        execute(f'wget --quiet -O {filename} {url}', directory=minmer_dir)

    # Finish up
    if update_timestamp or not os.path.exists(f'{minmer_dir}/minmer-updated.txt'):
        execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > minmer-updated.txt',
                directory=minmer_dir)


def setup_plsdb(outdir, keep_files=False, force=False):
    """Download precomputed PLSDB datasets."""
    url = 'https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/?zip'
    plsdb_dir = f'{outdir}/plasmid'
    if os.path.exists(plsdb_dir):
        if force:
            logging.info(f'--force, removing existing {plsdb_dir} setup')
            execute(f'rm -rf {plsdb_dir}')
        else:
            logging.info(f'{plsdb_dir} exists, skipping')
            return None

    execute(f'mkdir -p {plsdb_dir}')
    execute(f'wget --quiet -O plsdb.zip {url}', directory=plsdb_dir)
    execute('unzip plsdb.zip', directory=plsdb_dir)
    execute('ls > plsdb-orginal-names.txt', directory=plsdb_dir)

    # Rename files to generic prefix
    mash_file = os.path.basename(glob.glob(f'{plsdb_dir}/*.msh')[0])
    prefix = mash_file.replace('.msh', '')
    for plsdb_file in os.listdir(plsdb_dir):
        if plsdb_file.startswith(prefix) and prefix != 'plsdb':
            new_name = plsdb_file.replace(prefix, 'plsdb')
            execute(f'mv {plsdb_file} {new_name}', directory=plsdb_dir)

    # Clean up
    if not keep_files:
        execute('rm plsdb.zip', directory=plsdb_dir)

    # Finish up
    execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > plsdb-updated.txt',
            directory=plsdb_dir)


def create_summary(outdir, training_set=False):
    """Create a summary of available datasets in JSON format."""
    from collections import OrderedDict
    available_datasets = OrderedDict()

    # Ariba
    if os.path.exists(f'{outdir}/ariba'):
        available_datasets['ariba'] = []
        for db in sorted(os.listdir(f'{outdir}/ariba')):
            if db.endswith(".tar.gz"):
                name = db.replace(".tar.gz", "")
                available_datasets['ariba'].append({
                    'name': db,
                    'last_update': execute(
                        f'head -n 1 {outdir}/ariba/{name}-updated.txt', capture=True
                    ).rstrip()
                })

    # Minmers
    if os.path.exists(f'{outdir}/minmer/minmer-updated.txt'):
        available_datasets['minmer'] = {
            'sketches': [],
            'last_update': execute(
                f'head -n 1 {outdir}/minmer/minmer-updated.txt', capture=True
            ).rstrip()
        }
        for sketch in sorted(os.listdir(f'{outdir}/minmer')):
            if sketch != 'minmer-updated.txt':
                available_datasets['minmer']['sketches'].append(sketch)

    # PLSDB (plasmids)
    if os.path.exists(f'{outdir}/plasmid/plsdb-updated.txt'):
        available_datasets['plasmid'] = {
            'sketches': 'plsdb.msh',
            'blastdb': 'plsdb.fna',
            'last_update': execute(
                f'head -n 1 {outdir}/plasmid/plsdb-updated.txt', capture=True
            ).rstrip()
        }

    # Organism Specific
    if os.path.exists(f'{outdir}/species-specific'):
        available_datasets['species-specific'] = OrderedDict()
        for species in sorted(os.listdir(f'{outdir}/species-specific')):
            new_species = OrderedDict()
            new_species['mlst'] = []
            species_dir = f'{outdir}/species-specific/{species}'

            minmer = f'{species_dir}/minmer'
            if os.path.exists(f'{minmer}/refseq-genomes.msh'):
                new_species['minmer'] = {
                    'mash': f'species-specific/{species}/minmer/refseq-genomes.msh',
                    'last_updated': execute(
                        f'head -n 1 {minmer}/minmer-updated.txt',
                        capture=True
                    ).rstrip()
                }

            prokka = f'{species_dir}/annotation'
            new_species['annotation'] = { 'proteins': None, 'training_set': None, 'last_updated': None}
            if os.path.exists(f'{prokka}/proteins.faa'):
                new_species['annotation'] = {
                    'proteins': f'species-specific/{species}/annotation/proteins.faa',
                    'last_updated': execute(
                        f'head -n 1 {prokka}/proteins-updated.txt',
                        capture=True
                    ).rstrip()
                }

            if training_set:
                if not os.path.exists(prokka):
                    execute(f'mkdir -p {prokka}')
                execute(f'cp {training_set} {prokka}/prodigal.tf')
                new_species['annotation']['training_set'] = f'species-specific/{species}/annotation/prodigal.tf'

            if os.path.exists(f'{prokka}/genome_size.json'):
                with open(f'{prokka}/genome_size.json', 'r') as gs_fh:
                    json_data = json.load(gs_fh)
                    new_species['genome_size'] = json_data

            mlst = f'{species_dir}/mlst'
            if os.path.exists(f'{mlst}'):
                new_species['mlst'] = {}
                for schema in sorted(os.listdir(f'{mlst}')):
                    if os.path.exists(f'{mlst}/{schema}/{schema}-ariba.tar.gz'):
                        new_species['mlst'][schema] = {
                            'ariba': f'species-specific/{species}/mlst/{schema}/{schema}-ariba.tar.gz',
                            'blast': f'species-specific/{species}/mlst/{schema}/{schema}-blastdb.tar.gz',
                            'last_updated': execute(
                                f'head -n 1 {mlst}/{schema}/mlst-updated.txt', capture=True
                            ).rstrip()
                        }

            optionals = sorted([
                'reference-genomes', 'mapping-sequences', 'blast'
            ])
            new_species['optional'] = OrderedDict()
            for optional in optionals:
                # These are optional directories users can add data to
                optional_dir = f'species-specific/{species}/optional/{optional}'
                if not os.path.exists(optional_dir):
                    execute(f'mkdir -p {optional_dir}', directory=outdir)
                if optional == 'blast':
                    new_species['optional'][optional] = [
                        f'{optional_dir}/genes',
                        f'{optional_dir}/primers',
                        f'{optional_dir}/proteins',
                    ]
                    for blast_dir in new_species['optional'][optional]:
                        execute(f'mkdir -p {blast_dir}', directory=outdir)
                else:
                    new_species['optional'][optional] = f'{optional_dir}'

            available_datasets['species-specific'][species] = new_species

    with open(f'{outdir}/summary.json', 'w') as json_handle:
        logging.info(f'Writing summary of available datasets')
        json.dump(available_datasets, json_handle, indent=4)
        logging.debug(json.dumps(available_datasets, indent=4))


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture=False, stdout_file=None,
            stderr_file=None):
    """A simple wrapper around executor."""
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
    import textwrap
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Setup public datasets for Bactopia'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} outdir
              {PROGRAM} outdir --ariba 'vfdb_core'
              {PROGRAM} outdir --species 'Staphylococcus aureus' --include_genus
        ''')
    )

    parser.add_argument(
        'pubmlst', metavar="PUBMLST", type=str,
        help='Bactopia config file with PubMLST schema mappings for Ariba.'
    )

    parser.add_argument(
        'outdir', metavar="OUTPUT_DIRECTORY", type=str,
        help='Directory to write output.'
    )

    group1 = parser.add_argument_group('Ariba Reference Datasets')
    group1.add_argument(
        '--skip_ariba', action='store_true',
        help=('Skip setup of Ariba datasets')
    )
    group1.add_argument(
        '--ariba', metavar="STR", type=str, default='vfdb_core',
        help=('Setup Ariba datasets for a given reference or a list of '
              'references in a text file. (Default: vfdb_core)')
    )

    group2 = parser.add_argument_group('Bacterial Species')
    group2.add_argument(
        '--species', metavar="STR", type=str,
        help=('Download available MLST schemas and completed genomes for '
              'a given species or a list of species in a text file.')
    )

    group3 = parser.add_argument_group('Custom Prokka Protein FASTA')
    group3.add_argument(
        '--skip_prokka', action='store_true',
        help=('Skip creation of a Prokka formatted fasta for each species')
    )
    group3.add_argument(
        '--include_genus', action='store_true',
        help=('Include all genus members in the Prokka proteins FASTA')
    )
    group3.add_argument(
        '--limit', metavar="INT", type=int, default=0,
        help=('If available completed genomes exceeds a given limit, a random '
              'subsample will be taken.')
    )
    group3.add_argument(
        '--accessions', metavar="STR", type=str,
        help=('A list of RefSeq accessions to download.')
    )
    group3.add_argument(
        '--identity', metavar="FLOAT", type=float, default=0.9,
        help=('CD-HIT (-c) sequence identity threshold. (Default: 0.9)')
    )
    group3.add_argument(
        '--overlap', metavar="FLOAT", type=float, default=0.8,
        help=('CD-HIT (-s) length difference cutoff. (Default: 0.8)')
    )
    group3.add_argument(
        '--max_memory', metavar="INT", type=int, default=0,
        help=('CD-HIT (-M) memory limit (in MB). (Default: unlimited')
    )
    group3.add_argument(
        '--fast_cluster', action='store_true',
        help=("Use CD-HIT's (-g 0) fast clustering algorithm, instead of the "
              "accurate but slow algorithm.")
    )
    group3.add_argument(
        '--prodigal_tf', metavar="STR", type=str,
        help=("A pre-built Prodigal training file to add to the species "
              "annotation folder. Requires a single species (--species) and "
              "will replace existing training files.")
    )

    group4 = parser.add_argument_group('Minmer Datasets')
    group4.add_argument(
        '--skip_minmer', action='store_true',
        help='Skip download of pre-computed minmer datasets (mash, sourmash)'
    )

    group5 = parser.add_argument_group('PLSDB (Plasmid) BLAST/Sketch')
    group5.add_argument(
        '--skip_plsdb', action='store_true',
        help='Skip download of pre-computed PLSDB datbases (blast, mash)'
    )

    group6 = parser.add_argument_group('Helpful Options')
    group6.add_argument(
        '--cpus', metavar="INT", type=int, default=1,
        help=('Number of cpus to use. (Default: 1)')
    )
    group6.add_argument('--clear_cache', action='store_true',
                        help='Remove any existing cache.')

    group6.add_argument('--force', action='store_true',
                        help='Forcibly overwrite existing datasets.')
    group6.add_argument('--force_ariba', action='store_true',
                        help='Forcibly overwrite existing Ariba datasets.')
    group6.add_argument('--force_mlst', action='store_true',
                        help='Forcibly overwrite existing MLST datasets.')
    group6.add_argument('--force_prokka', action='store_true',
                        help='Forcibly overwrite existing Prokka datasets.')
    group6.add_argument('--force_minmer', action='store_true',
                        help='Forcibly overwrite existing minmer datasets.')
    group6.add_argument('--force_plsdb', action='store_true',
                        help='Forcibly overwrite existing PLSDB datasets.')
    group6.add_argument(
        '--keep_files', action='store_true',
        help=('Keep all downloaded and intermediate files.')
    )
    group6.add_argument(
        '--list_datasets', action='store_true',
        help=('List Ariba reference datasets and MLST schemas '
              'available for setup.')
    )

    group6.add_argument('--depends', action='store_true',
                        help='Verify dependencies are installed.')

    group7 = parser.add_argument_group('Adjust Verbosity')
    group7.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    group7.add_argument('--verbose', action='store_true',
                        help='Print debug related text.')
    group7.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.verbose))
    if args.depends:
        validate_requirements()
        sys.exit(0)
    else:
        validate_requirements()

    num_species = 0
    if args.species:
        num_species = validate_species(args.species)

    if args.include_genus:
        if not num_species:
            logging.error(f'Species (--species) not given, ignoring --include_genus')
            sys.exit(1)
        elif num_species > 1:
            logging.error(f'Only a single species (given {num_species}) can be used with --include_genus')
            sys.exit(1)

    if args.prodigal_tf:
        if not os.path.exists(args.prodigal_tf):
            logging.error(f'Unable to locate {args.prodigal_tf}, please verify path')
            sys.exit(1)
        elif not num_species:
            logging.error(f'A single species (--species) must be given to use --prodigal_tf')
            sys.exit(1)
        elif num_species > 1:
            logging.error(f'Only a single species (given {num_species}) can be used with --prodigal_tf')
            sys.exit(1)

    if args.accessions:
        if not os.path.exists(args.accessions):
            logging.error(f'Unable to locate {args.accessions}, please verify path')
            sys.exit(1)
        elif not num_species:
            logging.error(f'A single species (--species) must be given to use --accessions')
            sys.exit(1)
        elif num_species > 1:
            logging.error(f'Only a single species (given {num_species}) can be used with --accessions')
            sys.exit(1)
            
    ARIBA, PUBMLST = get_available_datasets(args.pubmlst, args.clear_cache)
    if args.list_datasets:
        list_datasets(ARIBA, PUBMLST)

    if not args.skip_ariba:
        if args.ariba:
            logging.info('Setting up Ariba datasets')
            setup_ariba(
                args.ariba, ARIBA, args.outdir, keep_files=args.keep_files,
                force=(args.force or args.force_ariba)
            )
        else:
            logging.info('No requests for an Ariba dataset, skipping')
    else:
        logging.info('Skipping Ariba dataset step')

    if not args.skip_minmer:
        logging.info('Setting up pre-computed Genbank/Refseq minmer datasets')
        setup_minmer(args.outdir, force=(args.force or args.force_minmer))
    else:
        logging.info('Skipping minmer dataset step')

    if not args.skip_plsdb:
        logging.info('Setting up pre-computed PLSDB (plasmids) datasets')
        setup_plsdb(args.outdir, keep_files=args.keep_files,
                    force=(args.force or args.force_plsdb))
    else:
        logging.info('Skipping PLSDB (plasmids) dataset step')

    # Organism datasets
    if args.species:
        species_dir = f'{args.outdir}/species-specific'
        logging.info('Setting up MLST datasets')
        setup_mlst(args.species, PUBMLST, species_dir,
                   force=(args.force or args.force_mlst))

        if not args.skip_prokka:
            logging.info('Setting up custom Prokka proteins')
            setup_prokka(
                args.species, PUBMLST, species_dir, cpus=args.cpus,
                include_genus=args.include_genus, limit=args.limit,
                user_accessions=args.accessions, identity=args.identity,
                overlap=args.overlap, max_memory=args.max_memory,
                fast_cluster=args.fast_cluster, keep_files=args.keep_files,
                force=(args.force or args.force_prokka)
            )
        else:
            logging.info('Skipping custom Prokka dataset step')
    else:
        logging.info('No requests for an species, skipping')

    create_summary(args.outdir, training_set=args.prodigal_tf)
