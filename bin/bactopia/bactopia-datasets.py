#! /usr/bin/env python3
import glob
import json
import logging
import os
import sys

from Bio import SeqIO
from executor import ExternalCommand, ExternalCommandFailed

PROGRAM = "bactopia datasets"
VERSION = "2.2.0"
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
    expected = ['pubmlst']
    if sum([k in data for k in expected]) != len(expected):
        logging.debug((f'Existing dataset cache ({CACHE_JSON}) is missing '
                       'expected fields, refreshing.'))
        data = {
            'pubmlst': pubmlst_schemas(pubmlst_file)
        }

        with open(CACHE_JSON, 'w') as cache_fh:
            logging.debug(f'Created dataset cache ({CACHE_JSON})')
            json.dump(data, cache_fh, indent=4, sort_keys=True)

    return data['pubmlst']


def validate_requirements():
    """Validate the required programs are available, if not exit (1)."""
    from shutil import which
    programs = {
        'makeblastdb': which('makeblastdb'),
        'cd-hit': which('cd-hit'),
        'wget': which('wget'),
        'unzip': which('unzip'),
        'gzip': which('gzip')
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


def validate_species(species, skip_spell_check=False):
    """Query input species against ENA to determine if it exists."""
    import requests
    ENDPOINT = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/any-name'
    checks = []

    if os.path.exists(species):
        with open(species, 'r') as handle:
            for line in handle:
                line = line.rstrip()
                if line:
                    checks.append(line)
    elif "," in species:
        checks = species.split(',')
    else:
        checks.append(species)
    
    species_key = {}
    for species in checks:
        if not skip_spell_check:
            species = species.strip()
            r = requests.get(f'{ENDPOINT}/{species}?limit=1')
            if r.status_code == requests.codes.ok:
                try:
                    json_data = r.json()
                    scientific_name = json_data[0]['scientificName'].replace('[', '').replace(']', '')
                    if scientific_name.lower() != species.lower():
                        # Error! Species/Organism found, but doesn't match input. This shouldn't
                        # (query is case-insensitive exact match) happen, but my grandma could "
                        # probably trigger it, so here it is!
                        logging.error((f'Input species ({species}) does not match return result '
                                    f'({scientific_name}), please check spelling.'))
                        sys.exit(1)
                    
                    species_key[species.lower()] = scientific_name
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
        else:
            logging.info(f'{species} skipped ENA Taxonomy spell check')
            species_key[species.lower()] = species
    return species_key


def validate_user_files(user_file, num_species, param):
    """Validate user input files."""
    if not os.path.exists(user_file):
        logging.error(f'Unable to locate {user_file}, please verify path')
        sys.exit(1)
    elif not num_species:
        logging.error(f'A single species (--species) must be given to use --{param}')
        sys.exit(1)
    elif num_species > 1:
        logging.error(f'Multiple species (n={num_species}) can not be used with --{param}')
        sys.exit(1)
    return True


def copy_user_files(source, destination):
    """Copy user files to a new directory."""
    if os.path.exists(source):
        if os.path.isdir(source):
            # Copy all contents in a directory
            logging.info(f'Copying {source}/* to {destination}/')
            execute(f'cp {source}/* {destination}/')
        else:
            # Copy a single file
            logging.info(f'Copying {source} to {destination}/')
            execute(f'cp {source} {destination}/')
    else:
        logging.error(f'{source} does not exist, skipping copy to {destination}.')


def pubmlst_schemas(pubmlst_file):
    """Read the PubMLST mappings and return a dict."""
    pubmlst = {}
    with open(pubmlst_file, 'rt') as pubmlst_fh:
        for line in pubmlst_fh:
            line = line.rstrip()
            if line and not line.startswith('pubmlst'):
                pubmlst_species, species, schema = line.split('\t')
                if species not in pubmlst:
                    pubmlst[species] = {}
                pubmlst[species][schema] = pubmlst_species
    return pubmlst


def available_datasets(pubmlst, missing=False):
    """Print available Ariba references, MLST schemas, and exit."""
    print_to = sys.stderr if missing else sys.stdout

    print("\nMLST schemas available from pubMLST.org:", file=print_to)
    for k,v in sorted(pubmlst.items()):
        if len(v) > 1:
            print(f'{k} ({len(v)} shemas)', file=print_to)
        else:
            print(f'{k}', file=print_to)
    sys.exit(1 if missing else 0)


def available_species(dataset_dir):
    """Print species already in the dataset, and exit."""
    summary_json = f'{dataset_dir}/summary.json'
    if os.path.exists(summary_json):
        species_list = []
        with open(summary_json, 'rt') as summary_fh:
            summary_data = json.load(summary_fh)
            for species in summary_data['species-specific'].keys():
                species_list.append(species.capitalize().replace('-', ' '))
        if species_list:
            print(','.join(species_list))
        else:
            print(f'No available species datasets, exiting...', file=sys.stderr)
        sys.exit(0)
    else:
        print(f'Please verify {dataset_dir} contains Bactopia Datasets, exiting...', file=sys.stderr)
        sys.exit(1)

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


def setup_mlst_request(request, available_schemas, species_key=None):
    """Return a list of mlst schemas to build."""
    requests = []
    if os.path.exists(request):
        with open(request, 'r') as handle:
            for line in handle:
                line = line.rstrip()
                if line:
                    requests.append(line)
    elif "," in request:
        for dataset in request.split(','):
            requests.append(dataset.capitalize().strip())
    else:
        requests.append(request.capitalize())

    schemas = []
    for species in requests:
        species = species_key[species.lower()]
        genus = species.split()[0]
        if species in available_schemas:
            for schema, pubmlst_name in available_schemas[species].items():
                schemas.append({'pubmlst': pubmlst_name, 'schema': schema, 'species': species})
        elif genus in available_schemas:
            # MLST schema is for a genus not just species
            for schema, pubmlst_name in available_schemas[genus].items():
                schemas.append({'pubmlst': pubmlst_name, 'schema': schema, 'species': species})
        else:
            logging.warning(f'{species} is not available from pubMLST.org, skipping')

    return schemas

def setup_mlst(request, available_datasets, pubmlst_json, outdir, force=False, species_key=None):
    """Setup MLST datasets for each requested schema."""
    import re
    requests = setup_mlst_request(request, available_datasets, species_key=species_key)
    if requests:
        for request in requests:
            schema = request['schema']
            species = request['species']

            species = re.sub(r'[ /()]', "-", species.lower())
            species = species.replace('--', '-').strip('-')
            mlst_dir = f'{outdir}/{species}/mlst'
            schema_dir = f'{mlst_dir}/{schema}'
            if os.path.exists(f'{mlst_dir}/{schema}-updated.txt'):
                if force:
                    logging.info(f'--force, removing existing {request["species"]} setup')
                    execute(f'rm -rf {mlst_dir}/{schema}-updated.txt')
                    execute(f'rm -rf {mlst_dir}/{schema}.tar.gz')
                else:
                    logging.info((f'{request["species"]} MLST Schema ({schema}) exists'
                                  ', skipping'))
                    continue
            elif force:
                logging.info(f'--force, removing existing {request["species"]} setup')
                execute(f'rm -rf {mlst_dir}')

            # Setup MLST dataset
            logging.info(f'Setting up {schema} MLST schema for {request["species"]}')
            species_request = request['pubmlst']
            blast_dir = f'{schema_dir}/blastdb'
            execute(f'mkdir -p {blast_dir}')

            # Download Profile
            profile = pubmlst_json[species_request]['profiles']
            execute(f'wget -O profile.txt {profile}', directory=blast_dir)

            # Download loci fastas and make blast db
            loci = profile = pubmlst_json[species_request]['loci']
            for locus, locus_url in loci.items():
                execute(f'wget -O {locus}.tfa {locus_url}', directory=blast_dir)
                execute(f'makeblastdb -in {locus}.tfa -dbtype nucl -out {locus}', directory=blast_dir)

            # Tarball directories
            execute(f'tar -zcvf {schema}.tar.gz {schema}/', directory=mlst_dir)
            execute(f'rm -rf {schema_dir}')

            # Finish up
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > {schema}-updated.txt',
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
                 cpus=1, species_key=None, assembly_level='complete'):
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
    requests = None
    if os.path.exists(request):
        requests = setup_requests(request, available_datasets, 'Prokka Proteins', skip_check=True)
    else:
        requests = setup_requests(request.capitalize(), available_datasets, 'Prokka Proteins', skip_check=True)
    if requests:
        for request in requests:
            species = re.sub(r'[ /()]', "-", request.lower())
            species = species.replace('--', '-').strip('-')
            prokka_dir = f'{outdir}/{species}/annotation'
            minmer_dir = f'{outdir}/{species}/minmer'
            clean_up = False
            genome_sizes = []
            skip_genome_size = False

            if os.path.exists(f'{prokka_dir}/{species}.faa'):
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
            logging.info(f'Downloading genomes (assembly level: {assembly_level})')
            genome_dir = f'{prokka_dir}/genomes'
            genus = species_key[request.lower()]
            execute(f'mkdir {genome_dir}')
            species_accession = []
            all_accessions = {}
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

                results = execute((f'ncbi-genome-download bacteria -g "{genus}" '
                                   f'-l {assembly_level} -F genbank -r 80 --dry-run'), capture=True, error_ok=True)
                
                if results:
                    for line in results.split('\n'):
                        if line and not line.startswith('Considering'):
                            accession, name = line.split('\t', 1)
                            all_accessions[accession] = name
                            if name.startswith(species_key[request.lower()]):
                                species_accession.append(accession)
                            accessions.append(accession)

                    if limit:
                        if len(accessions) > limit:
                            logging.info(f'Downloading {limit} genomes from a random subset of {len(accessions)} genomes.')
                            accessions = random.sample(accessions, limit)
                            contains_species = False
                            for accession in accessions:
                                if all_accessions[accession].startswith(species_key[request.lower()]):
                                    contains_species = True

                            if not contains_species:
                                if len(species_accession):
                                    logging.info(f'Random subset, does not include {species_key[request.lower()]} genomes, adding 1 to random subset.')
                                    accessions.append(random.sample(species_accession, 1)[0])
                        else:
                            logging.info(f'There are less available genomes than the given limit ({limit}), downloading all.')

                    if not len(species_accession):
                        logging.info(f'A completed genome does not exist for {species_key[request.lower()]}, skipping genome size statistics..')
                        skip_genome_size = True
                    
                    with open(accession_file, 'w') as accession_fh:
                        for accession in accessions:
                            accession_fh.write(f'{accession}\n')
                else:
                    logging.warning(f'No completed genomes found for "{genus}", skipping custom Prokka proteins')
                    continue

            execute((f'ncbi-genome-download bacteria -A {accession_file} '
                    f'-l {assembly_level} -o {prokka_dir}/genomes -F genbank -r 80 '
                    f'-m {prokka_dir}/ncbi-metadata.txt'))

            # Extract information from Genbank files
            genbank_files = execute(
                'find . -name "*.gbff.gz"', directory=prokka_dir, capture=True
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
                        if not skip_genome_size:
                            if record.annotations["organism"].lower().startswith(request.lower()):
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
            if not skip_genome_size:
                median_genome = int(median(genome_sizes))
                logging.info(
                    f'Median genome size: {median_genome} (n={total_genome})'
                )
            cdhit_cds = f'{prokka_dir}/{species}.faa'
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
                    'min': 0, 'median': 0, 'mean':0, 'max': 0, 'total': 0,
                    'description': 'No available completed genomes.'
                }
                if not skip_genome_size:
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
                f'sed -i "s={species}.faa:=after_cd-hit\t=" cdhit-stats.txt',
                directory=prokka_dir
            )
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > minmer-updated.txt',
                    directory=minmer_dir)

            # Clean up
            if not keep_files:
                execute(f'rm -rf {minmer} {passing_cds} {genome_dir}/')

    else:
        logging.info("No valid species to setup, skipping")


def setup_amr(outdir, force=False):
    """Download the latest antimicrobial resistance datasets."""
    datasets = ['amrfinderdb']
    amr_dir = f'{outdir}/antimicrobial-resistance'
    update_timestamp = False
    execute(f'mkdir -p {amr_dir}')

    for dataset in datasets:
        dataset_file = f'{amr_dir}/{dataset}.tar.gz'
        if os.path.exists(dataset_file):
            if force:
                logging.info(f'--force, removing existing {dataset_file} setup')
                execute(f'rm -f {dataset_file}')
                update_timestamp = True
            else:
                logging.info(f'{dataset_file} exists, skipping')
                continue

        if dataset == 'amrfinder':
            logging.info(f'Setting up latest AMRFinder+ database')
            prefix = 'amrfinderdb'
            execute(f'rm -rf {prefix} {prefix}-temp', directory=amr_dir)
            execute(f'mkdir -p {prefix} {prefix}-temp', directory=amr_dir)
            execute(f'amrfinder_update -d {prefix}-temp', directory=amr_dir)
            latest_db = os.readlink(f'{amr_dir}/{prefix}-temp/latest')
            execute(f'mv {latest_db}/* {prefix}/', directory=amr_dir)
            execute(f'tar -czvf {prefix}.tar.gz {prefix}/', directory=amr_dir)
            execute(f'rm -rf {prefix} {prefix}-temp', directory=amr_dir)
            execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > {prefix}-updated.txt', directory=amr_dir)
            logging.info(f'AMRFinder+ database saved to {amr_dir}/{prefix}.tar.gz')


def setup_minmer(outdir, force=False, skip_ssl_check=False):
    """Download precomputed Refseq (Mash) and Genbank (Sourmash) datasets."""
    datasets = {
        # Last updated: 2019-03-04
        'sourmash-genbank-k21.json.gz': 'https://osf.io/d7rv8/download',
        'sourmash-genbank-k31.json.gz': 'https://osf.io/4f8n3/download',
        'sourmash-genbank-k51.json.gz': 'https://osf.io/nemkw/download',
        'mash-refseq-k21.msh': (
            'https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh'
        )
    }
    opts = "--no-check-certificate" if skip_ssl_check else ""
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

        execute(f'wget {opts} --quiet -O {filename} {url}', directory=minmer_dir)

    # Finish up
    if update_timestamp or not os.path.exists(f'{minmer_dir}/minmer-updated.txt'):
        execute(f'date -u +"%Y-%m-%dT%H:%M:%SZ" > minmer-updated.txt',
                directory=minmer_dir)


def create_summary(outdir, training_set=False, reference=False, mapping=False, genes=False, proteins=False, primers=False):
    """Create a summary of available datasets in JSON format."""
    from collections import OrderedDict
    available_datasets = OrderedDict()

    available_datasets['antimicrobial-resistance'] = []
    available_datasets['minmer'] = {'sketches': [], 'last_update': None}
    available_datasets['plasmid'] = {'sketches': None, 'blastdb': None, 'last_update': None}

    # Antimicrobial Resistance
    if os.path.exists(f'{outdir}/antimicrobial-resistance'):
        for db in sorted(os.listdir(f'{outdir}/antimicrobial-resistance')):
            if db.endswith(".tar.gz"):
                if db != 'EMPTY.tar.gz':
                    name = db.replace(".tar.gz", "")
                    available_datasets['antimicrobial-resistance'].append({
                        'name': db,
                        'last_update': execute(
                            f'head -n 1 {outdir}/antimicrobial-resistance/{name}-updated.txt', capture=True
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

    # Organism Specific
    if os.path.exists(f'{outdir}/species-specific'):
        available_datasets['species-specific'] = OrderedDict()
        for species in sorted(os.listdir(f'{outdir}/species-specific')):
            if species.startswith("."):
                # Skip hidden files like .DS_Store
                continue
            new_species = OrderedDict()
            species_dir = f'{outdir}/species-specific/{species}'

            minmer = f'{species_dir}/minmer'
            new_species['minmer'] = {'mash': None, 'last_updated': None}
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
            if os.path.exists(f'{prokka}/{species}.faa'):
                new_species['annotation'] = {
                    'proteins': f'species-specific/{species}/annotation/{species}.faa',
                    'last_updated': execute(
                        f'head -n 1 {prokka}/proteins-updated.txt',
                        capture=True
                    ).rstrip()
                }

            if training_set or os.path.exists(f'{prokka}/prodigal.tf'):
                if not os.path.exists(prokka):
                    execute(f'mkdir -p {prokka}')
                if training_set:
                    execute(f'cp {training_set} {prokka}/prodigal.tf')
                new_species['annotation']['training_set'] = f'species-specific/{species}/annotation/prodigal.tf'

            new_species['genome_size'] = {'min': None, 'median': None, 'mean': None, 'max': None}
            if os.path.exists(f'{prokka}/genome_size.json'):
                with open(f'{prokka}/genome_size.json', 'r') as gs_fh:
                    json_data = json.load(gs_fh)
                    new_species['genome_size'] = json_data

            mlst = f'{species_dir}/mlst'
            new_species['mlst'] = {} 
            if os.path.exists(f'{mlst}'):
                for mlst_file in sorted(os.listdir(f'{mlst}')):
                    if mlst_file.endswith("tar.gz"):
                        schema = mlst_file.replace(".tar.gz", "")
                        new_species['mlst'][schema] = {
                            'data': f'species-specific/{species}/mlst/{schema}.tar.gz',
                            'last_updated': execute(f'head -n 1 {mlst}/{schema}-updated.txt', capture=True).rstrip()
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
            
            # Copy over the species specific files
            if reference:
                copy_user_files(reference, f'{species_dir}/optional/reference-genomes')

            if mapping:
                copy_user_files(mapping, f'{species_dir}/optional/mapping-sequences')

            if genes:
                copy_user_files(genes, f'{species_dir}/optional/blast/genes')

            if primers:
                copy_user_files(primers, f'{species_dir}/optional/blast/primers')

            if proteins:
                copy_user_files(proteins, f'{species_dir}/optional/blast/proteins')

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
            stderr_file=None, error_ok=False):
    """A simple wrapper around executor."""
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
    except ExternalCommandFailed as e:
        if "No downloads matched your filter" in e.error_message and error_ok:
            return None
        else:
            print(e)
            sys.exit(1)


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
              {PROGRAM}
              {PROGRAM} --species 'Staphylococcus aureus' --include_genus
        ''')
    )

    parser.add_argument(
        'bactopia', metavar="BACTOPIA_DIR", type=str,
        help='Directory where Bactopia repository is stored.'
    )

    parser.add_argument(
        '--outdir', metavar="STR", type=str, default="./datasets",
        help='Directory to write output. (Default ./datasets)'
    )

    group2 = parser.add_argument_group('Bacterial Species')
    group2.add_argument(
        '--species', metavar="STR", type=str,
        help=('Download available MLST schemas and completed genomes for '
              'a given species or a list of species in a text file.')
    )
    group2.add_argument(
        '--skip_mlst', action='store_true',
        help=('Skip setup of MLST schemas for each species')
    )
    group2.add_argument('--skip_spell_check', action='store_true',
                        help="The spelling of the species name will not be validated against ENA")

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
        '--assembly_level', default='complete', type=str,
        choices=['all', 'complete', 'chromosome', 'scaffold', 'contig', 'complete,chromosome', 'chromosome,complete'],
        help=('Assembly levels of genomes to download (Default: complete).')
    )
    group3.add_argument(
        '--limit', metavar="INT", type=int, default=1000,
        help=('If available completed genomes exceeds a given limit, a random '
              'subsample will be taken. (Default 1000)')
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

    group4 = parser.add_argument_group('Minmer Datasets')
    group4.add_argument(
        '--skip_minmer', action='store_true',
        help='Skip download of pre-computed minmer datasets (mash, sourmash)'
    )

    group6 = parser.add_argument_group('Antimicrobial Resistance Datasets')
    group6.add_argument(
        '--skip_amr', action='store_true',
        help='Skip download of antimicrobial resistance databases (e.g. AMRFinder+)'
    )

    group7 = parser.add_argument_group('Optional User Provided Datasets')
    group7.add_argument(
        '--prodigal_tf', metavar="STR", type=str,
        help=("A pre-built Prodigal training file to add to the species "
              "annotation folder. Requires a single species (--species) and "
              "will replace existing training files.")
    )

    group7.add_argument(
        '--reference', metavar="STR", type=str,
        help=("A reference genome (FASTA/GenBank (preferred)) file or directory "
              "to be added to the optional folder for variant calling. Requires "
              "a single species (--species).")
    )
    group7.add_argument(
        '--mapping', metavar="STR", type=str,
        help=("A reference sequence (FASTA) file or directory to be added to the "
              "optional folder for mapping. Requires a single species (--species).")
    )
    group7.add_argument(
        '--genes', metavar="STR", type=str,
        help=("A gene sequence (FASTA) file or directory to be added to the "
              "optional folder for BLAST. Requires a single species (--species).")
    )
    group7.add_argument(
        '--proteins', metavar="STR", type=str,
        help=("A protein sequence (FASTA) file or directory to be added to the "
              "optional folder for BLAST. Requires a single species (--species).")
    )
    group7.add_argument(
        '--primers', metavar="STR", type=str,
        help=("A primer sequence (FASTA) file or directory to be added to the "
              "optional folder for BLAST. Requires a single species (--species).")
    )
    group7.add_argument(
        '--force_optional', action='store_true',
        help='Overwrite any existing files in the optional folders'
    )

    group8 = parser.add_argument_group('Custom Options')
    group8.add_argument(
        '--cpus', metavar="INT", type=int, default=1,
        help=('Number of cpus to use. (Default: 1)')
    )
    group8.add_argument('--clear_cache', action='store_true',
                        help='Remove any existing cache.')

    group8.add_argument('--force', action='store_true',
                        help='Forcibly overwrite existing datasets.')
    group8.add_argument('--force_mlst', action='store_true',
                        help='Forcibly overwrite existing MLST datasets.')
    group8.add_argument('--force_prokka', action='store_true',
                        help='Forcibly overwrite existing Prokka datasets.')
    group8.add_argument('--force_minmer', action='store_true',
                        help='Forcibly overwrite existing minmer datasets.')
    group8.add_argument('--force_amr', action='store_true',
                        help='Forcibly overwrite existing antimicrobial resistance datasets.')
    group8.add_argument(
        '--keep_files', action='store_true',
        help=('Keep all downloaded and intermediate files.')
    )
    group8.add_argument(
        '--available_datasets', action='store_true',
        help='List MLST schemas available for setup.'
    )
    group8.add_argument(
        '--available_species', action='store_true',
        help=('List species available from current datasets.')
    )

    group8.add_argument('--depends', action='store_true',
                        help='Verify dependencies are installed.')

    group9 = parser.add_argument_group('Adjust Verbosity')
    group9.add_argument('--skip_ssl_check', action='store_true',
                        help="wget will run with --no-check-certificate to skip SSL checks")
    group9.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    group9.add_argument('--verbose', action='store_true',
                        help='Print debug related text.')
    group9.add_argument('--silent', action='store_true',
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

    with open(f'{args.bactopia}/data/pubmlst.json', 'rt') as fh:
        PUBMLST_JSON = json.load(fh)
    PUBMLST_TXT = get_available_datasets(f'{args.bactopia}/data/pubmlst.txt', args.clear_cache)
    if args.available_datasets:
        available_datasets(PUBMLST_TXT)

    if args.available_species:
        available_species(args.outdir)

    species_key = None
    num_species = 0
    if args.species:
        species_key = validate_species(args.species, args.skip_spell_check)
        num_species = len(species_key.keys())

    if args.include_genus:
        if not num_species:
            logging.warning(f'Species (--species) not given, ignoring --include_genus')

    if args.prodigal_tf:
        validate_user_files(args.prodigal_tf, num_species, 'prodigal_tf')

    if args.reference:
        validate_user_files(args.reference, num_species, 'reference')

    if args.mapping:
        validate_user_files(args.mapping, num_species, 'mapping')

    if args.genes:
        validate_user_files(args.genes, num_species, 'genes')

    if args.proteins:
        validate_user_files(args.proteins, num_species, 'proteins')

    if args.primers:
        validate_user_files(args.primers, num_species, 'primers')

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

    if not args.skip_minmer:
        logging.info('Setting up pre-computed Genbank/Refseq minmer datasets')
        setup_minmer(args.outdir, force=(args.force or args.force_minmer), skip_ssl_check=args.skip_ssl_check)
    else:
        logging.info('Skipping minmer dataset step')

    if not args.skip_amr:
        logging.info('Setting up antimicrobial resistance datasets')
        setup_amr(args.outdir, force=(args.force or args.force_amr))
    else:
        logging.info('Skipping antimicrobial resistance dataset step')

    # Organism datasets
    if args.species:
        species_dir = f'{args.outdir}/species-specific'

        if not args.skip_mlst:
            logging.info('Setting up MLST datasets')
            setup_mlst(
                args.species,
                PUBMLST_TXT,
                PUBMLST_JSON,
                species_dir,
                force=(args.force or args.force_mlst),
                species_key=species_key
            )

        if not args.skip_prokka:
            logging.info('Setting up custom Prokka proteins')
            setup_prokka(
                args.species, PUBMLST_TXT, species_dir, cpus=args.cpus,
                include_genus=args.include_genus, limit=args.limit,
                user_accessions=args.accessions, identity=args.identity,
                overlap=args.overlap, max_memory=args.max_memory,
                fast_cluster=args.fast_cluster, keep_files=args.keep_files,
                force=(args.force or args.force_prokka), species_key=species_key, 
                assembly_level=args.assembly_level
            )
        else:
            logging.info('Skipping custom Prokka dataset step')
    else:
        logging.info('No requests for an species, skipping')

    create_summary(args.outdir, training_set=args.prodigal_tf, reference=args.reference,
                mapping=args.mapping, genes=args.genes, proteins=args.proteins, primers=args.primers)
