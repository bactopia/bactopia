#! /usr/bin/env python3
"""
Query Taxon ID or Study accession against ENA and return a list of WGS results.

usage: bactopia search [-h] [--exact_taxon] [--outdir OUTPUT_DIRECTORY]
                       [--prefix PREFIX] [--limit INT] [--version]
                       STR

bactopia search - Search ENA for associated WGS samples

positional arguments:
  STR                   Taxon ID or Study accession

optional arguments:
  -h, --help            show this help message and exit
  --exact_taxon         Exclude Taxon ID descendents.
  --outdir OUTPUT_DIRECTORY
                        Directory to write output. (Default: .)
  --prefix PREFIX       Prefix to use for output file names. (Default: ena)
  --limit INT           Maximum number of results to return. (Default:
                        1000000)
  --version             show program's version number and exit

example usage:
  bactopia search PRJNA480016 --limit 20
  bactopia search 1280 --exact_taxon --limit 20'
  bactopia search "staphylococcus aureus" --limit 20

"""
import sys
VERSION = "1.4.8"
PROGRAM = "bactopia search"
ENA_URL = ('https://www.ebi.ac.uk/ena/data/warehouse/search?result=read_run&'
           'display=report')
FIELDS = [
    'study_accession', 'secondary_study_accession', 'sample_accession',
    'secondary_sample_accession', 'experiment_accession', 'run_accession',
    'submission_accession', 'tax_id', 'scientific_name',
    'instrument_platform', 'instrument_model', 'library_name',
    'library_layout', 'nominal_length', 'library_strategy',
    'library_source', 'library_selection', 'read_count',
    'base_count', 'center_name', 'first_public', 'last_updated',
    'experiment_title', 'study_title', 'study_alias', 'experiment_alias',
    'run_alias', 'fastq_bytes', 'fastq_md5', 'fastq_ftp', 'fastq_aspera',
    'fastq_galaxy', 'submitted_bytes', 'submitted_md5', 'submitted_ftp',
    'submitted_aspera', 'submitted_galaxy', 'submitted_format',
    'sra_bytes', 'sra_md5', 'sra_ftp', 'sra_aspera', 'sra_galaxy',
    'cram_index_ftp', 'cram_index_aspera', 'cram_index_galaxy',
    'sample_alias', 'broker_name', 'sample_title', 'first_created'
]

def ena_search(query, limit=1000000):
    """USE ENA's API to retreieve the latest results."""
    import requests
    # ENA browser info: http://www.ebi.ac.uk/ena/about/browser
    address = 'http://www.ebi.ac.uk/ena/data/warehouse/search'
    query = (
        f'"{query} AND library_source=GENOMIC AND '
        '(library_strategy=OTHER OR library_strategy=WGS OR '
        'library_strategy=WGA) AND (library_selection=MNase OR '
        'library_selection=RANDOM OR library_selection=unspecified OR '
        'library_selection="size fractionation")"'
    )
    result = 'result=read_run'
    display = 'display=report'
    limit = f'limit={limit}'
    url = f'{ENA_URL}&query={query}&{limit}&fields={",".join(FIELDS)}'
    response = requests.get(url)
    if not response.text:
        print(f'{query} did not return any results from ENA.', file=sys.stderr)
        sys.exit(1)

    return response.text.split('\n')


def parse_accessions(results, min_read_length=None, min_base_count=None):
    """Parse Illumina experiment accessions from the ENA results."""
    accessions = []
    filtered = {'min_base_count':0, 'min_read_length':0, 'technical':0, 'filtered': []}
    for line in results:
        if line.startswith(FIELDS[0]):
            continue
        else:
            col_vals = line.split('\t')
            if len(col_vals) == len(FIELDS):
                c = dict(zip(FIELDS, col_vals))
                if c['instrument_platform'] == "ILLUMINA":
                    passes = True
                    reason = []
                    if not c['fastq_bytes']:
                        passes = False
                        reason.append(f'Missing FASTQs')
                        filtered['technical'] += 1
                    else:
                        if min_read_length:
                            total_fastqs = len(c['fastq_bytes'].rstrip(';').split(';'))
                            read_length = int(float(c['base_count']) / (float(c['read_count']) * total_fastqs))
                            if read_length < min_read_length:
                                passes = False
                                reason.append(f'Failed mean read length ({read_length} bp) filter, expected > {min_read_length} bp')
                                filtered['min_read_length'] += 1

                        if min_base_count:
                            if float(c['base_count']) < min_base_count:
                                passes = False
                                reason.append(f'Failed base count ({c["base_count"]} bp) filter, expected > {min_base_count} bp')
                                filtered['min_base_count'] += 1

                    if passes:
                        accessions.append(c['experiment_accession'])
                    else:
                        filtered['filtered'].append({
                            'accession': c['experiment_accession'],
                            'reason': ';'.join(reason)
                        })

    return [list(set(accessions)), filtered]


def parse_query(query, exact_taxon=False):
    """Return the query based on if Taxon ID or BioProject/Study accession."""
    import re
    try:
        taxon_id = int(query)
        if exact_taxon:
            return ['taxon', f'tax_eq({taxon_id})']
        else:
            return ['taxon', f'tax_tree({taxon_id})']
    except ValueError:
        # It is a accession or scientific name
        # Test Accession
        # Thanks! https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html#accession-numbers
        if re.match(r'PRJ[E|D|N][A-Z][0-9]+|[E|D|S]RP[0-9]{6,}', query):
            return ['bioproject', f'(study_accession={query} OR secondary_study_accession={query})']
        elif re.match(r'SAM(E|D|N)[A-Z]?[0-9]+|(E|D|S)RS[0-9]{6,}', query):
            return ['biosample', f'(sample_accession={query} OR secondary_sample_accession={query})']
        elif re.match(r'(E|D|S)RR[0-9]{6,}', query):
            return ['run', f'(run_accession={query})']
        else:
            # Assuming it is a scientific name
            return ['taxon', f'tax_name("{query}")']


if __name__ == '__main__':
    import argparse as ap
    import random
    import textwrap

    parser = ap.ArgumentParser(
        prog='bactopia search',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Search ENA for associated WGS samples'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} PRJNA480016 --limit 20
              {PROGRAM} 1280 --exact_taxon --limit 20'
              {PROGRAM} "staphylococcus aureus" --limit 20
              {PROGRAM} SAMN01737350
              {PROGRAM} SRR578340
        ''')
    )
    parser.add_argument('query', metavar="STR", type=str,
                        help='Taxon ID or Study, BioSample, or Run accession')
    parser.add_argument(
        '--exact_taxon', action='store_true', help='Exclude Taxon ID descendents.'
    )
    parser.add_argument(
        '--outdir', metavar="OUTPUT_DIRECTORY", type=str, default=".",
        help='Directory to write output. (Default: .)'
    )
    parser.add_argument(
        '--prefix', metavar="PREFIX", type=str, default="ena",
        help='Prefix to use for output file names. (Default: ena)'
    )
    parser.add_argument(
        '--limit', metavar="INT", type=int, default=1000000,
        help='Maximum number of results to return. (Default: 1000000)'
    )

    parser.add_argument(
        '--biosample_subset', metavar="INT", type=int, default=0,
        help='If a BioSample has multiple Experiments, pick a random subset. (Default: Return All)'
    )

    parser.add_argument(
        '--min_read_length', metavar="INT", type=int,
        help='Filters samples based on minimum mean read length. (Default: No filter)'
    )
    parser.add_argument(
        '--min_base_count', metavar="INT", type=int,
        help='Filters samples based on minimum basepair count. (Default: No filter)'
    )
    parser.add_argument(
        '--min_coverage', metavar="INT", type=int,
        help='Filter samples based on minimum coverage (requires --genome_size)'
    )
    parser.add_argument(
        '--genome_size', metavar="INT", type=int,
        help='Genome size to estimate coverage (requires --coverage)'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    min_read_length = args.min_read_length
    min_base_count = args.min_base_count
    if args.min_coverage and args.genome_size:
        if args.min_base_count:
            print("--min_base_count cannot be used with --coverage/--genome_size. Exiting...",
                  file=sys.stderr)
            sys.exit(1)
        else:
            min_base_count = args.min_coverage * args.genome_size
    elif args.min_coverage or args.genome_size:
        print("--coverage and --genome_size must be used together. Exiting...",
              file=sys.stderr)
        sys.exit(1)

    query_type, query = parse_query(args.query, exact_taxon=args.exact_taxon)
    results = ena_search(query, limit=args.limit)
    accessions, filtered = parse_accessions(results, min_read_length=min_read_length,
                                            min_base_count=min_base_count)

    WARNING_MESSAGE = None
    if query_type == 'biosample' and args.biosample_subset > 0:
        if len(accessions) > args.biosample_subset:
            WARNING_MESSAGE = f'WARNING: Selected {args.biosample_subset} Experiment accession(s) from a total of {len(accessions)}'
            accessions = random.sample(accessions, args.biosample_subset)

    # Output the results
    results_file = f'{args.outdir}/{args.prefix}-results.txt'
    with open(results_file, 'w') as output_fh:
        for result in results:
            if result:
                output_fh.write(f'{result}\n')

    accessions_file = f'{args.outdir}/{args.prefix}-accessions.txt'
    with open(accessions_file, 'w') as output_fh:
        for accession in accessions:
            output_fh.write(f'{accession}\n')

    filtered_file = f'{args.outdir}/{args.prefix}-filtered.txt'
    with open(filtered_file, 'w') as output_fh:
        output_fh.write(f'accession\treason\n')
        for f in filtered['filtered']:
            output_fh.write(f'{f["accession"]}\t{f["reason"]}\n')

    with open(f'{args.outdir}/{args.prefix}-summary.txt', 'w') as output_fh:
        output_fh.write(f'QUERY: {query}\n')
        output_fh.write(f'LIMIT: {args.limit}\n')
        output_fh.write(f'RESULTS: {len(results) - 2} ({results_file})\n')
        output_fh.write(f'ILLUMINA ACCESSIONS: {len(accessions)} ({accessions_file})\n')

        if WARNING_MESSAGE:
            output_fh.write(f'\t{WARNING_MESSAGE}\n')

        if min_read_length or min_base_count:
            output_fh.write(f'FILTERED ACCESSIONS: {len(filtered["filtered"])}\n')
            if min_read_length:
                output_fh.write(f'\tFAILED MIN READ LENGTH ({min_read_length} bp): {filtered["min_read_length"]}\n')
            if min_base_count:
                output_fh.write(f'\tFAILED MIN BASE COUNT ({min_base_count} bp): {filtered["min_base_count"]}\n')
        else:
            output_fh.write(f'FILTERED ACCESSIONS: no filters applied\n')

        output_fh.write(f'\tMISSING FASTQS: {filtered["technical"]}\n')
