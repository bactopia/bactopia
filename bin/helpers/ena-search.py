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
VERSION = "1.2.1"
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


def parse_accessions(results):
    """Parse Illumina experiment accessions from the ENA results."""
    accessions = []
    for line in results:
        if line.startswith(FIELDS[0]):
            continue
        else:
            col_vals = line.split('\t')
            if len(col_vals) == len(FIELDS):
                c = dict(zip(FIELDS, col_vals))
                if c['instrument_platform'] == "ILLUMINA":
                    accessions.append(c['experiment_accession'])
    return list(set(accessions))


def parse_query(query, exact_taxon=False):
    """Return the query based on if Taxon ID or BioProject/Study accession."""
    import re
    try:
        taxon_id = int(query)
        if exact_taxon:
            return f'tax_eq({taxon_id})'
        else:
            return f'tax_tree({taxon_id})'
    except ValueError:
        # It is a accession or scientific name
        # Test Accession
        # Thanks! https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html#accession-numbers
        if re.match(r'PRJ[E|D|N][A-Z][0-9]+|[E|D|S]RP[0-9]{6,}', query):
            return f'(study_accession={query} OR secondary_study_accession={query})'
        else:
            # Assuming it is a scientific name
            return f'tax_name("{query}")'


if __name__ == '__main__':
    import argparse as ap
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
        ''')
    )
    parser.add_argument('query', metavar="STR", type=str,
                        help='Taxon ID or Study accession')
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
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    query = parse_query(args.query, exact_taxon=args.exact_taxon)
    results = ena_search(query, limit=args.limit)
    accessions = parse_accessions(results)

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

    with open(f'{args.outdir}/{args.prefix}-summary.txt', 'w') as output_fh:
        output_fh.write(f'QUERY: {query}\n')
        output_fh.write(f'LIMIT: {args.limit}\n')
        output_fh.write(f'RESULTS: {len(results) - 2} ({results_file})\n')
        output_fh.write(f'ILLUMINA ACCESSIONS: {len(accessions)} ({accessions_file})\n')
