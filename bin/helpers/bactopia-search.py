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
import os
import sys
VERSION = "1.7.1"
PROGRAM = "bactopia search"
ENA_URL = ('https://www.ebi.ac.uk/ena/portal/api/search')
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


def ena_search(query, is_accession, limit=1000000):
    """USE ENA's API to retreieve the latest results."""
    import requests
    import time

    # ENA browser info: http://www.ebi.ac.uk/ena/about/browser
    query_original = query
    data = {
        'dataPortal': 'ena',
        'dccDataOnly': 'false',
        'download': 'false',
        'result': 'read_run',
        'format': 'tsv',
        'limit': limit,
        'fields': ",".join(FIELDS)
    }
    if is_accession:
        data['includeAccessions'] = query
    else:
        data['query'] = (
            f'"{query} AND library_source=GENOMIC AND '
            '(library_strategy=OTHER OR library_strategy=WGS OR '
            'library_strategy=WGA) AND (library_selection=MNase OR '
            'library_selection=RANDOM OR library_selection=unspecified OR '
            'library_selection="size fractionation")"'
        )

    headers = {'accept': '*/*', 'Content-type': 'application/x-www-form-urlencoded'}
    response = requests.post(ENA_URL, headers=headers, data=data)
    time.sleep(1)
    if not response.text:
        print(f'WARNING: {query_original} did not return any results from ENA.', file=sys.stderr)
        return [[], []]
    else:
        results = response.text.rstrip().split('\n')
        return [results[0], results[1:]]


def parse_accessions(results, min_read_length=None, min_base_count=None):
    """Parse Illumina experiment accessions from the ENA results."""
    accessions = []
    filtered = {'min_base_count': 0, 'min_read_length': 0, 'technical': 0, 'filtered': []}
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


def is_biosample(accession):
    """Check if input accession is a BioSample."""
    import re
    if re.match(r'SAM(E|D|N)[A-Z]?[0-9]+|(E|D|S)RS[0-9]{6,}', accession):
        return True
    return False


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks?page=1&tab=votes#tab-top
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def parse_query(q, accession_limit, exact_taxon=False):
    """Return the query based on if Taxon ID or BioProject/Study accession."""
    import re
    queries = []
    if os.path.exists(q):
        with open(q, 'r') as handle:
            for line in handle:
                line = line.rstrip()
                if line:
                    queries.append(line)
    elif "," in q:
        queries = q.split(',')
    else:
        queries.append(q)

    results = []
    accessions = []

    for query in queries:
        try:
            taxon_id = int(query)
            if exact_taxon:
                results.append(['taxon', f'tax_eq({taxon_id})'])
            else:
                results.append(['taxon', f'tax_tree({taxon_id})'])
        except ValueError:
            # It is a accession or scientific name
            # Test Accession
            # Thanks! https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html#accession-numbers
            if re.match(r'PRJ[E|D|N][A-Z][0-9]+|[E|D|S]RP[0-9]{6,}', query):
                accessions.append(query)
            elif re.match(r'SAM(E|D|N)[A-Z]?[0-9]+|(E|D|S)RS[0-9]{6,}', query):
                results.append(['biosample', f'(sample_accession={query} OR secondary_sample_accession={query})'])
            elif re.match(r'(E|D|S)RR[0-9]{6,}', query):
                accessions.append(query)
            else:
                # Assuming it is a scientific name
                results.append(['taxon', f'tax_name("{query}")'])

    # Split the accessions into set number
    for chunk in chunks(accessions, accession_limit):
        results.append(['accession', ','.join(chunk)])

    return results


if __name__ == '__main__':
    import argparse as ap
    import datetime
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
              {PROGRAM} SAMN01737350,SRR578340
              {PROGRAM} accessions.txt
        ''')
    )
    parser.add_argument('query', metavar="STR", type=str,
                        help=('Taxon ID or Study, BioSample, or Run accession (can also be comma '
                              'separated or a file of accessions)')
    )
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
        help='Maximum number of results (per query) to return. (Default: 1000000)'
    )
    parser.add_argument(
        '--accession_limit', metavar="INT", type=int, default=5000,
        help='Maximum number of accessions to query at once. (Default: 5000)'
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
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

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
    
    if args.biosample_subset > 0:
        if not is_biosample(args.query):
            print("--biosample_subset requires a single BioSample. Input query: {args.query} is not a BioSample. Exiting...",
                  file=sys.stderr)
            sys.exit(1)

    today = datetime.datetime.now().replace(microsecond=0).isoformat()
    results = []
    result_header = None
    accessions = []
    filtered = {'min_base_count':0, 'min_read_length':0, 'technical':0, 'filtered': {}}
    summary = []
    queries = parse_query(args.query, args.accession_limit, exact_taxon=args.exact_taxon)
    i = 1
    results_file = f'{args.outdir}/{args.prefix}-results.txt'
    accessions_file = f'{args.outdir}/{args.prefix}-accessions.txt'
    filtered_file = f'{args.outdir}/{args.prefix}-filtered.txt'
    for query_type, query in queries:
        is_accession = True if query_type == 'accession' else False
        query_header, query_results = ena_search(query, is_accession, limit=args.limit)
        results = list(set(results + query_results))
        if not result_header:
            result_header = query_header
        query_accessions, query_filtered = parse_accessions(query_results, min_read_length=min_read_length,
                                                            min_base_count=min_base_count)
        if len(query_accessions):
            WARNING_MESSAGE = None
            if query_type == 'biosample' and args.biosample_subset > 0:
                if len(query_accessions) > args.biosample_subset:
                    WARNING_MESSAGE = f'WARNING: Selected {args.biosample_subset} Experiment accession(s) from a total of {len(query_accessions)}'
                    query_accessions = random.sample(query_accessions, args.biosample_subset)
            accessions = list(set(accessions + query_accessions))
            filtered['min_base_count'] += query_filtered['min_base_count']
            filtered['min_read_length'] += query_filtered['min_read_length']
            filtered['technical'] += query_filtered['technical']
            for filtered_sample in query_filtered['filtered']:
                filtered['filtered'][filtered_sample['accession']] = filtered_sample['reason']
        else:
            if query_results:
                WARNING_MESSAGE = f'WARNING: {query} did not return any Illumina results from ENA.'
            else:
                WARNING_MESSAGE = f'WARNING: {query} did not return any results from ENA.'

        # Create Summary
        query_string = query
        if query_type == 'accession':
            total_accessions = len(query.split(','))
            if total_accessions > 5:
                query_string = f"{total_accessions} accessions were queried"
            else:
                query_string = query
        if len(queries) > 1:
            summary.append(f'QUERY ({i} of {len(queries)}): {query_string}')
            i += 1
        else:
            summary.append(f'QUERY: {query_string}')
        summary.append(f'DATE: {today}')
        summary.append(f'LIMIT: {args.limit}')
        summary.append(f'RESULTS: {len(query_results)} ({results_file})')
        summary.append(f'ILLUMINA ACCESSIONS: {len(query_accessions)} ({accessions_file})')

        if WARNING_MESSAGE:
            summary.append(f'\t{WARNING_MESSAGE}')

        if min_read_length or min_base_count:
            summary.append(f'FILTERED ACCESSIONS: {len(filtered["filtered"])}')
            if min_read_length:
                summary.append(f'\tFAILED MIN READ LENGTH ({min_read_length} bp): {query_filtered["min_read_length"]}')
            if min_base_count:
                summary.append(f'\tFAILED MIN BASE COUNT ({min_base_count} bp): {query_filtered["min_base_count"]}')
        else:
            summary.append(f'FILTERED ACCESSIONS: no filters applied')

        summary.append(f'\tMISSING FASTQS: {filtered["technical"]}')
        summary.append("")

    # Output the results
    with open(results_file, 'w') as output_fh:
        output_fh.write(f'{result_header}\n')
        for result in results:
            if result:
                output_fh.write(f'{result}\n')

    with open(accessions_file, 'w') as output_fh:
        for accession in accessions:
            output_fh.write(f'{accession}\n')

    with open(filtered_file, 'w') as output_fh:
        output_fh.write(f'accession\treason\n')
        for accession, reason in filtered['filtered'].items():
            output_fh.write(f'{accession}\t{reason}\n')

    with open(f'{args.outdir}/{args.prefix}-summary.txt', 'w') as output_fh:
        output_fh.write('\n'.join(summary))
