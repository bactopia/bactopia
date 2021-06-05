#! /usr/bin/env python3
"""
"""
PROGRAM = "select-references"
VERSION = "1.7.1"


def use_eutils(accession):
    from Bio import Entrez
    import time
    import json
    Entrez.email = "robert.petit@emory.edu"
    Entrez.tool = "BactopiaSelectReferences"
    accession = accession.split('.')[0]
    handle = Entrez.esearch(db="assembly", term=accession, retmax="500")
    record = Entrez.read(handle, validate=False)
    time.sleep(1) # Be kind to NCBI

    handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]))
    record = Entrez.read(handle, validate=False)
    time.sleep(1) # Be kind to NCBI

    records = []
    excluded = set()
    for assembly in record['DocumentSummarySet']["DocumentSummary"]:
        if assembly["ExclFromRefSeq"]:
            # PGAP can cause some Assemblies to eventually become excluded from RefSeq
            # https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
            for reason in assembly["ExclFromRefSeq"]:
                excluded.add(reason)
        else:
            records.append(assembly["AssemblyAccession"])

    if excluded:
        return [','.join(list(excluded)), True]
    else:
        return [sorted(records, reverse=True)[0], False]


def use_http(accession):
    """
    Use urllib to get a link.
    Example GCF_001548295: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/548/295/ 

    Need to extract "<a href="GCF_001548295.1_ASM154829v1/">GCF_001548295.1_ASM154829v1/</a>"
    """
    import re
    import requests
    accession, version = accession.split('.')
    db, digits = accession.split("_")
    digits_split = '/'.join(re.findall('.{1,3}', digits))
    url = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{digits_split}'
    
    r = requests.get(url)
    current_accession = []
    if r.status_code == 200: 
        # Success
        links = re.findall("href=[\"\'](.*?)[\"\']", r.text)
        for link in links:
            if link.startswith(accession):
                t_db, t_version, t_extra = link.split("_", 2)
                current_accession.append(f"{t_db}_{t_version}")

        if len(current_accession) == 1:
            return [current_accession[0], False, None, None]
        else:
            if not len(current_accession):
                return [current_accession, False, True, "Unable to parse and accession"]
            else:
                return [sorted(current_accession, reverse=True)[0], False, None, None]
        
    else:
        return [accession, True, False, f"Accession does not exist at {url}, status code {r.status_code}"]


def check_assembly_version(accession):
    try:
        return use_eutils(accession)
    except Exception as e:
        if e.msg == "Bad Gateway":
            print("NCBI servers are down, trying fallback.", file=sys.stderr)
            current_accession, excluded, has_error, reason = use_http(accession)
            if has_error:
                print(f"Fallback failed. Reason: {reason}", file=sys.stderr)
                sys.exit(42)
            else:
                return [current_accession, excluded]
        else:
            sys.exit(1)


if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import random
    import sys
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Select references based on Mash distance'
        )
    )

    parser.add_argument(
        'mash', metavar="FILE", type=str,
        help='Text file of Mash distances.'
    )
    parser.add_argument(
        'total', metavar="INT", type=int,
        help='Total number of references to download.'
    )
    parser.add_argument(
        '--random_tie_break', action='store_true',
        help=(
            'Select random random genome on matching Mash distances. '
            '(Default: Earliest accession'
        )
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    mash_distances = defaultdict(list)
    with open(args.mash, 'rt') as mash_fh:
        for line in mash_fh:
            reference, distance = line.rstrip().split('\t')
            mash_distances[distance].append(reference)

    remaining = args.total
    for distance, references in sorted(mash_distances.items()):
        if args.random_tie_break:
            random.shuffle(references)
        else:
            references = sorted(references)

        for reference in references:
            if reference:
                print(use_http(reference))
                current_accession, excluded = check_assembly_version(reference)
                if excluded:
                    print(
                        f'Skipping {reference}, it no longer in RefSeq. Reason: {current_accession}',
                        file=sys.stderr
                    )
                else:
                    difference = False if reference == current_accession else True
                    print(f'{reference}\t{distance}\t{current_accession}\t{difference}')
                    remaining -= 1
                    if not remaining:
                        break

        if not remaining:
            break
