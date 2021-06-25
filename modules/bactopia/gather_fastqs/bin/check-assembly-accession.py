#! /usr/bin/env python3
"""
"""
PROGRAM = "check-assembly-accession"
VERSION = "1.6.0"


def check_assembly_version(accession):
    from Bio import Entrez
    import time
    import json
    Entrez.email = "robert.petit@emory.edu"
    Entrez.tool = "BactopiaCheckAssemblyAccession"

    handle = Entrez.esearch(db="assembly", term=accession, retmax="500")
    record = Entrez.read(handle, validate=False)
    time.sleep(1)  # Be kind to NCBI

    if len(record["IdList"]):
        handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]))
        record = Entrez.read(handle, validate=False)

        time.sleep(1)  # Be kind to NCBI

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
    else:

        return [f'No records found for {accession}', True]


if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import random
    import sys
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Verifies NCBI Assembly accession is latest and still available'
        )
    )

    parser.add_argument(
        'reference', metavar="STR", type=str,
        help='NCBI Assembly accession to be tested.'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    reference = args.reference.split('.')[0]
    current_accession, excluded = check_assembly_version(reference)
    if excluded:
        print(
            f'Skipping {reference}. Reason: {current_accession}',
            file=sys.stderr
        )
    else:
        print(f'Using {current_accession} for {args.reference}', file=sys.stderr)
        print(current_accession)
