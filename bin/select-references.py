#! /usr/bin/env python3
"""
"""
PROGRAM = "select-references"
VERSION = "1.2.2"

def check_assembly_version(accession):
    from Bio import Entrez
    import time
    Entrez.email = "robert.petit@emory.edu"
    Entrez.tool = "BactopiaSelectReferences"

    handle = Entrez.esearch(db="assembly", term=accession.split(".")[0], retmax="500")
    record = Entrez.read(handle)
    time.sleep(1) # Be kind to NCBI

    handle = Entrez.esummary(db="assembly", id=",".join(record["IdList"]))
    record = Entrez.read(handle)
    time.sleep(1) # Be kind to NCBI

    records = []
    for assembly in record['DocumentSummarySet']["DocumentSummary"]:
        records.append(assembly["AssemblyAccession"])

    return sorted(records, reverse=True)[0]


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
                current_accession = check_assembly_version(reference)
                difference = False if reference == current_accession else True
                print(f'{reference}\t{distance}\t{current_accession}\t{difference}')
                remaining -= 1
                if not remaining:
                    break

        if not remaining:
            break
