#! /usr/bin/env python3
"""
"""
PROGRAM = "select-references"
VERSION = "1.2.2"

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

    count = 0
    for distance, references in sorted(mash_distances.items()):
        count += 1
        if len(references) > 1:
            if args.random_tie_break:
                print(f'{random.choice(references)}')
            else:
                print(f'{sorted(references)[0]}')
        else:
            print(references[0])

        if count == args.total:
            break
