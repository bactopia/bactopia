#! /usr/bin/env python3
"""
"""
import json

PROGRAM = "merge-blast-json"
VERSION = "1.5.5"

def read_json(json_file):
    json_data = None
    with open(json_file, 'rt') as json_fh:
        json_data = json.load(json_fh)
    return json_data

if __name__ == '__main__':
    import argparse as ap
    import os
    import sys
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Merge the BLAST results of multi-FASTA runs'
        )
    )

    parser.add_argument(
        'blast', metavar="FILE", type=str,
        help='Directory containing JSON files'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    jsons = [f'{args.blast}/{f.name}' for f in os.scandir(args.blast) if f.name.endswith('.json')]
    merged_json = None
    for json_file in jsons:
        json_data = read_json(json_file)
        if merged_json:
            # Bactopia uses parallel so only one fasta entry will ever be queried hence [0]
            merged_json['BlastOutput2'].append(json_data['BlastOutput2'][0])
        else:
            merged_json = json_data

    print(json.dumps(merged_json, indent=4))
