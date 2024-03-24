#! /usr/bin/env python3
"""
Create a simple before and after report from scrubbing
"""
PROGRAM = "scrubber-summary"
VERSION = "3.0.1"
import json
import sys

def read_json(json_file):
    json_data = None
    with open(json_file, 'rt') as json_fh:
        json_data = json.load(json_fh)
    return json_data

if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Create a simple before and after report from scrubbing'
        )
    )
    parser.add_argument('sample', metavar="SAMPLE", type=str,
                        help='the name of the sample')
    parser.add_argument('original', metavar="ORIGINAL", type=str,
                        help='A summary of the original FASTQs in JSON format')
    parser.add_argument('scrubbed', metavar="SCRUBBED", type=str,
                        help='A summary of the scrubbed FASTQs in JSON format')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    original_json = read_json(args.original)
    scrubbed_json = read_json(args.scrubbed)

    # Summary
    cols = [
        'sample',
        'original_read_total',
        'scrubbed_read_total',
        'host_read_total',
    ]
    results = [
        args.sample,
        str(original_json['qc_stats']['read_total']),
        str(scrubbed_json['qc_stats']['read_total']),
        str(original_json['qc_stats']['read_total'] - scrubbed_json['qc_stats']['read_total']),
    ]

    print('\t'.join(cols))
    print('\t'.join(results))
