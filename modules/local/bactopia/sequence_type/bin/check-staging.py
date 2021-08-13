#! /usr/bin/env python3
"""
Sometimes with AWS, files might fail to download but not cause an error.
This script checks to verify all expected inputs are staged.
"""
PROGRAM = "check-staging"
VERSION = "1.6.0"


if __name__ == '__main__':
    import argparse as ap
    import os
    import sys
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Verifies inputs for a process are available.'
        )
    )

    parser.add_argument('--fq1', metavar="STR", type=str, help='Either SE or R1 Fastq.')
    parser.add_argument('--fq2', metavar="STR", type=str, help='Either SE or R1 Fastq.')
    parser.add_argument('--extra', metavar="STR", type=str, help='Extra files')
    parser.add_argument('--genome_size', metavar="STR", type=str, help='Genome size text file')
    parser.add_argument('--assembly', metavar="STR", type=str, help='Genome assembly.')
    parser.add_argument('--is_single', action='store_true', help='Input FASTQ is single end')
    parser.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    
    if not args.is_single and args.fq2 == "null":
        # This is an issue, both files are not present
        sys.exit(80)

    if args.fq1:
        if not os.path.exists(args.fq1):
            sys.exit(81)

    if args.fq2:
        if not os.path.exists(args.fq2):
            sys.exit(82)

    if args.extra:
        if args.extra != "empty.fna.gz":
            if not os.path.exists(args.extra):
                sys.exit(90)

    if args.genome_size:
        if not os.path.exists(args.genome_size):
            sys.exit(91)

    if args.assembly:
        if not os.path.exists(args.assembly):
            sys.exit(92)
