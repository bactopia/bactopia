#! /usr/bin/env python3
"""
Read a directory and prepare a FOFN of FASTQs
"""
VERSION = "1.0.1"
PROGRAM = "bactopia prepare"

if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import glob
    import os
    import sys

    parser = ap.ArgumentParser(
        prog='bactopia prepare',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Read a directory and prepare a FOFN of FASTQs'
        )
    )
    parser.add_argument('path', metavar="STR", type=str,
                        help='Directory where FASTQ files are stored')
    parser.add_argument(
        '-e', '--ext', metavar='STR', type=str,
        default=".fastq.gz",
        help='Extension of the FASTQs. Default: .fastq.gz'
    )
    parser.add_argument(
        '-s', '--sep', metavar='STR', type=str,
        default="_",
        help='Split FASTQ name on the last occurrence of the separator. Default: _'
    )

    parser.add_argument(
        '--pattern', metavar='STR', type=str,
        default="*.fastq.gz",
        help='Glob pattern to match FASTQs. Default: *.fastq.gz'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    abspath = os.path.abspath(args.path)
    FASTQS = {}
    for fastq in glob.glob(f'{abspath}/*{args.pattern}'):
        fastq_name = os.path.basename(fastq).replace(args.ext, "")

        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(args.sep, 1)
        sample_name = split_vals[0]
        single_end = False if len(split_vals) == 2 else True
        if sample_name not in FASTQS:
            FASTQS[sample_name] = {'single_end': False, 'fastqs': []}
        FASTQS[sample_name]['fastqs'].append(fastq)
        FASTQS[sample_name]['single_end'] = single_end

    print("sample\tr1\tr2")
    errors = []
    for sample, vals in sorted(FASTQS.items()):
        fastqs = vals['fastqs']
        if len(fastqs) > 2:
            errors.append(
                f'ERROR: "{sample}" has more than two different FASTQ files, please check.'
            )
        elif len(fastqs) != 1 and vals['single_end']:
            errors.append(
                f'ERROR: "{sample}" might be single end, but multiple FASTQs matched, please check.'
            )

        fastq_string = ""
        if len(fastqs) == 1:
            fastq_string = f'{fastqs[0]}\t'
        else:
            fastq_string = "\t".join(sorted(fastqs))

        fastq_string = "\t".join(sorted(fastqs))
        print(f'{sample}\t{fastq_string}')

    for error in errors:
        print(error, file=sys.stderr)
