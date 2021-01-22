#! /usr/bin/env python3
"""
Sometimes with AWS, files might fail to download but not cause an error.
This script checks to verify all expected inputs are staged.
"""
PROGRAM = "check-staging"
VERSION = "1.6.0"
import sys


def read_json(json_file):
    import json
    json_data = None
    with open(json_file, 'rt') as json_fh:
        json_data = json.load(json_fh)
    return json_data


def write_error(filename, error_msg):
    print(error_msg, file=sys.stderr)
    with open(filename, "wt") as fh_out:
        fh_out.write(error_msg)
    return 1


def check_reads(fq1, sample, min_reads, fq2=None):
    error = 0
    total_reads = fq1 + fq2 if fq2 else fq1

    if total_reads < min_reads:
        error_msg = (f"{sample} FASTQ(s) contain {total_reads} total reads. This does not \n"
                    f"exceed the required minimum {min_reads} read count. Further analysis is \n"
                    "discontinued.\n")
        error += write_error(f'{sample}-low-read-count-error.txt', error_msg)

    if fq2:
        if fq1 != fq2:
            # different number of reads in the pair
            error_msg = (f"{sample} FASTQs have different read counts (R1: {fq1}, R2: {fq2}). Please \n"
                        "investigate these FASTQs. Further analysis is discontinued.\n")
            error += write_error(f'{sample}-different-read-count-error.txt', error_msg)

    return error


def check_basepairs(fq1, sample, min_basepairs, fq2=None, min_proportion=None):
    error = 0
    total_bp= fq1 + fq2 if fq2 else fq1

    if total_bp < min_basepairs:
        error_msg = (f"{sample} FASTQ(s) contain {total_bp} total basepairs. This does not \n"
                    f"exceed the required minimum {min_basepairs} bp. Further analysis is \n"
                    "discontinued.\n")
        error += write_error(f'{sample}-low-sequence-depth-error.txt', error_msg)
            
    if fq2:
        proportion = float(fq1) / float(fq2) if fq1 < fq2 else float(fq2) / float(fq1)
        if proportion < min_proportion:
            # More basepairs in one sample that exceeds minimum proportion
            error_msg = (f"{sample} FASTQs failed to meet the minimum shared basepairs ({min_proportion}). \n"
                        f"They shared {proportion:.4f} basepairs, with R1 having {fq1} bp and \n"
                        f"R2 having {fq2} bp. Further analysis is discontinued.\n")
            error += write_error(f'{sample}-low-basepair-proportion-error.txt', error_msg)

    return error


if __name__ == '__main__':
    import argparse as ap
    import os
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Verifies inputs for a process are available.'
        )
    )

    parser.add_argument('--sample', metavar="STR", type=str, help='Name of the input sample.')
    parser.add_argument('--fq1', metavar="STR", type=str, help='Stats for SE or R1 FASTQ in JSON format.')
    parser.add_argument('--fq2', metavar="STR", type=str, help='Stats for R2 FASTQ in JSON format.')
    parser.add_argument('--min_proportion', metavar="FLOAT", type=float, 
                        help='The proportion of sequenced basepairs that R1 and R2 must be')
    parser.add_argument('--min_reads', metavar="INT", type=int, help='Minimum number of reads.')
    parser.add_argument('--min_basepairs',metavar="INT", type=int, help='Minimum number of seqeunced basepairs')
    parser.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    
    error = 0
    if args.fq1 and args.fq2:
        # Paired end
        r1 = read_json(args.fq1)
        r2 = read_json(args.fq2)
        error += check_reads(r1["qc_stats"]["read_total"], args.sample, args.min_reads, 
                             fq2=r2["qc_stats"]["read_total"])
        error += check_basepairs(r1["qc_stats"]["total_bp"], args.sample, args.min_basepairs, 
                                 fq2=r2["qc_stats"]["total_bp"], min_proportion=args.min_proportion)
        
    else:
        se = read_json(args.fq1)
        error += check_reads(se["qc_stats"]["read_total"], args.sample, args.min_reads)
        error += check_basepairs(se["qc_stats"]["total_bp"], args.sample, args.min_basepairs)

    sys.exit(error)
