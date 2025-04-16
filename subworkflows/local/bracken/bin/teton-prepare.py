#! /usr/bin/env python3
"""
teton-prepare.py
"""
PROGRAM = "teton-prepare"
VERSION = "3.1.0"
import sys
from pathlib import Path


def is_local(filename):
    if filename.startswith('gs://') or filename.startswith('s3://') or filename.startswith('az://') or filename.startswith('https://'):
        return False
    return True


if __name__ == '__main__':
    import argparse as ap
    import pandas as pd

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Update the Bracken abundances with unclassified counts'
        )
    )

    parser.add_argument('prefix', metavar="PREFIX", type=str,
                        help='Prefix to use for output files')
    parser.add_argument('sizemeup', metavar="SIZEMEUP_OUTPUT", type=str,
                        help='The output from sizemeup')
    parser.add_argument('run_type', metavar="RUN_TYPE", type=str,
                        help='The runtype of the sample')
    parser.add_argument('fastqs', metavar="FASTQS", type=str,
                        help='A comma separated list of fastq files')
    parser.add_argument('outdir', metavar="OUTDIR", type=str,
                        help='The output directory for the Teton')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    sample_sheet = {
        'sample': args.prefix,
        'runtype': args.run_type,
        'genome_size': 0,
        'species': "UNKNOWN_SPECIES",
        'r1': "",
        'r2': "",
        'extra': ""
    }

    """
    Example TSV sizemeup output
    we only want species and size columns

    name    tax_id  category        size    source  method
    Staphylococcus aureus   1280    bacteria        2800000 ncbi    manually-set
    """
    df = pd.read_csv(args.sizemeup, sep="\t")
    df = df[['name', 'size', 'category']]
    category = df['category'][0]
    sample_sheet['genome_size'] = str(df['size'][0])
    sample_sheet['species'] = str(df['name'][0])

    # Sort out fastqs
    fastqs = args.fastqs.split(",")
    outdir = f"{args.outdir}/{args.prefix}/teton/tools/scrubber"

    if args.run_type == "paired-end":
        sample_sheet['r1'] = f"{outdir}/{fastqs[0]}"
        sample_sheet['r2'] = f"{outdir}/{fastqs[1]}"
    elif args.run_type == "single-end" or args.run_type == "ont":
        sample_sheet['r1'] = f"{outdir}/{fastqs[0]}"
    elif args.run_type == "hybrid" or args.run_type == "short-polish":
        sample_sheet['r1'] = f"{outdir}/{fastqs[0]}"
        sample_sheet['r2'] = f"{outdir}/{fastqs[1]}"
        sample_sheet['extra'] = f"{outdir}/{fastqs[2]}"

    # Verify FASTQs exist if they are on local storage
    for key in ['r1', 'r2', 'extra']:
        if sample_sheet[key]:
            if is_local(sample_sheet[key]):
                if Path(sample_sheet[key]).exists():
                    sample_sheet[key] = Path(sample_sheet[key]).resolve()
                elif Path(f"../../../{sample_sheet[key]}").exists():
                    # --outdir may have been relative, let's adjust based on being in a Nextflow work dir
                    sample_sheet[key] = Path(f"../../../{sample_sheet[key]}").resolve()
                else:
                    print(f"Error: {sample_sheet[key]} does not exist", file=sys.stderr)
                    sys.exit(1)

    """
    Example output
    sample  runtype genome_size     species r1      r2      extra
    sample01        paired-end      0       UNKNOWN_SPECIES fastqs/sample01_R1.fastq.gz fastqs/sample01_R2.fastq.gz

    write headers then, write only bacteria samples to STDOUT
    """
    if category == "bacteria":
        with open(f"{args.prefix}.bacteria.tsv", 'w') as fh:
            print("\t".join(sample_sheet.keys()), file=fh)
            print("\t".join([str(x) for x in sample_sheet.values()]), file=fh)

        with open(f"{args.prefix}.nonbacteria.tsv", 'w') as fh:
            print("\t".join(sample_sheet.keys()), file=fh)
    else:
        with open(f"{args.prefix}.bacteria.tsv", 'w') as fh:
            print("\t".join(sample_sheet.keys()), file=fh)

        with open(f"{args.prefix}.nonbacteria.tsv", 'w') as fh:
            print("\t".join(sample_sheet.keys()), file=fh)
            print("\t".join([str(x) for x in sample_sheet.values()]), file=fh)
