#! /usr/bin/env python3
"""
bracken-to-excel.py
"""
PROGRAM = "bracken-to-excel"
VERSION = "2.2.2"
import sys

if __name__ == '__main__':
    import argparse as ap
    import pandas as pd

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Write Bracken abundances to an Excel file'
        )
    )

    parser.add_argument('prefix', metavar="PREFIX", type=str,
                        help='Prefix to use for output files')
    parser.add_argument('bracken_abundances', metavar="BRACKEN_ABUNDANCES", type=str,
                        help='The merged Bracken output with abundances')
    parser.add_argument('--limit', metavar='INT', type=int, default=5,
                        help='Limit the result to the top N rows.')
    parser.add_argument(
        '--include_unclassified', action='store_true',
        help='Include results for unclassified reads.'
    )
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    """
    Read Bracken abundances
    name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
    Pseudomonas aeruginosa	287	S	591369	1309432	1900801	0.96175
    Pseudomonas sp. Y5-11	2749808	S	11682	14249	25931	0.01312
    Stutzerimonas stutzeri	316	S	5810	1135	6945	0.00351
    """
    bracken = pd.read_csv(args.bracken_abundances, sep='\t')
    samples = bracken["sample"].unique()

    # Write Excel file, with each sample getting its own sheet
    with pd.ExcelWriter(f"{args.prefix}.xlsx") as writer:  
        for sample in samples:
            sheet_name = sample
            if len(sample) > 31:
                sheet_name = sample[:31]
            df = bracken[bracken["sample"] == sample]
            if not args.include_unclassified:
                df = df[df["name"] != "unclassified"]
            
            df.sort_values(by='fraction_total_reads', ascending=False).head(args.limit).to_excel(writer, sheet_name=sheet_name, index=False)
            print(sample)
