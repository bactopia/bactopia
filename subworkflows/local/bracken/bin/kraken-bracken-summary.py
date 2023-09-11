#! /usr/bin/env python3
"""
kraken-bracken-summary.py
"""
PROGRAM = "kraken-bracken-summary"
VERSION = "2.2.0"
import sys

def kraken2_unclassified_count(kraken2_report):
    """
      0.23	4500	4500	U	0	unclassified
    """
    with open(kraken2_report, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split("\t")
            if cols[3] == "U":
                return int(cols[2])

def braken_root_count(bracken_report):
    """
    100.00	1976389	0	R	1	root
    """
    with open(bracken_report, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split("\t")
            if cols[3] == "R":
                return float(cols[1])

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
    parser.add_argument('kraken2_report', metavar="KRAKEN2_REPORT", type=str,
                        help='The Kraken2 report')
    parser.add_argument('bracken_report', metavar="BRACKEN_REPORT", type=str,
                        help='The BRacken updated Kraken2 report')
    parser.add_argument('bracken_abundances', metavar="BRACKEN_ABUNDANCES", type=str,
                        help='The Bracken output with abundances')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    unclassified_count = kraken2_unclassified_count(args.kraken2_report)
    total_count = unclassified_count + braken_root_count(args.bracken_report)

    """
    Read Bracken abundances
    name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
    Pseudomonas aeruginosa	287	S	591369	1309432	1900801	0.96175
    Pseudomonas sp. Y5-11	2749808	S	11682	14249	25931	0.01312
    Stutzerimonas stutzeri	316	S	5810	1135	6945	0.00351
    """
    bracken = pd.read_csv(args.bracken_abundances, sep='\t')
    bracken['fraction_total_reads'] = bracken['new_est_reads'] / total_count
    bracken = bracken.sort_values(by='fraction_total_reads', ascending=False)

    # Write top two and unclassified
    cols = [
        'sample',
        'bracken_primary_species',
        'bracken_primary_species_abundance',
        'bracken_secondary_species',
        'bracken_secondary_species_abundance',
        'bracken_unclassified_abundance'
    ]
    results = [
        args.prefix,
        bracken['name'].iloc[0] if bracken['fraction_total_reads'].iloc[0] >= 0.01 else "No primary abundance > 1%",
        "{0:.5f}".format(bracken['fraction_total_reads'].iloc[0]) if bracken['fraction_total_reads'].iloc[0] >= 0.01 else "",
        bracken['name'].iloc[1] if bracken['fraction_total_reads'].iloc[1] >= 0.01 else "No secondary abundance > 1%",
        "{0:.5f}".format(bracken['fraction_total_reads'].iloc[1]) if bracken['fraction_total_reads'].iloc[1] >= 0.01 else "",
        "{0:.5f}".format(unclassified_count / total_count)
    ]
    with open("{0}.bracken.tsv".format(args.prefix), "wt") as fh_out:
        fh_out.write("{}\n".format('\t'.join(cols)))
        fh_out.write("{}\n".format('\t'.join(results)))

    # Add unclassified to data table and re-sort
    unclassified = pd.DataFrame.from_dict({
        'name': ['unclassified'],
        'taxonomy_id': [0],
        'taxonomy_lvl': ['U'],
        'kraken_assigned_reads': [unclassified_count],
        'added_reads': [0],
        'new_est_reads': [unclassified_count],
        'fraction_total_reads': [unclassified_count / total_count]
    })
    bracken = pd.concat([bracken, unclassified], axis=0)
    bracken = bracken.sort_values(by='fraction_total_reads', ascending=False)
    bracken.to_csv("{0}.bracken.adjusted.abundances.txt".format(args.prefix), sep='\t', float_format='%.5f', index=False)
