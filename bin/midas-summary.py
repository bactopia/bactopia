#! /usr/bin/env python3
"""
midas-uniq.py

The public database for MIDAS includes different strains for some species. This script will limit
resolution to the species-level.

For example:
species_id	count_reads	coverage	relative_abundance
Pseudomonas_aeruginosa_57148	3608	53.91975497702909	0.9606821142074781
Pseudomonas_aeruginosa_55861	170	2.206776556776557	0.03931788579252203

Becomes:
species_id	count_reads	coverage	relative_abundance
Pseudomonas_aeruginosa	3778	53.91975497702909	0.9999

The coverages are not changed, and only the highest is kept.
"""
PROGRAM = "midas-uniq"
VERSION = "3.0.0"
import sys

if __name__ == '__main__':
    import argparse as ap
    import pandas as pd

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Update the MIDAS report for species-level resolution'
        )
    )
    parser.add_argument('prefix', metavar="PREFIX", type=str,
                        help='Prefix to be used for output files')
    parser.add_argument('midas', metavar="REPORT", type=str,
                        help='An output MIDAS report')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    """
    Example MIDAS
    species_id	count_reads	coverage	relative_abundance
    Pseudomonas_aeruginosa_57148	3608	53.91975497702909	0.9606821142074781
    Pseudomonas_aeruginosa_55861	170	2.206776556776557	0.03931788579252203
    Anaerovibrio_lipolyticus_60416	0	0.0	0.0
    Bacillus_bogoriensis_60417	0	0.0	0.0
    """
    midas = pd.read_csv(args.midas, sep='\t')

    # Use "representatives" later
    midas.rename(columns={'species_id':'representatives'}, inplace=True)

    # Split to genus and species then merge them with a space
    midas['genus'] = midas['representatives'].apply(lambda x: x.split('_')[0])
    midas['species'] = midas['representatives'].apply(lambda x: x.split('_')[1])
    midas["species_id"] = midas["genus"].astype(str) + " " + midas["species"]
    midas.drop(columns=['genus', 'species'], inplace=True)

    # Reorder columns 
    midas = midas[['species_id', 'count_reads', 'coverage', 'relative_abundance', 'representatives']]

    # Group by species and aggregate results
    # representatives now includes the original species_id value for each row that was merged
    midas_uniq = midas.groupby(midas['species_id'], as_index=False).aggregate({
        'count_reads': 'sum',
        'coverage': 'max',
        'relative_abundance': 'sum',
        'representatives':  ','.join,
    })
    midas_uniq = midas_uniq.sort_values(by='relative_abundance', ascending=False)

    # round to 5 places after decimal to match Bracken output
    midas_uniq.to_csv("{0}.midas.adjusted.abundances.txt".format(args.prefix), sep='\t', float_format='%.5f', index=False)

    # Summary
    cols = [
        'sample',
        'midas_primary_species',
        'midas_primary_species_abundance',
        'midas_secondary_species',
        'midas_secondary_species_abundance'
    ]
    results = [
        args.prefix,
        midas_uniq['species_id'].iloc[0] if midas_uniq['relative_abundance'].iloc[0] >= 0.01 else "No primary abundance > 1%",
        "{0:.5f}".format(midas_uniq['relative_abundance'].iloc[0]) if midas_uniq['relative_abundance'].iloc[0] >= 0.01 else "",
        midas_uniq['species_id'].iloc[1] if midas_uniq['relative_abundance'].iloc[1] >= 0.01 else "No secondary abundance > 1%",
        "{0:.5f}".format(midas_uniq['relative_abundance'].iloc[1]) if midas_uniq['relative_abundance'].iloc[1] >= 0.01 else ""
    ]
    with open("{0}.midas.tsv".format(args.prefix), "wt") as fh_out:
        fh_out.write("{}\n".format('\t'.join(cols)))
        fh_out.write("{}\n".format('\t'.join(results)))
