#! /usr/bin/env python3
"""
usage: phyloflash-summary [-h] [--extension STR] [--prefix STR] [--version]
                          DIR

phyloflash-summary (v1.2.4) - Remove duplicates from input FASTAs

positional arguments:
  DIR              Directory containing phyloFlash reports (JSON format)

optional arguments:
  -h, --help       show this help message and exit

Helpers:
  --extension STR  Extension to glob inputs on. (Default: phyloFlash.json)
  --version        show program's version number and exit

example usage:
  phyloflash-summary ./
"""
PROGRAM = "phyloflash-summary"
VERSION = "1.7.1"

def read_fasta(fasta):
    """ Kudos: https://www.biostars.org/p/710/ """
    from itertools import groupby
    fasta_fh = open(fasta)
    faiter = (x[1] for x in groupby(fasta_fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

def format_taxon(taxonomy):
    """ Format the taxonomy string. """
    genus = taxonomy.split(';')[-2]
    organism = taxonomy.split(';')[-1]
    return [genus, organism]

if __name__ == '__main__':
    import argparse as ap
    import json
    import glob
    import sys
    import textwrap
    from os.path import basename
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Remove duplicates from input FASTAs',
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} ./
        ''')
    )

    parser.add_argument(
        'input_dir', metavar="DIR", type=str,
        help='Directory containing phyloFlash reports (JSON format)'
    )

    group5 = parser.add_argument_group('Helpers')
    group5.add_argument(
        '--extension', metavar="STR", type=str, default="phyloFlash.json",
        help='Extension to glob inputs on. (Default: phyloFlash.json)'
    )
    group5.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    results = {}
    col_names = ['sample', 'assembly_species', 'assembly_taxon', 'top_mapped_taxon', 'top_mapped_reads', 'chao1_estimate',
                 'assembly_hit', 'input_reads', 'mapped_reads', 'ntu_observed_once', 'ntu_observed_twice',
                 'ntu_observed_three_plus', 'message']
    print('\t'.join(col_names))
    for json_report in glob.glob(f'{args.input_dir}/*.{args.extension}'):
        message = ""
        assembly_taxon = ""

        with open(json_report, 'rt') as json_fh:
            json_data = json.load(json_fh)
            sample = json_data['sample_name']
            message = []
            results[sample] = {
                'sample': sample,
                'input_reads': json_data['input_reads'],
                'mapped_reads': json_data['mapped_reads'],
                'chao1_estimate': json_data['chao1_estimate'],
                'ntu_observed_once': json_data['observed_once'],
                'ntu_observed_twice': json_data['observed_twice'],
                'ntu_observed_three_plus': json_data['observed_three_plus'],
                'assembly_species': '',
                'assembly_taxon': '',
                'assembly_hit': '',
                'message': ''
            }

            results[sample]['top_mapped_taxon'] = json_data['ntu_mapping']['results'][0]['NTU'].split(';')[-1]
            results[sample]['top_mapped_reads'] = json_data['ntu_mapping']['results'][0]['reads']

            if 'ssu_assembly' in json_data:
                if len(json_data['ssu_assembly']['results']) > 1:
                    assembly_taxon = []
                    assembly = []
                    assembly_hit = []
                    for hit in json_data['ssu_assembly']['results']:
                        genus, organism = format_taxon(hit['taxonomy'])
                        assembly_taxon.append(genus)
                        assembly.append(organism)
                        assembly_hit.append(hit['dbHit'])

                    if len(set(assembly_taxon)) == 1:
                        results[sample]['assembly_taxon'] = ';'.join(list(set(assembly_taxon)))
                        results[sample]['assembly_species'] = ';'.join(list(set(assembly)))
                        results[sample]['assembly_hit'] = ';'.join(list(set(assembly_hit)))
                    else:
                        results[sample]['assembly_taxon'] = ';'.join(assembly_taxon)
                        results[sample]['assembly_species'] = ';'.join(assembly)
                        results[sample]['assembly_hit'] = ';'.join(assembly_hit)
                    message.append("WARNING: Multiple SSUs were assembled by SPAdes")
                else:
                    genus, organism = format_taxon(json_data['ssu_assembly']['results'][0]['taxonomy'])
                    results[sample]['assembly_taxon'] = genus
                    results[sample]['assembly_species'] = organism
                    results[sample]['assembly_hit'] = json_data['ssu_assembly']['results'][0]['dbHit']
            else:
                message.append("WARNING: missing SPAdes assembly of the SSU")

            if results[sample]['assembly_taxon']:
                if results[sample]['assembly_taxon'] != results[sample]['top_mapped_taxon'] :
                    message.append(
                        f"WARNING: Assembly taxon ({results[sample]['assembly_taxon']}) does not match "
                        f"top mapped taxon ({results[sample]['top_mapped_taxon'] })" 
                    )
            results[sample]['message'] = ';'.join(message)

    for sample, vals in sorted(results.items()):
        line = []
        for col in col_names:
            line.append(str(vals[col]))
        print('\t'.join(line))  
