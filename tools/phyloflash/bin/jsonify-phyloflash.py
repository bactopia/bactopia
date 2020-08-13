#! /usr/bin/env python3
"""
usage: jsonify-phyloflash [-h] [--version] STR

jsonify-phyloflash (v1.2.4) - Summarize the phyloFlash results

positional arguments:
  STR         The final (non-HTML) report of the phyloFlash results

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

example usage:
  jsonify-phyloflash ./
"""
PROGRAM = "jsonify-phyloflash"
VERSION = "1.4.4"


def read_phyloflash(phyloflash):
    """ Parse the phyloFlash report """
    from collections import OrderedDict
    phyloflash_results = OrderedDict()
    phyloflash_results['runtime_info'] = OrderedDict()
    RENAME_KEY = {
        'Library name': 'sample_name',
        'Forward read file': 'r1',
        'Reverse read file': 'r2',
        'Current working directory': 'cwd',
        'Minimum mapping identity': 'min_indentity',
        'Input PE-reads': 'input_reads',
        'Mapped SSU reads': 'mapped_reads',
        'Mapping ratio': 'mapped_ratio',
        'Detected median insert size': 'insert_size_median',
        'Insert size standard deviation': 'insert_size_std',
        'Used insert size:': 'insert_size_used',
        'Ratio of assembled SSU reads': 'assembled_read_ratio',
        'CPUs used': 'cpus',
        'NTUs observed once': 'observed_once',
        'NTUs observed twice': 'observed_twice',
        'NTUs observed three or more times': 'observed_three_plus',
        'NTU Chao1 richness estimate': 'chao1_estimate'
    }

    with open(phyloflash, 'rt') as phyloflash_fh:
        section = None
        status = None
        cols = None
        for line in phyloflash_fh:
            line = line.rstrip()
            if line:
                if line.startswith('List of NTUs in order of abundance'):
                    section = 'ntu_mapping'
                    status = 'header'
                    phyloflash_results[section] = OrderedDict((
                        ('description', 'List of NTUs in order of abundance (min. 3 reads mapped)'),
                        ('results', [])
                    ))
                elif line.startswith('SSU assembly based taxa'):
                    section = 'ssu_assembly'
                    status = 'header'
                    phyloflash_results[section] = OrderedDict((
                        ('description', 'SSU assembly based taxa'),
                        ('results', [])
                    ))
                elif line.startswith('Taxonomic affiliation of unassembled'):
                    section = 'unassembled'
                    status = 'result'
                    cols = ['NTU', 'reads']
                    phyloflash_results[section] = OrderedDict((
                        ('description', 'Taxonomic affiliation of unassembled reads (min. 3 reads mapped)'),
                        ('results', [])
                    ))
                elif line.startswith('---'):
                    section = None
                    status = None
                    cols = None
                elif line.startswith('phyloFlash'):
                    program, version = line.split("-")[0].rstrip().split()
                    phyloflash_results['runtime_info'] = {
                        'name': program,
                        'version': version,
                        'description': line.split("-")[1].lstrip()
                    }
                elif line.startswith('Read mapping based higher taxa'):
                    continue
                elif status == 'header':
                    cols = line.split('\t')
                    status = 'result'
                elif status == 'result':
                    phyloflash_results[section]['results'].append(dict(zip(cols, line.split('\t'))))
                else:
                    key, val = line.split('\t')
                    key = key.replace(':', '')
                    renamed_key = None
                    if key in RENAME_KEY:
                        renamed_key = RENAME_KEY[key]
                    else:
                        renamed_key = key.lower()

                    if renamed_key in ['command', 'cpus', 'r1', 'r2', 'cwd', 'min_indentity']:
                        phyloflash_results['runtime_info'][renamed_key] = val
                    else:
                        phyloflash_results[renamed_key] = val 

    return phyloflash_results


if __name__ == '__main__':
    import argparse as ap
    import json
    import sys
    import textwrap

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Summarize the phyloFlash results',
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} ./
        ''')
    )

    parser.add_argument(
        'phyloflash_report', metavar="STR", type=str,
        help='The final (non-HTML) report of the phyloFlash results'
    )
    parser.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    results = read_phyloflash(args.phyloflash_report)
    print(json.dumps(results, indent=4))
