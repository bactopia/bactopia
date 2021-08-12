#! /usr/bin/env python3
"""
usage: format-16s-fasta [-h] [--outdir OUTPUT_DIRECTORY] [--extension STR]
                        [--prefix STR] [--version]
                        DIR

format-16s-fasta (v1.2.4) - Remove duplicates from input FASTAs

positional arguments:
  DIR                   Directory containing 16s FASTAs

optional arguments:
  -h, --help            show this help message and exit

Helpers:
  --outdir OUTPUT_DIRECTORY
                        Directory to write output. (Default: ./)
  --extension STR       Extension to glob inputs on. (Default: fasta)
  --prefix STR          Prefix to use for output files. (Default: 16s)
  --version             show program's version number and exit

example usage:
  format-16s-fasta ./
"""
PROGRAM = "format-16s-fasta"
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

if __name__ == '__main__':
    import argparse as ap
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
        help='Directory containing 16s FASTAs'
    )

    group5 = parser.add_argument_group('Helpers')
    group5.add_argument(
        '--outdir', metavar="OUTPUT_DIRECTORY", type=str, default="./",
        help='Directory to write output. (Default: ./)'
    )
    group5.add_argument(
        '--extension', metavar="STR", type=str, default="fasta",
        help='Extension to glob inputs on. (Default: fasta)'
    )
    group5.add_argument(
        '--prefix', metavar="STR", type=str, default="16s",
        help='Prefix to use for output files. (Default: 16s)'
    )
    group5.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    seqs = {}
    matches = {}
    for fasta in glob.glob(f'{args.input_dir}/*.{args.extension}'):
        sample = basename(fasta).split('.')[0]
        match = None
        taxon = None
        fasta_seqs = read_fasta(fasta)
        i = 1
        for header, seq in read_fasta(fasta):
            if header.startswith(sample):
                if sample in seqs:
                    # Capture multiple copies
                    seqs[f'{sample}_{i}'] = seq
                    i += 1
                else:
                    seqs[sample] = seq
            else:
                taxon = header.split(';')[-1]
                match = header.split()[0]
                accession = match.split('.')[0]
                header_taxon = taxon
                if len(header_taxon.split()) > 1:
                    header_taxon = ' '.join(header_taxon.split()[0:2])
                match_header = f'{header_taxon} ({accession})'
                seqs[match_header] = seq
        matches[sample] = [match, taxon]

    # Write merged fasta
    with open(f'{args.outdir}/{args.prefix}-merged.fasta', 'w') as fasta_out:
        for header, seq in seqs.items():
            fasta_out.write(f'>{header}\n')
            fasta_out.write(f'{seq}\n')

    # Write matches (since duplicates were removed)
    with open(f'{args.outdir}/{args.prefix}-matches.txt', 'w') as match_out:
        match_out.write("sample\taccesion\ttaxon\n")
        for sample, match in matches.items():
            match_out.write(f"{sample}\t{match[0]}\t{match[1]}\n")
