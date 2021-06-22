#! /usr/bin/env python3
"""
usage: mask-consensus [-h] [--mincov INT] [--version]
                      SAMPLE REFERENCE SUBS_FASTA SUBS_VCF COVERAGE

mask-consensus - Snippy consensus (subs) with coverage masking.

positional arguments:
  SAMPLE        Sample name
  REFERENCE     Reference name
  SUBS_FASTA    Input "consensus.subs.fa" FASTA file
  SUBS_VCF      Input ".subs.vcf" VCF file
  COVERAGE      Per-base coverage of alignment

optional arguments:
  -h, --help    show this help message and exit
  --mincov INT  Minimum required coverage to not mask.
  --version     show program's version number and exit
"""
PROGRAM = "mask-consensus"
VERSION = "1.6.0"
import sys


def read_coverage(coverage):
    """Read the per-base coverage input."""
    import re
    accession = None
    length = None
    first_line = True
    coverages = {}
    with open(coverage, 'rt') as coverage_fh:
        for line in coverage_fh:
            line = line.rstrip()
            if line.startswith('##'):
                # ##contig=<ID=NZ_CP020108,length=5407749>
                contig = re.search(r'contig=<ID=(.*),length=([0-9]+)>', line)
                if contig:
                    accession = contig.group(1)
                    length = contig.group(2)
                    coverages[accession] = {'length':int(length), 'positions': []}
                else:
                    print(f'{line} is an unexpected format.', file=sys.stderr)
                    sys.exit(1)
            else:
                if line:
                    coverages[accession]['positions'].append(int(line))

    for accession, vals in coverages.items():
        if len(vals['positions']) != vals['length']:
            print(f'Observed bases ({len(vals["positions"])} in {accession} not expected length ({vals["length"]}).', file=sys.stderr)
            sys.exit(1)

    return coverages


def read_vcf(vcf):
    """Get positions with a substitution."""
    subs = {}
    with open(vcf, 'rt') as vcf_fh:
        for line in vcf_fh:
            if not line.startswith("#"):
                line = line.split('\t')
                # 0 = accession, 1 = position
                if line[0] not in subs:
                    subs[line[0]] = {}
                subs[line[0]][line[1]] = True
    return subs


def read_fasta(fasta):
    """Parse the input FASTA file."""
    from Bio import SeqIO
    seqs = {}
    with open(fasta, 'r') as fasta_fh:
        for record in SeqIO.parse(fasta_fh,'fasta'):
            seqs[record.name] = str(record.seq)
    return seqs


def mask_sequence(sequence, coverages, subs, mincov):
    """Mask positions with low or no coverage in the input FASTA."""
    masked_seqs = {}
    
    for accession, vals in coverages.items():
        bases = []
        coverage = vals['positions']
        for i, cov in enumerate(coverage):
            if cov >= mincov:
                # Passes
                if accession in subs:
                    if str(i+1) in subs[accession]:
                        # Substitution
                        bases.append(sequence[accession][i].lower())
                    else:
                        # Same as reference
                        bases.append(sequence[accession][i])
                else:
                    # No SNPs, Same as reference
                    bases.append(sequence[accession][i])
            elif cov:
                # Low coverage
                bases.append("N")
            else:
                # 0 coverage
                bases.append('n')

        if len(bases) != len(sequence[accession]):
            print(f'Masked sequence ({len(bases)} for {accession} not expected length ({len(sequence[accession])}).',
                file=sys.stderr)
            sys.exit(1)
        else:
            masked_seqs[accession] = bases

    return masked_seqs


def format_header(sample, reference, accession, length):
    """Return a newly formatted header."""
    title = f'Pseudo-seq with called substitutions and low coverage masked'
    return f'>gnl|{accession}|{sample} {title} [assembly_accession={reference}] [length={length}]'


def chunks(s, n):
    """
    Produce `n`-character chunks from `s`.
    https://stackoverflow.com/questions/7111068/split-string-by-count-of-characters
    """
    for start in range(0, len(s), n):
        yield s[start:start+n]


if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Snippy consensus (subs) with coverage masking.'
        )
    )
    parser.add_argument('sample', metavar="SAMPLE", type=str,
                        help='Sample name')
    parser.add_argument('reference', metavar="REFERENCE", type=str,
                        help='Reference name')
    parser.add_argument('fasta', metavar="SUBS_FASTA", type=str,
                        help='Input "consensus.subs.fa" FASTA file')
    parser.add_argument('vcf', metavar="SUBS_VCF", type=str,
                        help='Input ".subs.vcf" VCF file')
    parser.add_argument('coverage', metavar="COVERAGE", type=str,
                        help='Per-base coverage of alignment')
    parser.add_argument('--mincov', metavar='INT', type=int, default=10,
                        help='Minimum required coverage to not mask.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    coverages = read_coverage(args.coverage)
    sub_positions = read_vcf(args.vcf)
    seqs = read_fasta(args.fasta)
    masked_seqs = mask_sequence(seqs, coverages, sub_positions, args.mincov)
    for accession, seq in masked_seqs.items():
        header = format_header(args.sample, args.reference, accession, len(seq))
        print(header)
        for chunk in chunks(seq, 60):
            print("".join(chunk))
