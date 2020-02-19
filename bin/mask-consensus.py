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
  COVERAGE      Directory where BLAST databases are stored

optional arguments:
  -h, --help    show this help message and exit
  --mincov INT  Minimum required coverage to not mask.
  --version     show program's version number and exit
"""
PROGRAM = "mask-consensus"
VERSION = "1.3.0"
import sys


def read_coverage(coverage):
    """Read the per-base coverage input."""
    import re
    accession = None
    length = None
    first_line = True
    coverages = []
    with open(coverage, 'rt') as coverage_fh:
        for line in coverage_fh:
            line = line.rstrip()
            if first_line:
                # ##contig=<ID=NZ_CP020108,length=5407749>
                contig = re.search(r'contig=<ID=(.*),length=([0-9]+)>', line)
                if contig:
                    accession = contig.group(1)
                    length = contig.group(2)
                else:
                    print(f'{line} is an unexpected format.', file=sys.stderr)
                    sys.exit(1)
                first_line = False
            else:
                coverages.append(int(line))

    if len(coverages) != int(length):
        print(f'Observed bases ({len(coverages)} not expected length ({length}).', file=sys.stderr)
        sys.exit(1)

    return [accession, coverages]


def read_vcf(vcf):
    """Get positions with a substitution."""
    subs = {}
    with open(vcf, 'rt') as vcf_fh:
        for line in vcf_fh:
            if not line.startswith("#"):
                subs[line.split('\t')[1]] = True
    return subs


def read_fasta(fasta):
    """Parse the input FASTA file."""
    from Bio import SeqIO
    for fasta in SeqIO.parse(open(fasta),'fasta'):
        return str(fasta.seq)


def mask_sequence(sequence, coverages, subs, mincov):
    """Mask positions with low or no coverage in the input FASTA."""
    bases = []
    for i, cov in enumerate(coverages):
        if cov >= mincov:
            # Passes
            if str(i+1) in subs:
                # Substitution
                bases.append(sequence[i].lower())
            else:
                # Same as reference
                bases.append(sequence[i])
        elif cov:
            # Low coverage
            bases.append("N")
        else:
            # 0 coverage
            bases.append('n')

    if len(bases) != len(sequence):
        print(f'Masked sequence ({len(bases)} not expected length ({len(sequence)}).',
              file=sys.stderr)
        sys.exit(1)

    return bases


def format_header(sample, reference, accession, length):
    """Return a newly formatted header."""
    title = f'Pseudo-seq with called substitutions and low coverage masked'
    return f'>gnl|{accession}|{sample} {title} assembly_accession={reference} length={length}'


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
                        help='Directory where BLAST databases are stored')
    parser.add_argument('--mincov', metavar='INT', type=int, default=10,
                        help='Minimum required coverage to not mask.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    accession, coverages = read_coverage(args.coverage)
    sub_positions = read_vcf(args.vcf)
    seq = read_fasta(args.fasta)
    masked_seq = mask_sequence(seq, coverages, sub_positions, args.mincov)
    header = format_header(args.sample, args.reference, accession, len(coverages))
    # Print output
    print(header)
    for chunk in chunks(masked_seq, 60):
        print("".join(chunk))
