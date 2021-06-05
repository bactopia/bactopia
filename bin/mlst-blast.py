#! /usr/bin/env python3
"""
usage: mlst-blast.py [-h] [--cpu INT] [--quiet] [--compressed]
                     FASTA BLAST_DIR OUTPUT

Determine MLST via BLAST

positional arguments:
  FASTA         Input FASTA file to determine MLST
  BLAST_DIR     Directory where BLAST databases are stored
  OUTPUT        File to output results to

optional arguments:
  -h, --help    show this help message and exit
  --cpu INT     Number of processors to use.
  --quiet       Do not output each command.
  --compressed  Input FASTA is Gzipped.
"""
PROGRAM = "mlst-blast"
VERSION = "1.7.1"


def pipe_command(cmd_1, cmd_2, stdout=False, stderr=False, verbose=True,
                 shell=False):
    """
    Execute a single command and return STDOUT and STDERR.

    If stdout or stderr are given, output will be written to given file name.
    """
    import subprocess
    if verbose:
        print('{0} | {1}'.format(' '.join(cmd_1), ' '.join(cmd_2)))
    stdout = open(stdout, 'w') if stdout else subprocess.PIPE
    stderr = open(stderr, 'w') if stderr else subprocess.PIPE
    p1 = subprocess.Popen(cmd_1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd_2, stdin=p1.stdout, stdout=stdout, stderr=stderr)
    p1.stdout.close()
    return p2.communicate()


def blast_alleles(input_file, blast, blastn_results, num_cpu,
                  verbose=True, compressed=False):
    """Blast assembled contigs against MLST blast database."""
    from collections import OrderedDict
    import glob
    import json
    from os.path import basename, splitext

    outfmt = "6 sseqid bitscore slen length nident mismatch pident evalue"
    results = {}
    loci = {}
    perfect_matches = []
    total_loci = 0
    for tfa in sorted(glob.glob(f'{blast}/*.tfa')):
        total_loci += 1
        blastdb = splitext(tfa)[0]
        allele = basename(blastdb)
        loci[allele] = True
        blastn = pipe_command(
            ['zcat' if compressed else 'cat', input_file],
            ['blastn', '-db', blastdb, '-query', '-', '-outfmt', outfmt,
             '-max_target_seqs', '10000', '-num_threads', num_cpu,
             '-evalue', '10000', '-ungapped', '-dust', 'no',
             '-word_size', '28'], verbose=verbose
        )
        max_bitscore = 0
        top_hits = []
        not_first = False
        for hit in blastn[0].decode("utf-8").split('\n'):
            if hit:
                cols = hit.split('\t')
                if len(cols) > 1:
                    if float(cols[1]) > max_bitscore and not_first:
                        max_bitscore = float(cols[1])

                    if cols[2] == cols[3] and cols[2] == cols[4]:
                        # perfect match
                        cols.append('perfect_match')
                        top_hits.append(cols)
                        break
                    else:
                        if float(cols[1]) == max_bitscore:
                            cols.append(
                                'has_snps' if cols[2] == cols[3] else 'partial'
                            )
                            top_hits.append(cols)
                        else:
                            break
        top_hit = []
        if not top_hits:
            # Did not return a hit
            top_hit = ['0'] * 10
            top_hit[0] = '{0}.0'.format(allele)
        elif len(top_hits) == 1:
            # Had only a single top hit
            top_hit = top_hits[0]
            top_hit.append(1)
        else:
            min_allele = 1000000
            for hit in top_hits:
                allele_number = int(hit[0].split('.')[1])
                if allele_number < min_allele:
                    # Give priority to the earliest allele on record
                    min_allele = allele_number
                    top_hit = hit
            top_hit.append(len(top_hits))

        results[allele] = OrderedDict((
            ('sseqid', top_hit[0]),
            ('bitscore', top_hit[1]),
            ('slen', top_hit[2]),
            ('length', top_hit[3]),
            ('nident', top_hit[4]),
            ('mismatch', top_hit[5]),
            ('pident', top_hit[6]),
            ('evalue', top_hit[7]),
            ('match_type', top_hit[8]),
            ('shared_bitscore', top_hit[9])
        ))
        if top_hit[8] == 'perfect_match':
            perfect_matches.append(top_hit[0])

    results['ST'] = OrderedDict((
        ('st', 'ND'), ('perfect_matches', len(perfect_matches))
    ))

    # Read Profile
    profile = {}
    extra = OrderedDict()
    extra_cols = OrderedDict()
    with open(f'{blast}/profile.txt', 'r') as profile_fh:
        for line in profile_fh:
            cols = line.rstrip('\n').split('\t')
            if line.startswith('ST'):
                col_names = cols
                for col_name in col_names:
                    if col_name != "ST" and col_name not in loci:
                        extra_cols[col_name] = True
            else:
                ST = None
                alleles = []
                extra = OrderedDict()
                for i, val in enumerate(cols):
                    if val:
                        col_name = col_names[i]
                        if col_name == "ST":
                            ST = val
                        elif col_name in loci:
                            alleles.append(f'{col_name}.{val}')
                        else:
                            extra[col_name] = val
                profile[';'.join(sorted(alleles))] = {'st': ST, 'extra': extra}

    for extra_col in extra_cols:
        if extra_col not in results['ST']:
            results['ST'][extra_col] = ''

    if len(perfect_matches) == total_loci:
        pattern = ';'.join(sorted(perfect_matches))
        if pattern in profile:
            results['ST']['st'] = profile[pattern]['st']
            for extra_col in extra_cols:
                if extra_col in profile[pattern]['extra']:
                    results['ST'][extra_col] = profile[pattern]['extra'][extra_col]
        else:
            results['ST']['st'] = 'Novel'

    with open(blastn_results, 'w') as json_fh:
        json.dump(results, json_fh, indent=4, separators=(',', ': '))


if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog='mlst-blast.py',
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Determine MLST via BLAST'
    )
    parser.add_argument('fasta', metavar="FASTA", type=str,
                        help='Input FASTA file to determine MLST')
    parser.add_argument('blast', metavar="BLAST_DIR", type=str,
                        help='Directory where BLAST databases are stored')
    parser.add_argument('output', metavar="OUTPUT", type=str,
                        help='File to output results to')
    parser.add_argument('--cpu', metavar='INT', type=int, default=1,
                        help='Number of processors to use.')
    parser.add_argument('--quiet', action='store_true',
                        help='Do not output each command.')
    parser.add_argument('--compressed', action='store_true',
                        help='Input FASTA is Gzipped.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    blast_alleles(args.fasta, args.blast, args.output, str(args.cpu),
                  verbose=not args.quiet, compressed=args.compressed)
