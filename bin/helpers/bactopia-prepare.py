#! /usr/bin/env python3
"""
usage: bactopia prepare [-h] [-f STR] [-a STR] [--fastq_seperator STR]
                        [--fastq_pattern STR] [--pe1_pattern STR]
                        [--pe2_pattern STR] [--assembly_pattern STR] [-r]
                        [--long_reads] [--merge] [--prefix STR] [--version]
                        STR

bactopia prepare - Read a directory and prepare a FOFN of
FASTQs/FASTAs

positional arguments:
  STR                   Directory where FASTQ files are stored

optional arguments:
  -h, --help            show this help message and exit
  -f STR, --fastq_ext STR
                        Extension of the FASTQs. Default: .fastq.gz
  -a STR, --assembly_ext STR
                        Extension of the FASTA assemblies. Default: .fna.gz
  --fastq_seperator STR
                        Split FASTQ name on the last occurrence of the
                        separator. Default: _
  --fastq_pattern STR   Glob pattern to match FASTQs. Default: *.fastq.gz
  --pe1_pattern STR     Designates difference first set of paired-end reads.
                        Default: ([Aa]|[Rr]1) (R1, r1, 1, A, a)
  --pe2_pattern STR     Designates difference second set of paired-end reads.
                        Default: ([Bb]|[Rr]2) (R2, r2, 2, AB b)
  --assembly_pattern STR
                        Glob pattern to match assembly FASTAs. Default:
                        *.fna.gz
  -r, --recursive       Directories will be traversed recursively
  --long_reads          Single-end reads should be treated as long reads
  --merge               Flag samples with multiple read sets to be merged by
                        Bactopia
  --prefix STR          Replace the absolute path with a given string.
                        Default: Use absolute path
  --version             show program's version number and exit
"""
VERSION = "1.7.1"
PROGRAM = "bactopia prepare"


def search_path(path, pattern, recursive=False):
    from pathlib import Path
    if recursive:
        return Path(path).rglob(pattern)
    else:
        return Path(path).glob(pattern)


def get_path(fastq, abspath, prefix):
    fastq_path = str(fastq.absolute())
    if prefix:
        return fastq_path.replace(abspath, prefix.rstrip("/"))
    return fastq_path


if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import glob
    import os
    import re
    import sys

    parser = ap.ArgumentParser(
        prog='bactopia prepare',
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Read a directory and prepare a FOFN of FASTQs/FASTAs'
        )
    )
    parser.add_argument('path', metavar="STR", type=str,
                        help='Directory where FASTQ files are stored')
    parser.add_argument(
        '-f', '--fastq_ext', metavar='STR', type=str,
        default=".fastq.gz",
        help='Extension of the FASTQs. Default: .fastq.gz'
    )
    parser.add_argument(
        '-a', '--assembly_ext', metavar='STR', type=str,
        default=".fna.gz",
        help='Extension of the FASTA assemblies. Default: .fna.gz'
    )
    parser.add_argument(
        '--fastq_seperator', metavar='STR', type=str,
        default="_",
        help='Split FASTQ name on the last occurrence of the separator. Default: _'
    )

    parser.add_argument(
        '--fastq_pattern', metavar='STR', type=str,
        default="*.fastq.gz",
        help='Glob pattern to match FASTQs. Default: *.fastq.gz'
    )

    parser.add_argument(
        '--pe1_pattern', metavar='STR', type=str, default="[Aa]|[Rr]1",
        help='Designates difference first set of paired-end reads. Default: ([Aa]|[Rr]1) (R1, r1, 1, A, a)'
    )

    parser.add_argument(
        '--pe2_pattern', metavar='STR', type=str, default="[Bb]|[Rr]2",
        help='Designates difference second set of paired-end reads. Default: ([Bb]|[Rr]2) (R2, r2, 2, AB b)'
    )

    parser.add_argument(
        '--assembly_pattern', metavar='STR', type=str,
        default="*.fna.gz",
        help='Glob pattern to match assembly FASTAs. Default: *.fna.gz'
    )

    parser.add_argument(
        '-r', '--recursive', action='store_true',
        help='Directories will be traversed recursively'
    )

    parser.add_argument(
        '--long_reads', action='store_true',
        help='Single-end reads should be treated as long reads'
    )

    parser.add_argument(
        '--merge', action='store_true',
        help='Flag samples with multiple read sets to be merged by Bactopia'
    )

    parser.add_argument(
        '--prefix', metavar='STR', type=str,
        help='Replace the absolute path with a given string. Default: Use absolute path'
    )

    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    abspath = os.path.abspath(args.path)
    SAMPLES = {}

    # Match FASTQS
    for fastq in search_path(abspath, args.fastq_pattern, recursive=args.recursive):
        fastq_name = fastq.name.replace(args.fastq_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(args.fastq_seperator, 1)
        sample_name = split_vals[0]
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {'pe': {'r1': [], 'r2': []}, 'se': [], 'assembly': []}

        if len(split_vals) == 1:
            # single-end
            SAMPLES[sample_name]['se'].append(get_path(fastq, abspath, args.prefix))
        else:
            # paired-end
            pe1 = re.compile(args.pe1_pattern)
            pe2 = re.compile(args.pe2_pattern)
            if pe1.match(split_vals[1]):
                SAMPLES[sample_name]['pe']['r1'].append(get_path(fastq, abspath, args.prefix))
            elif pe2.match(split_vals[1]):
                SAMPLES[sample_name]['pe']['r2'].append(get_path(fastq, abspath, args.prefix))
            else:
                print(f'ERROR: Could not determine read set for "{fastq_name}".', file=sys.stderr)
                print(f'ERROR: Found {split_vals[1]} expected (R1: {args.pe1_pattern} or R2: {args.pe2_pattern})', file=sys.stderr)
                print(f'ERROR: Please use --pe1_pattern and --pe2_pattern to correct and try again.', file=sys.stderr)
                sys.exit(1)

    # Match assemblies
    for assembly in glob.glob(f'{abspath}/**/*{args.assembly_pattern}', recursive=args.recursive):
        sample_name = os.path.basename(assembly).replace(args.assembly_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {'pe': [], 'se': [], 'assembly': []}
        SAMPLES[sample_name]['assembly'].append(get_path(assembly, abspath, args.prefix))

    FOFN = []
    for sample, vals in sorted(SAMPLES.items()):
        r1_reads = vals['pe']['r1']
        r2_reads = vals['pe']['r2']
        se_reads = vals['se']
        assembly = vals['assembly']
        errors = []
        is_single_end = False
        multiple_read_sets = False
        pe_count = len(r1_reads) + len(r2_reads)

        # Validate everything
        if len(assembly) > 1:
            # Can't have multiple assemblies for the same sample
            errors.append(f'ERROR: "{sample}" cannot have more than two assembly FASTA, please check.')
        elif len(assembly) == 1 and (pe_count or len(se_reads)):
            # Can't have an assembly and reads for a sample
            errors.append(f'ERROR: "{sample}" cannot have assembly and sequence reads, please check.')

        if len(r1_reads) != len(r2_reads):
            # PE reads must be a pair
            errors.append(f'ERROR: "{sample}" must have equal paired-end read sets (R1 has {len(r1_reads)} and R2 has {len(r2_reads)}, please check.')
        elif pe_count > 2:
            # PE reads must be a pair
            if args.merge:
                multiple_read_sets = True
            else:
                errors.append(f'ERROR: "{sample}" cannot have more than two paired-end FASTQ, please check.')

        if args.long_reads:
            if not pe_count and len(se_reads):
                # Long reads must also have short PE reads
                print(f'WARNING: "{sample}" does not have paired-end reads, treating as single-end short reads, please verify.', file=sys.stderr)
                is_single_end = True
        else:
            if len(se_reads) > 1:
                # Can't have multiple SE reads
                if args.merge:
                    multiple_read_sets = True
                else:
                    errors.append(f'ERROR: "{sample}" has more than two single-end FASTQs, please check.')
            elif pe_count and len(se_reads):
                # Can't have SE and PE reads unless long reads
                errors.append(f'ERROR: "{sample}" has paired and single-end FASTQs, please check.')

        if errors:
            print('\n'.join(errors), file=sys.stderr)
        else:
            runtype = ''
            r1 = ''
            r2 = ''
            extra = ''

            if assembly:
                runtype = 'assembly'
                extra = assembly[0]

            if pe_count:
                if multiple_read_sets:
                    if args.long_reads:
                        runtype = 'hybrid-merge-pe'
                    else:
                        runtype = 'merge-pe'
                    r1 = ','.join(sorted(r1_reads))
                    r2 = ','.join(sorted(r2_reads))
                else:
                    runtype = 'paired-end'
                    r1 = r1_reads[0]
                    r2 = r2_reads[0]

            if se_reads:
                if args.long_reads and not is_single_end:
                    runtype = 'hybrid'
                    extra = se_reads[0]
                else:
                    if multiple_read_sets:
                        runtype = 'merge-se'
                        r1 = ','.join(se_reads)
                    else:
                        runtype = 'single-end'
                        r1 = se_reads[0]

            FOFN.append([sample, runtype, r1, r2, extra])

    if FOFN:
        print('sample\truntype\tr1\tr2\textra')
        for line in FOFN:
            print('\t'.join(line))
