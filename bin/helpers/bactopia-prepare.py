#! /usr/bin/env python3
"""
usage: bactopia prepare [-h] [-f STR] [-a STR] [--fastq_seperator STR]
                        [--fastq_pattern STR] [--assembly_pattern STR]
                        [--long_reads] [--version]
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
  --assembly_pattern STR
                        Glob pattern to match assembly FASTAs. Default:
                        *.fna.gz
  --long_reads          Single-end reads should be treated as long reads
  --version             show program's version number and exit
"""
VERSION = "1.4.2"
PROGRAM = "bactopia prepare"

if __name__ == '__main__':
    import argparse as ap
    from collections import defaultdict
    import glob
    import os
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
        '--assembly_pattern', metavar='STR', type=str,
        default="*.fna.gz",
        help='Glob pattern to match assembly FASTAs. Default: *.fna.gz'
    )

    parser.add_argument(
        '--long_reads', action='store_true',
        help='Single-end reads should be treated as long reads'
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
    for fastq in glob.glob(f'{abspath}/*{args.fastq_pattern}'):
        fastq_name = os.path.basename(fastq).replace(args.fastq_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(args.fastq_seperator, 1)
        sample_name = split_vals[0]
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {'pe': [], 'se': [], 'assembly': []}

        if len(split_vals) == 1:
            # single-end
            SAMPLES[sample_name]['se'].append(fastq)
        else:
            # paired-end
            SAMPLES[sample_name]['pe'].append(fastq)

    # Match assemblies
    for assembly in glob.glob(f'{abspath}/*{args.assembly_pattern}'):
        sample_name = os.path.basename(assembly).replace(args.assembly_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {'pe': [], 'se': [], 'assembly': []}
        SAMPLES[sample_name]['assembly'].append(assembly)

    FOFN = []
    for sample, vals in sorted(SAMPLES.items()):
        pe_reads = vals['pe']
        se_reads = vals['se']
        assembly = vals['assembly']
        errors = []
        is_single_end = False

        # Validate everything
        if len(assembly) > 1:
            # Can't have multiple assemblies for the same sample
            errors.append(f'ERROR: "{sample}" cannot have more than two assembly FASTA, please check.')
        elif len(assembly) == 1 and (len(pe_reads) or len(se_reads)):
            # Can't have an assembly and reads for a sample
            errors.append(f'ERROR: "{sample}" cannot have assembly and sequence reads, please check.')

        if len(pe_reads) == 1:
            # PE reads must be a pair
            errors.append(f'ERROR: "{sample}" must have two paired-end FASTQ, please check.')
        elif len(pe_reads) > 2:
            # PE reads must be a pair
            errors.append(f'ERROR: "{sample}" cannot have more than two paired-end FASTQ, please check.')

        if args.long_reads:
            if not len(pe_reads) and len(se_reads):
                # Long reads must also have short PE reads 
                print(f'WARNING: "{sample}" does not have paired-end reads, treating as single-end short reads, please verify.', file=sys.stderr)
                is_single_end = True
        else:
            if len(se_reads) > 1:
                # Can't have multiple SE reads
                errors.append(f'ERROR: "{sample}" has more than two single-end FASTQs, please check.')
            elif len(pe_reads) and len(se_reads):
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

            if pe_reads:
                runtype = 'paired-end'
                r1, r2 = sorted(pe_reads)

            if se_reads:
                if args.long_reads and not is_single_end:
                    runtype = 'hybrid'
                    extra = se_reads[0]
                else:
                    runtype = 'single-end'
                    r1 = se_reads[0]

            FOFN.append([sample, runtype, r1, r2, extra])

    if FOFN:
        print('sample\truntype\tr1\tr2\textra')
        for line in FOFN:
            print('\t'.join(line))
