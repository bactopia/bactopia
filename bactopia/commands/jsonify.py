import bactopia
PROGRAM = 'bactopia jsonify'
VERSION = bactopia.__version__


def main():
    import argparse as ap
    import json
    import os
    import sys
    import textwrap
    from bactopia.parse import parse_bactopia_files

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=f'{PROGRAM} (v{VERSION}) - Create a JSON summary of Bactopia outputs',
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} SAMPLE_DIRECTORY
        ''')
    )

    parser.add_argument(
        'sample', metavar="SAMPLE_DIRECTORY", type=str,
        help='Sample directory containing Bactopia output.'
    )

    parser.add_argument(
        '--outdir', metavar="OUTPUT_DIRECTORY", type=str, default="./",
        help='Directory to write output. (Default: ./)'
    )

    parser.add_argument(
        '--prefix', metavar="STR", type=str, default="bactopia",
        help='Prefix to use for output files. (Default: bactopia)'
    )
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing reports.')
    parser.add_argument('--depends', action='store_true',
                        help='Verify dependencies are installed.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')
    parser.add_argument('--verbose', action='store_true',
                        help='Increase the verbosity of output.')
    parser.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    sample_name = os.path.basename(args.sample.rstrip("/"))
    sample_prefix = os.path.dirname(args.sample.rstrip("/"))
    json_file = f'{args.outdir}/{sample_name}.json'
    if os.path.exists(f'{args.outdir}/{sample_name}.json') and not args.force:
        print(f"Existing JSON for {sample_name} found in {args.outdir}. Will not overwirte unless --force is used. Exiting.",
              file=sys.stderr)
        sys.exit(1)
    
    with open(json_file, 'wt') as json_fh:
        json.dump(parse_bactopia_files(sample_prefix, sample_name), json_fh)


if __name__ == '__main__':
    main()
