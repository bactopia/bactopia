# bactopia/bactopia: Changelog

## v1.2.4 bactopia/bactopia "Rabbit Charm" - 2019/12/20
### `Added`
- `--conda_test` to be used for conda build test
- `--skip_fastq_check` to skip check that input FASTQs meet minimum requirements
- Undocumented parameters to the usage

### `Fixed`
- snippy not working with samtools 1.10
- NXF_HOME variable is exported to the conda env location
- speed of checking if conda environments are built

## v1.2.3 bactopia/bactopia "Tropical Punches +1" - 2019/12/19

### `Added`
- select-references selects latest assembly accession version (BioPython/Entrez)
- select-references skips assemblies that have been excluded from RefSeq
- test to for paired-end related errors (e.g. different read counts)
- `--min_genome_size` and `--max_genome_size` parameter for estimated genome sizes
    - Check is also made after assembly
- `update-version.sh` improvements
- Better genome size estimates using Mash for high and low coverage sequences
- Script to update conda environments
- added `--conda_help` to be used for conda build test

### `Fixed`
- `--random_tie_break` always true
- not using latest assembly accession for ncbi-genome-download
- usage of assemblies that have been excluded from RefSeq
- allowing PE reads with different read counts to be processed (hint... they fail pretty quickly)
- failure to stop analysis of sample with low read counts
- coverage reported as 'inf'
- references to cgmlst in the setup datasets
- non-explicit patterns in publishDir
- low coverage/read errors after QC were not put in root dir
- snippy not working with samtools 1.10

## v1.2.2 bactopia/bactopia "Tropical Punches" - 2019/10/22

### `Added`
- Size of "work" directory to the execution summary
- User controlled overwrites of existing output files
- Check for unknown parameters at runtime
- FASTQ downloads from SRA (via fastq-dl and fasterq-dump) 
- Documentation updates
- Script for building containers

### `Fixed`
- `bactopia` command now explicitly states which tag to use for Nextflow run
- Version info not updated in Dockerfile and Singularity
- Duplicated QC'd FASTQs
- nextflow: docker "Memory limited without swap" error

### `Removed`
- cgmlst support in `bactopia datasets`
- setup.py left over from pre-conda config

## v1.2.1 bactopia/bactopia "Fruit Punches" - 2019/10/17

### `Added`
- `bactopia build` to build Conda environments
- Version info pulled from nextflow.config
- Set default values resource allocations
- Documentation on new changes
- Automatic building of Conda environments, if none exist
- `--nfdir` to determine where bactopia is being run from

### `Fixed`
- Never ending typos
- `--datasets` now, not `--dataset` <-(Typo)
- path for outputting Nextflow reports
- Typo in antimicrobial_resistance.sh (task.cpus not cpus)
- `--species` is now consistent between `bactopia` and `bactopia datasets`
- Bug when checking if specific species dataset exists, but no species datasets exist
- Cleaned up version update script
- Cleaned up usage

### `Removed`
- `--max_cpus` ability to limit total cores used, access to config is being deprecated in Nextflow
- `--max_cpus` since it is redudant to `--cpus` now

## v1.2.0 bactopia/bactopia "Beestinger" - 2019/10/16

### `Added`
- `--compress` to compress certain outputs, default uncompressed
- Species name check in `bactopia datasets`
- Use requests package instead of urllib3
- Added `bactopia search` to query ENA for list of Illumina accessions
- Documentation 
    - Feedback edits
    - Output overview
    - Additional program acknowledgements
    - bibtex of citations
    - missing parameters to usage
    - info for `--genome_size` parameter
    - `bactopia search` usage
    - Workflow overview
- blastdbcmd compatible seqid to assembly fasta
    - allows search for entries with sample name
- Mask low coverage regions in consensus (subs only) fasta
- Added --dry_run to build conda envs one at a time (prevent parallel issues)
- Added Singularity recipes
- Added SLURM config

### `Fixed`
- Never ending typos
- `bactopia datasets` lowercase species names not found in MLST schemas
- `bactopia version` no longer calls nextflow
- SEQUENCE_TYPE channel groups FASTQ and assembly
- MINMER_QUERY channel groups FASTQ and signature
- Ariba MLST always running with `--noclean`
- Bugs related `--compress`
- Reduced size of per-base coverage outputs
- Removed `-parse_seqids` from makeblastdb command, caused blast queries to fail
- genomeCoverageBed failing on empty BAM files

### `Removed`
- `--clean_cache` function

## v1.1.0 bactopia/bactopia "Wooden Sword +1" - 2019/09/19

### `Added`
- NCBI's amrfinder
- Dockerfile for main bactopia install
- Completed documentation!

### `Fixed`
- insertion_sequences inputs are not now grouped into single channel
- Unintended FASTQ duplication via poor publishDir pattern

## v1.0.1 bactopia/bactopia "Wooden Sword" - 2019/09/12

### `Added`
- README.md documentation

### `Fixed`
- call_variants_auto bug fixed
- documentation
- version numbers

## v1.0.0 bactopia/bactopia "Wooden Sword" - 2019/09/04
- Initial release of bactopia/bactopia
