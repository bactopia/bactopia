# bactopia/bactopia: Changelog

## v1.2.2 bactopia/bactopia "Tropical Punches" - 2019/10/??

### `Added`
- Size of "work" directory to the execution summary
- User controlled overwrites of existing output files
- Check for unknown parameters at runtime

### `Fixed`
- `bactopia` command now explicitly states which tag to use for Nextflow run
- Version info not updated in Dockerfile and Singularity
- Duplicated QC'd FASTQs

### `Removed`
- cgmlst support in `bactopia datasets`


## v1.2.1 bactopia/bactopia "Fruit Punches" - 2019/10/17

### `Added`
- `bactopia build` to build Conda environments
- Version info pulled from nextflow.config
- Set default values resource allocations
- Documentation on new changes
- Automatic building of Conda environments, if none exist
- `--nfdir` to determine where bactopia is being run from

### `Fixed`
- Neverending typos
- `--datasets` now, not `--dataset` <-(Typo)
- path for outputing Nextflow reports
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
- Neverending typos
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
