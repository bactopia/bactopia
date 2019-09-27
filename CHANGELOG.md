# bactopia/bactopia: Changelog

## v1.1.???? bactopia/bactopia "Beestinger" - 2019/??/??

### `Added`
- `--compress` to compress certain outputs, default uncompressed
- Species name check in `bactopia datasets`
- Use requests package instead of urllib3
- Added `bactopia search` to query ENA for list of Illumina accessions
- Documentation 
    - Feedback edits
    - Output overview
    - Additional program acknowledgements
    - Added missing parameters to usage
    - Added info for `--genome_size` parameter

### `Fixed`
- Neverending typos
- `bactopia datasets` lowercase species names not for in MLST schemas
- `bactopia version` no longer calls nextflow
- Removed `--clean_cache` function
- SEQUENCE_TYPE channel groups FASTQ and assembly
- Ariba MLST always running with `--noclean`

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
