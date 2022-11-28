---
title: changelog
description: A full list of Bactopia releases and a description of the changes.
---
# Changelog

## v2.2.0 bactopia/bactopia "Tornado Tempo" - 2022/11/28

### `Added`
- `--dragonflye_opts` to pass additional params to dragonflye
- `--use_fastp` to use fastp for QCing reads
- `bactopia datasets` no longer depends on `ariba`
- `--skip_spell_check` in `bactopia datasets`
- updated organisms available from PubMLST
- custom profile for `wsvl`
- updated citations
- use `mambaforge` for docker builds
- Support tarball inputs for large databases
- Bactopia Tools (`bactopia --wf <NAME>`)
    - `pasty` - in silico serogrouping of _Pseudomonas aeruginosa_ isolates
    - `pbptyper` - In silico Penicillin Binding Protein (PBP) typer for _Streptococcus pneumoniae_ assemblies
    - `shigeifinder` - Serotyping Shigella and EIEC assemblies
- bump program versions in modules
    - `bakta` -> 1.5.1 
    - `bbmap` -> 39.01
    - `checkm-genome` -> 1.2.2
    - `csvtk` -> 0.25.0
    - `dragonflye` -> 1.0.13
    - `eggnog-mapper` -> 2.1.9
    - `fastani` -> 1.33
    - `fastq-scan` -> 1.0.1
    - `gamma` -> 2.2
    - `gtdbtk` -> 2.1.1
    - `iqtree` -> 2.2.0.3
    - `kleborate` -> 2.2.0
    - `mlst` -> 2.23.0
    - `mob_suite` -> 3.1.0
    - `mykrobe` -> 0.12.0
    - `nanoplot` -> 1.40.2
    - `ncbi-amrfinderplus`-> 3.10.45
    - `panaroo` -> 1.3.0
    - `pirate` -> 1.0.5
    - `rasusa` -> 0.7.0
    - `rgi` -> 6.0.1
    - `shigatyper` -> 2.0.3
    - `snpeff` -> 5.1
    - `sourmash` -> 4.5.0
    - `sra-tools` -> 3.0.0
    - `tb-profiler` -> 4.4.0

### `Fixed`
- respect short polish hybrid assembly requests
- missing variable in `bakta_download` module
- restructure to match `nf-core/modules`
- compatibility with Nextflow >=22.10
- separate versions for Linux and OSX in `assemble_genomes`
- `bactopia datasets` error messages switched to warnings

### `Enhancements to OSS`
- bump Bakta to 1.5.0 [nf-core/modules#2000](https://github.com/nf-core/modules/pull/2000)
- add recipe for pasty [bioconda/bioconda-recipes#35930](https://github.com/bioconda/bioconda-recipes/pull/35930)
- add recipe for pbptyper [bioconda/bioconda-recipes#36222](https://github.com/bioconda/bioconda-recipes/pull/36222)
- add module for pasty [nf-core/modules#2003](https://github.com/nf-core/modules/pull/2003)
- add module for pbptyper [nf-core/modules#2005](https://github.com/nf-core/modules/pull/2005)
- add module for shigeifinder [nf-core/modules#2523](https://github.com/nf-core/modules/pull/2523)

## v2.1.1 bactopia/bactopia "Dragon Twister" - 2022/07/09

### `Added`
- Add Panaroo to the `pangenome` Bactopia Tools
- Add `genotyphi` and `seroba` to Merlin
- bump `fastq-dl` to v1.1.1
- Bactopia Tools (`bactopia --wf <NAME>`)
    - `busco` - Assembly completeness based on evolutionarily informed expectations
    - `genotyphi` - Salmonella Typhi genotyping with Mykrobe outputs
    - `gubbins` - Rapid phylogenetic analysis of recombinant bacterial sequences
    - `mcroni` - Sequence variation in mobilized colistin resistance (mcr-1) genes
    - `mykrobe` - Antimicrobial resistance detection for specific species
    - `panaroo` - Pipeline for pangenome investigations
    - `plasmidfinder` - Plasmid identification from assemblies
    - `seroba` - Serotyping of Streptococcus pneumoniae from sequence reads
    - `snippy` - Rapid variant calling from Illumina sequence reads with optional core-SNP phylogeny

### `Fixed`
- improve GTDB db handling
- Properly capture assemblies with 0
- user provided adapters and phix fastas
- `--available_datasets` in the `bactopia datasets`
- name collisions in the Ariba Bactopia Tool
- bare double-quotes in Ariba reports
- updated tests for new version

### `Enhancements to OSS`
- Add database to Seroba recipe [bioconda/bioconda-recipes#35378](https://github.com/bioconda/bioconda-recipes/pull/35378)
- Update ariba dependencies [bioconda/bioconda-recipes#35383](https://github.com/bioconda/bioconda-recipes/pull/35383)
- patch pymummer recipe [bioconda/bioconda-recipes#35379](https://github.com/bioconda/bioconda-recipes/pull/35379)
- Update PlasmidFinder for better container support [bioconda/bioconda-recipes#35314](https://github.com/bioconda/bioconda-recipes/pull/35314)
- add module for ariba [nf-core/modules#1731](https://github.com/nf-core/modules/pull/1731)
- add module for mcroni [nf-core/modules#1750](https://github.com/nf-core/modules/pull/1750)
- add module for Mykrobe and GenoTyphi [nf-core/modules#1731](https://github.com/nf-core/modules/pull/1818)
- add module for plasmidfinder [nf-core/modules#1773](https://github.com/nf-core/modules/pull/1773)
- add module for Seroba [nf-core/modules#1816](https://github.com/nf-core/modules/pull/1816)
- add module for snippy/core [nf-core/modules#1855](https://github.com/nf-core/modules/pull/1855)
- let tempfile determine temp dir location [sanger-pathogens/seroba#68](https://github.com/sanger-pathogens/seroba/pull/68)
- allow the user to specify temp dir or use the system default [sanger-pathogens/pymummer#36](https://github.com/sanger-pathogens/pymummer/pull/36)

## v2.1.0 bactopia/bactopia "Silver Mist" - 2022/06/08

### `Added`
- `--use_bakta` to replace Prokka with Bakta in main workflow
- Split Bakta into two processes `download` and `run`
- added `shigatyper` to Merlin
- custom profile for `arcc`
- `bactopia prepare` error message when nothing found, and `--examples` for example use cases
- renamed `--fastqs` to `--samples` to better reflect its usage (fastqs and assemblies)
- add `--check_samples` to validate user generated FOFNs
- `--short_polish` for short read polishing for long-read assemblies
- Moved Ariba to Bactopia Tools
- bump versions in Bactopia envs
- bump Bactopia Tool versions
    - `amrfinder+` version to 3.10.30
    - `gtdb` version to 2.1.0
    - `mlst` version to 2.22.0
- Bactopia Tools (`bactopia --wf <NAME>`)
    - `ariba` - Gene identification by local assemblies
    - `gamma` - Identification, classification, and annotation of translated gene matches
    - `shigatyper`: Shigella serotype from Illumina or Oxford Nanopore reads

### `Fixed`
- @nickjhathaway Fixed assembly_level in `bactopia datasets` and auto `call_variants`
- `bactopia search` now supports ONT reads
- `vcf-annotator` unable to write to `/dev/stdout`
- typo in `bactopia download` when using singularity
- `bactopia download` will not rebuild envs its already built in a run
- `--include_genus` with `bactopia datasets` now works with multiple samples
- Parameters to copy optional species-specific files in `bactopia datasets`
- Don't use `AUTO` for IQ-TREE cpus
- `gtdb` download working again
- empty FASTQs after adapter and phix removal
- low read count error for ONT samples
- Removed `sequence_type` environment

## v2.0.3 bactopia/bactopia "King Cobra" - 2022/03/23

### `Added`
- logo and new metro style workflow
- added `test` profile for easy testing
- `mashdist` now outputs a merged TSV
- print missing required parameters in Bactopia Tools
- bump bakta to v1.4.0

### `Fixed`
- typo in `bakta` subworkflow
- `ismapper` not available in `bactopia-tools.nf`
- `kraken2` not available in `bactopia-tools.nf`
- Traits file being required for `pangenome`
- `mashtree` module not collecting all FASTA
- tests for ISMapper

## v2.0.2 bactopia/bactopia "Black Mamba" - 2022/02/22

### `Added`
- `bactopia download` which merges the `bactopia build` and `bactopia pull` commands
- modules can now use prebuilt envs
- `--available_species` to print the species in a bactopia dataset
- `--bbduk_opts` to fine tune bbduk
- set `quay` as default docker registry
- Let Nextflow handle stdout and stderr
- added `hpsuissero`, `legsta`, `sistr`, `ssuissero` to Merlin
- Bactopia Tools (`bactopia --wf <NAME>`)
    - `amrfinderplus`: Identify antimicrobial resistance in genes or proteins
    - `abricate`: Mass screening of contigs for antimicrobial and virulence genes
    - `checkm`: Assess the assembly quality of your samples
    - `hpsuissero`: Serotype prediction of Haemophilus parasuis assemblies
    - `kraken2`: Taxonomic classification of sequence reads
    - `legsta`: Typing of Legionella pneumophila assemblies
    - `mlst`: Scan contig files against PubMLST typing schemes
    - `mobsuite`: Reconstruct and annotate plasmids in bacterial assemblies
    - `rgi`: Predict antibiotic resistance from assemblies
    - `sistr`: Serovar prediction of Salmonella assemblies
    - `ssuissero`: Serotype prediction of Streptococcus suis assemblies

### `Fixed`
- `toInteger` on `null` genome size
- Index out of bounds in mlst-blast.py (caused by non allele columns)
- FOFN will fail if sample name is `null`
- improve available workflow logic
- PIRATE not exporting gene_presence_absence.csv correctly
- channel imports for `scoary` workflow
- `argument type mismatch` when using `--use_roary`
- missing channel when `agrvate` not ran as `staphtyper`
- float required options to the top when missing
- broken tests due to changes
- `bactopia datasets` including hidden directories
- Bioperl pinnings for prokka and pirate

## v2.0.1 bactopia/bactopia "Red Roc" - 2022/01/06

### `Added`
- Merlin for auto selection of species-specific tools based
- `--ask_merlin` to execute the Merlin subworkflow in Bactopia
- Schema compatible with NF-Tower
- `--skip_qc_plots` to skip FastQC and Nanoplot
- Disable CheckM by default, must use `--run_checkm`
- `--cluster_opts` for passing additional options to SLURM
- inclusion of local files (assembly, gff) in some subworkflows
- workflow for dev conda builds
- all subworkflows can be imported on other workflows
- Bactopia Tools (`bactopia --wf <NAME>`)
    - `gtdb`: Identify marker genes and assign taxonomic classifications
    - `mashdist`: Calculate Mash distances between sequences
    - `merlin`: MinmER assisted species-specific bactopia tool seLectIoN

### `Fixed`
- Resource over allocation for java related programs (-Xmx)
- ordering of skipped steps when `--datasets` not used
- logic in container selection for nf-core modules
- channel manipulation in `ncbigenomedownload`
- missing logic for `fastani`
- version captures in subworkflows
- spell check for species with `[` or `]` in their names
- version outputs for Staphtyper subworkflow
- use process names in version.yml
- `CSVTK_CONCAT` empty channels
- `options.suffix` now used as `prefix` in modules

## v2.0.0 bactopia/bactopia "Red Hawk" - 2021/12/05
With Bactopia v2 comes __a lot__ of changes! I would like to extend a huge thanks to Davi Marcon and Abhinav Sharma for their work initially converting Bactopia to DSL2. Your efforts were the momentum I needed to get the ball rolling on Bactopia v2. Thank you very much for taking your time to make such a siginificant contribution!

### `Added`
- support for Nanopore reads
- `staphopia` as a named pipeline (alias for `bactopia --wf staphopia`) for _S. aureus_ genomes
- `bactopia/bactopia-tests` repo with test data
- walkthrough for testing
- `bactopia-datasets/staphylococcus_aureus` repo with curatated _S. aureus_ datasets
- per-module testing via `pytest` (100+ tests and 7000+ outputs tested)
- per-module `meta.yml` and `params.json` for auto-building docs site
- framefork for adding new Bactopia Tools
- 19 total Bactopia Tools (`bactopia --wf <NAME>`)
    - Subworkflows (3)
        - `eggnog`: Functional annotation of proteins using orthologous groups and phylogenies
        - `pangenome`: Pangenome analysis with optional core-genome phylogeny
        - `staphtyper`: Determine the agr, spa and SCCmec types for _Staphylococcus aureus_ genomes
    - Modules (16):
        - `agrvate`: Rapid identification of _Staphylococcus aureus_ agr locus type and agr operon variants.
        - `bakta`: Rapid annotation of bacterial genomes and plasmids
        - `ectyper`: In-silico prediction of _Escherichia coli_ serotype
        - `emmtyper`: emm-typing of _Streptococcus pyogenes_ assemblies
        - `fastani`: fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
        - `hicap`: Identify cap locus serotype and structure in your _Haemophilus influenzae_ assemblies
        - `ismapper`: Identify insertion sites positions in bacterial genomes
        - `kleborate`: Screening Klebsiella genome assemblies for MLST, sub-species, and other related genes of interest
        - `lissero`: Serogroup typing prediction for _Listeria monocytogenes_
        - `mashtree`: Quickly create a tree using Mash distances
        - `meningotype`: Serotyping of _Neisseria meningitidis_
        - `ngmaster`: Multi-antigen sequence typing for _Neisseria gonorrhoeae_
        - `seqsero2`: Salmonella serotype prediction from reads or assemblies
        - `spatyper`: Computational method for finding spa types in _Staphylococcus aureus_
        - `staphopiasccmec`: Primer based SCCmec typing of _Staphylococcus aureus_ genomes
        - `tbprofiler`: Detect resistance and lineages of _Mycobacterium tuberculosis_ genomes
- Use `mamba` instead of conda for env building
- Reduced total Conda envs/Docker containers down to 7 (previously 12 not including bactopia tools)
- default to compressed outputs (`--skip_compression` to output uncompressed outputs)
- Tutorial outputs made available
- update github actions

### `Fixed`
- Cache issue causing `-resume` to fail
- amrfinder+ database not compatible error
- incorrectly parsed system memory

### Adapted from `nf-core`
- nf-core pytest setup
- nf-core/modules for bactopia tools
    - Bactopia v2 release contributed 20+ modules to nf-core/modules
- nf-core/tools arg parser
    - adapted to import params and usage based on config file

### Process Consolidation
- `makeblastdb` -> `assemble_genome`
- `call_variants`, `download_reference` -> `call_variants`
- `fastq_status`, `estiamte_genome_size` -> `gather_samples`
- `count_31mers` -> `minmer_sketch`

### `Removed`
- `bactopia tools` -> Handled by Nextflow now (`bactopia --wf <NAME>`)
- `bactopia versions` -> Program versions are output ever run now.

## v1.7.1 bactopia/bactopia "Shellbuster" - 2021/06/04
### `Added`
- bumped GTDB to v1.5.0
- added soft ulimit for `staph-typer`

### `Fixed`
- Matched PIRATE's parameter syntax for the tools
- `staph-typer` now uses GetBaseName

### 'Removed'
- PLSDB references in `bactopia datasets`

## v1.7.0 bactopia/bactopia "Chocobo Wand" - 2021/04/27
### `Added`
- Bactopia Tool `staph_typer` for agr, spa, and sccmec typing
- `--min_coverage` parameter to filter based on min coverage

### 'Removed'
- `plasmid_blast` no longer apart of main workflow

## v1.6.5 bactopia/bactopia "Z's Trident" - 2021/03/30
### `Added`
- version pins to process envs

### `Fixed`
- syntax for sourmash 4.0

## v1.6.4 bactopia/bactopia "Trident +1" - 2021/03/26
### `Added`
- added Python3.6+ to all environments

## v1.6.3 bactopia/bactopia "Trident" - 2021/03/25
### `Added`
- extra fields to `mlst-blast.py` outputs
- added Python3 to `qc_reads` environment

### `Fixed`
- rstrip on empty extra fields in mlst profile
- different BLAST+ software versions mismatch
- tbb pinnings
- `--help` and `--version` for `bactopia tools`

## v1.6.2 bactopia/bactopia "Fuscina" - 2021/03/19
### `Added`
- inputs are checked to be gzipped (this does not include FOFN)
- `--skip_amr` to skip AMRFinder+ analysis
- new  `bactopia tool` for `hicap`
- `unicycler` can be used for Illumina reads only (`--assembler unicycler`)

### `Fixed`
- AMRFinder+ software and database version mismatch
- check-fastqs.py syntax errors with prints
- `ismapper` tool processing of include/exclude files

## v1.6.1 bactopia/bactopia "Obelisk" - 2021/02/22
### `Fixed`
- sample names with "." in them breaking auto variant calling
- contig naming incompatible with GenBank

## v1.6.0 bactopia/bactopia "Harpoon" - 2021/01/22
### `Added`
- `bactopia pull` to pre-build Singularity images
- `--singularity_cache` parameter to set location of image downloads
- `--registry` to choose Docker registry to use (DockerHub, GitHub, Quay)
- `--max_downloads` sets maximum number of downloads (FASTQ/assembly) allowed at once
- `--min_time` sets the minimum amount of time a job should be given
- `bactopia search` now uses POST requests, and groups accessions into single query
- strip HTML from FASTA headers used in BLAST
- Dockerfiles now have conda.md5 label to determine if rebuild is necessary
    - MD5 is updated in Dockerfile when env is updated
- AMRFinder+ database is now provided by `bactopia datasets`
- Parameterized profile (slurm, awsbatch, etc...) variables
- `bactopia build`  
    - will retry in case of HTTP connection issues
    - `--include_tool` will build Bactopia Tool environments
- GitHub Actions
    - build Docker containers on new release (or manual trigger)
    - test that the Conda environment yamls are still valid
    - test bactopia with conda on Linux and OSX
    - test bactopia on Linux with Docker and Singularity

### `Fixed`
- redundant environment version files
- failed FASTQ/Assembly downloads no longer stop whole run
- `--max_retry` is honored now
- antimicrobial_resistance process honors `amrdir` variable
- change directory `antimicrobial_resistance` to `antimicrobial-resistance`
- rename `check_staging.py` to `check-staging.py` for consistency
- Bactopia not producing valid exit code on failure

### `Removed`
- `--containerPath` variable is replaced by `--singularity_cache`
- Native Singularity recipes, will now convert Docker to Singularity
- docs are now on bactopia.github.io repo

## v1.5.6 bactopia/bactopia "Metal Slime Earring" - 2021/01/13
### `Added`
- tweaks to the CI (via GitHub Actions)
- docker containers use quay.io now
- docker containers now use conda-pack
- `--nfconfig` will skip the conda environment build step
- input accessions are checked to be Assembly or Experiment accessions
- improved version increment script
- executor profiles can be configured by parameters

### `Fixed`
- phyloflash and download_reference environment errors
- environment path in Bactopia Tools Dockerfile and Singularity recipes
- moved version from conda yaml to conda version file
- streamlined Docker recipes
- undefined `--ftp_only` message
- typo in singularity profile
- stderr logged to file is also printed to screen
- `qc_reads` memory used now determined by base config

## v1.5.5 bactopia/bactopia "Avenger's Earring" - 2021/01/04
### `Added`
- `--prefix` option for `bactopia prepare`
- date is included in `bactopia summary` output

### `Fixed`
- removed usage personal (rpetit3) conda channel
    - aspera connect no longer supported
    - shovill-se is now used from Bioconda
- updated conda environments (phyloflash broken)

## v1.5.4 bactopia/bactopia "Nemesis Earring" - 2020/12/17
### `Fixed`
- quoted arguments being broken up (e.g. `--species "Staphylococcus aureus` -> `--species Staphylococcus`)
- `mashtree` tool failure to download with `--accessions`
- remove ncbi-genome-download header when using `--dry-run`
- undefined `name` variable in plasmid_blast
- custom work dir causing two `-w` parameters
- `bactopia search` results now contains all (illumina and non-illumina)

## v1.5.3 bactopia/bactopia "Morion Earring" - 2020/12/04
### `Added`
- Changelog moved to docs
- recursive search for `bactopia prepare`
- allow multiple FASTQs per sample (FASTQs are merged)

### `Fixed`
- unable to run `bactopia datasets` without parameters
- PLSDB blast results in invalid JSON format
- Error message for unaccepted `run_type`

## v1.5.2 bactopia/bactopia "Physical Earring" - 2020/11/18

### `Fixed`
- `--skip_qc` causing "file not found"
- `qc_reads` not honoring FINAL_BP and FINAL_READS checks

## v1.5.1 bactopia/bactopia "Astral Earring" - 2020/11/17

### `Fixed`
- `bactopia tools` not a valid project name
- `bactopia tools` `--cleanup_workdir` unrecognized variable

## v1.5.0 bactopia/bactopia "Cassie Earring" - 2020/11/12
### `Added`
- Conda environments will check if in sync with latest version now
- md5sums of all conda envs
- Verify species-specific datasets exist
- separate work dir for bactopia and bactopia tools
- `--cleanup_workdir` to delete work directory after successful run
- default values for `bactopia datasets` summary.json
- Fallback to NCBI Assembly when eUtils is down
- Additional pre-process QC checks
- OSX/Linux conda envs for Bactopia Tools
- Documentation edits and updates
- Bactopia and Nextflow versions are now output for logging purposes
- option to skip QC step (`--skip_qc`)
- `bactopia datasets` can now specify assembly level
- `bactopia tools` now use reusable conda envs
- `bactopia tools` for Roary and PIRATE can now include local assemblies

### `Fixed`
- Warn user if no completed genomes are available
- use of `--genera` for ncbi-genome-download
- improved genome_size handling
- explicit file passing for AWS Batch
- Memory estimates now floored
- PLSDB blast not being executed

## v1.4.11 bactopia/bactopia "Metamorph Ring" - 2020/09/19
### `Added`
- `bactopia build` checks if each environment is built before building
- Can specify `bactopia build` to build a specific environment
- Removed build numbers in Conda environment yamls
- Created separate Conda yamls for Linux and Mac
- NCBI assembly accessions will retrieve the most current version (e.g. .1, .2, .3, etc...)

### `Fixed`
- `bactopia datasets` trailing whitespace in species names
- `bactopia datasets` random subsample missing specified species when `--limit` and `--include_genus` used
- GitLab CI OSX compatibility
- Adaptive resource allocations
- Datasets are checked for existence
- Variant calls against references with multiple chromosomes

## v1.4.10 bactopia/bactopia "Jelly Ring" - 2020/08/25
### `Added`
- `card` is back as a default Ariba dataset
- Added timestamps to `versions` files

### `Fixed`
- `bactopia search` not creating `--outdir`
- `gtdb` tool not using prefix in outdir naming
- `pirate` tool using pangenome alignment instead of core
- Use of `scratch` causing logs to fail

## v1.4.9 bactopia/bactopia "Toreador's Ring" - 2020/08/23
### `Added`
- Support for multiple accession
    - `bactopia search` (SRA)
    - Bactopia Tools: `fastani`, `mashtree`, `pirate`, `roary` (Assembly)

### `Fixed`
- Undefined variable in `mapping_query.sh`
- ENA API endpoint for `bactopia search`
- Updated GTDB-TK to 1.3.0 to support latest downloads
- FastANI tool merge_results in no longer a separate step
    - ANI is now one-to-many calculations
- `--reassemble` misapplied

## v1.4.8 bactopia/bactopia "Shikaree Ring" - 2020/08/20
### `Added`
- Versions are logged for Bactopia
- STDOUT/STDERR logs are kept for each sample
    - can be skipped using `--skip_logs`

### `Fixed`
- Long sample names breaking Prokka annotation
- Syntax errors in Bactopia tools
- null values being tested as integers
- Ariba card and mlst downloads not working
- missing parameter in GTDB Bactopia tool

## v1.4.7 bactopia/bactopia "Serket Ring" - 2020/08/17
### `Added`
- `--no_cache` to skip caching ncbi assembly info

## v1.4.6 bactopia/bactopia "Astral Ring" - 2020/08/17
### `Added`
- Option to rebuild conda envs to default location
- Updated fastq-dl for sra-toolkit forced interaction workaround

## v1.4.5 bactopia/bactopia "Bomb Queen Ring" - 2020/08/13
### `Fixed`
- `--min_basepairs` and `--min_reads` not being honored

## v1.4.4 bactopia/bactopia "Vilma's Ring" - 2020/08/13
### `Fixed`
- `annotate_genome` name clashes
- Updated `fastq-dl` version to support new ENA API endpoint

## v1.4.3 bactopia/bactopia "Sattva Ring" - 2020/08/13
### `Added`
- `--skip_ariba` option in `bactopia datasets`

### `Fixed`
- `bactopia versions` and `bactopia citations` improper execution
- Convert spaces to tabs in citation doc
- Corrected CheckM name in program version info file
- CARD no longer default Ariba dataset download for `bactopia datasets`

## v1.4.2 bactopia/bactopia "Tamas Ring" - 2020/08/10
### `Added`
- added requirement checks of `--accessions` in `bactopia datasets`
- improved ENA spell check for species name

### `Fixed`
- file of accessions not working with `bactopia datasets`
- Dockerfile and Singularity being missed by `update-version.sh`

## v1.4.1 bactopia/bactopia "Rajas Ring" - 2020/08/06
### `Added`
- Links to publication (woohoo!)
- Can pass a Prodigal training file to `bactopia datasets`

### `Fixed`
- Typos in the Docs
- blast_primers now uses `blastn` and blast_genes uses `megablast`
- validExitStatus deprecation warning

## v1.4.0 bactopia/bactopia "Archer's Ring" - 2020/07/01
### `Added`
- New Bactopia Tools
    - `eggnog` for functional annotation using eggNOG-mapper
    - `mashtree` to create a tree using Mash distances
    - `pirate` to create pangenome using PIRATE
    - `ismapper` for insertion site discovery
    - Documentation for new tools and tweaks to existing
- BTs `roary` and `pirate` can now be run on just completed genomes
- Can limit number of completed genomes downloaded where applicable
- `bactopia datasets` can provide list of RefSeq accessions to download
- `bactopia search` can now use BioSample and Run accessions
- `bactopia search` can select a subset of Experiments associated with a BioSample
- Support for organisms with multiple MLST schemas
- Assembly QC via QUAST and CheckM
- Assemblies (local or NCBI Assembly accession) as inputs for Bactopia
- Long reads as supplementary to paired end reads for hybrid assembly
- Tools versions are locked to the minor version, not the patch

### `Fixed`
- `summary` will now determine absolute path of inputs
- `fastani` improved user reference import
- went back a version on `call_variants` openjdk version
- all bactopia tools now put nextflow info in the same folder as outputs
- Typos in docs
- Bactopia Tools now check existence of include and exclude files
- Lots more documentation
- Updated citations/tools used by Bactopia
- Did I mention typos?

### `Removed`
- ISMapper as part of the main pipeline (its now a tool)
- `insertion-sequences` in bactopia datasets

## v1.3.1 bactopia/bactopia "Emperor Hairpin" - 2020/04/20
### `Added`
- `summary` tool now gives reason for rank
- `summary` tools now splits `failed` into `exclude` and `qc-failure`
- Better documentation on how `--cpus` works in Nextflow
- Efficiency info when executed on standard profile
- split `blast_query` into `blast_genes`, `blast_primers` and `blast_proteins`
- `mapping_query` now creates multifasta of fastas at maps at once then splits per-base coverage into separate files
- `--nfconfig` users can provide their own Nextflow config file
- `fastani` users can provide their own reference now
- `bactopia versions` will print versions for tools used by Bactopia
- `bactopia citations` will print citations for tools and datasets used by Bactopia
- `bactopia search` can filter based on total bases, mean read length, and missing FASTQs
- blast queries results are only JSON format for easy parsing later
- added `--compliant` option for Prokka annotation

### `Fixed`
- build-containers.sh not working with Bactopia Tools
- Bactopia Tools container tools missing environment.yml
- Typo in `fastani` usage
- Samples with multiple QC errors counted for each error
- Incorrect ISMapper version
- typo in `summary` SLURM profile
- `gtdb` Singularity container not mounting path to GTDB database
- `roary` missing `rename` in containers
- `blast_query` results overwriting one another
- `build-containers.sh` now creates a "latest" tag
- `bactopia tool roary` outputs results based on the given prefix
- renamed `--addgenes` to `--nogenes`
- updated ASA³P citation
- Typos in Bactopia Tools docs
- Link in README.md

## v1.3.0 bactopia/bactopia "Leaping Boots" - 2020/02/19
### `Added`
- bactopia tools (BT) framework
    - docs for each tool
    - subcommand to execute tools `bactopia tools`
    - `fastani` - pairwise average nucleotide identity
    - `gtdb` - assigning objective taxonomic classifications
    - `phyloflash` - 16s assembly, alignment and tree
    - `roary` - pan-genome and core genome tree
    - `summary` - summary of results
- `--include` and `--exclude` to modify which samples used in BT analysis
- `update-version.sh` improvements
- can now set how Nextflow publishes outputs (copy, symlink, etc...) via `--publish_mode`
- Warning if output directory already exists and require `--force` to overwrite
- option (`--rfam`) to turn on ncRNA annotation in Prokka
- reduced "README.md" contents, instead point to documentation
- Updated acknowledgements and bibtex for citations

### `Fixed`
- nextflow.config version out of sync
- `--available_datasets` accessing not existent variable
- `--available_datasets` is tested before requiring inputs
- Make use of tbl2asn-forever
- adjusted how Bactopia is executed, `nextflow run` no longer pulls from github

## v1.2.4 bactopia/bactopia "Rabbit Charm" - 2019/12/20
### `Added`
- `--conda_help` to be used for conda build test
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
