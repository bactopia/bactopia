[![Anaconda-Server Badge](https://anaconda.org/bioconda/bactopia/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/bactopia/badges/downloads.svg)](https://anaconda.org/bioconda/bactopia)
[![Gitter](https://badges.gitter.im/bactopia/bactopia.svg)](https://gitter.im/bactopia/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# Overview

Bactopia is an extensive workflow for processing Illumina sequencing of bacterial genomes. The goal of Bactopia is process your data with a broad set of tools, so that you can get to the fun part of analyses quicker! 

Bactopia was inspired by [Staphopia](https://staphopia.emory.edu/), a workflow we (Tim Read and myself) released that targets *Staphylococcus aureus* genomes.
Using what we learned from Staphopia and user feedback, Bactopia was developed from scratch with usability, portability, and speed in mind from the start.

Bactopia uses [Nextflow](https://www.nextflow.io/) to manage the workflow, allowing for support of many types of environments (e.g. cluster or cloud). Bactopia allows for the usage of many public datasets as well as your own datasets to further enhance the analysis of your seqeuncing. Bactopia only uses software packages available from
[Bioconda](https://bioconda.github.io/) (or other
[Anaconda channels](https://anaconda.org/)) to make installation
as simple as possible for *all* users.

# Documentation
Documentation for Bactopia is available at https://bactopia.github.io/. The documentation includes much of the information here, but also a tutorial replicating [Staphopia](https://staphopia.emory.edu) and a complete overview of the usage.

# Quick Start
```
conda create -y -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
bactopia datasets datasets

# Paired-end
bactopia --R1 ${SAMPLE}_R1.fastq.gz --R2 ${SAMPLE}_R2.fastq.gz --sample ${SAMPLE} \
         --dataset datasets/ --outdir ${OUTDIR}

# Single-End
bactopia --SE ${SAMPLE}.fastq.gz --sample ${SAMPLE} --dataset datasets/ --outdir ${OUTDIR}

# Multiple Samples
bactopia prepare directory-of-fastqs/ > fastqs.txt
bactopia --fastqs fastqs.txt --dataset datasets --outdir ${OUTDIR}

# Single ENA/SRA Experiment
bactopia --accession SRX000000 --dataset datasets --outdir ${OUTDIR}

# Multiple ENA/SRA Experiments
bactopia --accessions accessions.txt --dataset datasets --outdir ${OUTDIR}
```

# Installation
Bactopia has a **a lot** of tools built into its workflow. As you can imagine, all these tools lead to numerous dependencies, and navigating dependencies can often turn into a very frustrating process. With this in mind, from the onset Bactopia was developed to only include programs that are installable using [Conda](https://conda.io/en/latest/).

Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. In other words, it makes it super easy to get the tools you need installed! The [official Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is a good starting point for getting started with Conda. Bactopia has been tested using the
[Miniconda installer](https://conda.io/en/latest/miniconda.html),
but the [Anaconda installer](https://www.anaconda.com/distribution/) should
work the same.

Once you have Conda all set up, you are ready to create an environment for
Bactopia. To do so, you can use the following command:

```
conda create -n bactopia -c conda-forge -c bioconda bactopia
```

After a few minutes you will have a new conda environment suitably named *bactopia*. To activate this environment, you will can use the following command:

```
conda activate bactopia
```

And voilà, you are all set to get started processing your data!

# Basic Usage For Bactopia
Bactopia is a wrapper around many different tools. Each of these tools may (or may not) have there own configurable parameters for you to tweak. In order to facilitate getting started with Bactopia, this section has been limited to discussion of only a few parameters. However, if you are interested in the full list of configurable parameters in Bactopia, please check out the [Complete Usage](docs/usage-complete.md) section.

## Usage
```
bactopia

Required Parameters:
    ### For Procesessing Multiple Samples
    --fastqs STR            An input file containing the sample name and
                                absolute paths to FASTQs to process

    ### For Processing A Single Sample
    --R1 STR                First set of reads for paired end in compressed (gzip)
                                FASTQ format

    --R2 STR                Second set of reads for paired end in compressed (gzip)
                                FASTQ format

    --SE STR                Single end set of reads in compressed (gzip) FASTQ format

    --sample STR            The name of the input sequences

    ### For Downloading from ENA
    --accessions            An input file containing ENA/SRA experiement accessions to
                                be processed

    --accession             A single ENA/SRA Experiment accession to be processed


Dataset Parameters:
    --datasets DIR          The path to available datasets that have
                                already been set up

    --species STR           Determines which species-specific dataset to
                                use for the input sequencing

Optional Parameters:
    --coverage INT          Reduce samples to a given coverage
                                Default: 100x

    --genome_size INT       Expected genome size (bp) for all samples
                                Default: Mash Estimate

    --outdir DIR            Directory to write results to
                                Default .

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --max_cpus INT          The maximum number of processors this workflow
                                should have access to at any given moment
                                Default: 1

    --cpus INT              Number of processors made available to a single
                                process. If greater than "--max_cpus" it
                                will be set equal to "--max_cpus"
                                Default: 1

Useful Parameters:
    --available_datasets    Print a list of available datasets found based
                                on location given by "--datasets"

    --example_fastqs        Print example of expected input for FASTQs file

    --check_fastqs          Verify "--fastqs" produces the expected inputs

    --clean_cache           Removes 'work' and '.nextflow' logs. Caution, if used,
                                the Nextflow run cannot be resumed.

    --keep_all_files        Keeps all analysis files created. By default, intermediate
                                files are removed. This will not affect the ability
                                to resume Nextflow runs, and only occurs at the end
                                of the process.

    --version               Print workflow version information

    --help                  Show this message and exit

    --help_all              Show a complete list of adjustable parameters
```

## FASTQ Inputs
Bactopia has multiple approaches to specify your input sequences. You can make use of your local FASTQs or download FASTQs from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena). Which approach really depends on what you need to achieve! The following sections describe methods to process single samples, multiple samples, downloading samples from the ENA.

### Local
#### Single Sample
When you only need to process a single sample at a time, Bactopia allows that! You only have to the sample name (`--sample`) and the whether the read set is paired-end (`--R1` and `--R2`) or a single-end (`--SE`). 

##### Use --R1, --R2 for Paired-End FASTQs
`bactopia --sample my-sample --R1 /path/to/my-sample_R1.fastq.gz --R2 /path/to/my-sample_R2.fastq.gz`

##### Use --SE for Single-End FASTQs
`bactopia --sample my-sample --SE /path/to/my-sample.fastq.gz`

#### Multiple Samples
For multiple samples, you must create a file with information about the inputs, a *file of filenames* (FOFN). This file specifies sample names and location of FASTQs to be processed. Using this information, paired-end or single-end information can be extracted as well as naming output files.

While this is an additional step for you, the user, it helps to avoid potential pattern matching errors. 

Most importantly, by taking this approach, you can process hundreds of samples in a single command. There is also the added benefit of knowing which FASTQs were analysed and their location at a later time!

##### Use --fastqs for Multiple Samples
`bactopia --fastqs my-samples.txt`

##### The FOFN Format
You can use the `--example_fastqs` to get an example of the expected structure for the input FASTQs FOFN.

```
bactopia --example_fastqs
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [naughty_borg] - revision: 0416ba407c
Printing example input for "--fastqs"

sample  r1      r2
test001 /path/to/fastqs/test_R1.fastq.gz        /path/to/fastqs/test_R2.fastq.gz
test002 /path/to/fastqs/test.fastq.gz
```

The expected structure is a **tab-delimited** table with three columns:

1. `sample`: A unique prefix, or unique name, to be used for naming output files
2. `r1`: If paired-end, the first pair of reads, else the single-end reads
3. `r2`: If paired-end, the second pair of reads

These three columns are used as the header for the file. In other words, all input FOFNs require their first line to be:
```
sample  r1      r2
```

All lines after the header line, contain unique sample names and location(s) to associated FASTQ file(s). Absolute paths should be used to prevent any *file not found* errors due to the relative path changing.

In the example above, two samples would be processed by Bactopia. Sample `test001` has two FASTQs and would be processed as pair-end reads. While sample `test002` only has a single FASTQ and would be processed as single-end reads.

##### Generating A FOFN
A script named `prepare-fofn` has been included to help aid (hopefully!) the process of creating a FOFN for your samples. This script will attempt to find FASTQ files in a given directory and output the expected FOFN format. It will also output any potential issues associated with the pattern matching.

This is currently an experimental function. There are likely bugs to be ironed out. Please be sure to give the resulting FOFN a quick look over.

###### Usage
```
bactopia prepare [-h] [-e STR] [-s STR] [--pattern STR] [--version] STR

bactopia prepare - Read a directory and prepare a FOFN of FASTQs

positional arguments:
  STR                Directory where FASTQ files are stored

optional arguments:
  -h, --help         show this help message and exit
  -e STR, --ext STR  Extension of the FASTQs. Default: .fastq.gz
  -s STR, --sep STR  Split FASTQ name on the last occurrence of the separator.
                     Default: _
  --pattern STR      Glob pattern to match FASTQs. Default: *.fastq.gz
  --version          show program's version number and exit
```

###### Examples
Here is an example using the default parameters. In the example, sample *SRR00000* has more than 2 FASTQs matched to it, which is recognized as an error.

```
bactopia prepare tests/dummy-fastqs/
sample  r1      r2
SRR00000        /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000.fastq.gz /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_1.fastq.gz       /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_2.fastq.gz
ERROR: "SRR00000" has more than two different FASTQ files, please check.
```

After tweaking the `--pattern` parameter a little bit. The error is corrected and sample *SRR00000* is properly recognized as a paired-end read set.

```
bactopia prepare tests/dummy-fastqs/ --pattern *_[12].fastq.gz
sample  r1      r2
SRR00000        /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_1.fastq.gz       /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_2.fastq.gz
```

There are a number of ways to tweak the pattern. Just please be sure to give a quick look over of the resulting FOFN.

##### Validating FOFN
When a FOFN is given, the first thing Bactopia does is verify all FASTQ files are found. If everything checks out, each sample will then be processed, otherwise a list of samples with errors will be output to STDERR. 

If you would like to only validate your FOFN (and not run the full pipeline), you can use the `--check_fastqs` parameter.

###### Without Errors
```
bactopia --check_fastqs --fastqs example-data/good-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit3/bactopia/bactopia` [astonishing_colden] - revision: 96c6a1a7ae
Printing what would have been processed. Each line consists of an array of
three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]

Found:
[test001, false, [/home/rpetit3/bactopia/tests/fastqs/test_R1.fastq.gz, /home/rpetit3/bactopia/tests/fastqs/test_R2.fastq.gz]]
[test002, true, [/home/rpetit3/bactopia/tests/fastqs/test.fastq.gz]]
```
Each sample has passed validation and is put into a three element array:

1. sample - the name for this sample
2. is_single_end - the reads are single-end (true) or paired-end (false)
3. fastq_array - the fastqs associated with the sample

This array is then automatically queued up for proccessing by Nextflow.

###### With errors
```
bactopia --check_fastqs --fastqs tests/data/bad-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit3/bactopia/bactopia` [kickass_mestorf] - revision: 222a5ad8b1
LINE 4:ERROR: Please verify /home/rpetit3/bactopia/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit3/bactopia/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit3/bactopia/test/fastqs/test002_R2.fastq.gz exists, and try again
Sample name "test002" is not unique, please revise sample names
The header line (line 1) does not follow expected structure.
Verify sample names are unique and/or FASTQ paths are correct
See "--example_fastqs" for an example
Exiting
```

In the above example, there are mulitple errors. Lines 4 and 5 (`LINE 4:ERROR` or `LINE 5:ERROR`) suggest that based on the given paths the FASTQs do not exist. The sample name `test002` has been used multiple times, and must be corrected. There is also an issue with the header line that must be looked into.

### European Nucleotide Archive
There are a lot of publicly avilable sequences, and you might want to include some of those in your analysis! If that sounds like you, Bactopia has that built in for you! You can give a single *Experiment* accession (`--accession`) or a file where each line is a single *Experiment* accession (`--accessions`). Bactopia will then query ENA to determine *Run* accession(s) associated with the given Experiment accession and proceed download (from ENA) corresponding FASTQ files. After the download is completed, it will be processed through Bactopia.

##### Use --accession for a Single Experiment Accession
`bactopia --accession SRX476958`

##### Use --accessions for Multiple Experiment Accessions
```
bactopia search PRJNA480016 --limit 5
bactopia --accessions ena-accessions.txt \
         --dataset datasets/ \
         --species staphylococcus-aureus \
         --coverage 100 \
         --genome_size median \
         --max_cpus 8 \
         --cpus 2 \
         --outdir ena-multiple-samples
```

## `--max_cpus` & `--cpus`
When Nextflow executes, it uses all available cpus to queue up processes. As you might imagine, if you are on a single server with multiple users, this approach of using all cpus might annoy other users! (Whoops sorry!) To circumvent this feature, two parmeters have been included `--max_cpus` and `--cpus`.

```
    --max_cpus INT          The maximum number of processors this workflow
                                should have access to at any given moment
                                Default: 1

    --cpus INT              Number of processors made available to a single
                                process. If greater than "--max_cpus" it
                                will be set equal to "--max_cpus"
                                Default: 1
```

What `--max_cpus` does is specify to Nextflow the maximum number of cpus it is allowed to occupy at any given time. `--cpus` on the other hand, specifies how many cpus any given step (qc, assembly, annotation, etc...) can occupy. 

By default `--max_cpus` is set to 1 and if `--cpus` is set to a value greater than `--max_cpus` it will be set equal to `--max_cpus`. This appoach errs on the side of caution, by not occupying all cpus on the server without the user's consent!

## `--clean_cache`
Bactopia will keep Nextflow's *work* cache even after successfully completing. While the cache is maintained Bactopia is resumable using the `-resume` parameter. This does however introduce a potentential storage overhead. The cache will contain multiple intermediate files (e.g. uncompressed FASTQs, BAMs, etc...) for each sample that was processed. In other words, it can get pretty large!

If you would like to clean up the cache you can use `--clean_cache`. This will remove the cache **only** after a successful execution (e.g. everything thing finished without errors). This is accomplished by removing the `work` directory created by Nextflow. As you might have guessed, by using removing the cache, Bactopia will no longer be resumeable.

At the end of the day, you can always decide to not use `--clean_cache` and manually remove the `work` directory when you feel it is safe!

## `--keep_all_files`
In some processes, Bactopia will delete large intermediate files (e.g. multiple uncompressed FASTQs) **only** after a process successfully completes. Since this a per-process function, it does not affect Nextflow's ability to resume (`-resume`)a workflow. You can deactivate this feature using `--keep_all_files`. Please, keep in mind the *work* directory is already large, this will make it 2-3 times larger.

# Build Datasets
Bactopia can make use of many existing public datasets, as well as private datasets. The process of downloading, building, and (or) configuring these datasets for Bactopia has been automated for the user.

**Highly recommended to complete this step!**

*This step is completely optional, but it is highly recommended that you do not. By skipping this step of setting up public datasets, Bactopia will be limited to analyses like quality control, assembly, and 31-mer counting.*

## Available Datasets
Some datasets included are applicable to all bacterial species and some are specific to a bacterial species. If specified at runtime, Bactopia will recognize the datasets and execute the appropriate analyses.

### General
**[Ariba's `getref` Reference Datasets](https://github.com/sanger-pathogens/ariba/wiki/Task:-getref)**  
Allows reference datasets (resistance, virulence, and plamids) to be automatically downloaded and configured for usage by Ariba  

**[RefSeq Mash Sketch](https://mash.readthedocs.io/en/latest/data.html)**  
~100,000 genomes and plasmids from NCBI RefSeq, used to give an idea of what is your sequencing data (e.g. Are the sequences what you expected?)  

**[GenBank Sourmash Signatures](https://sourmash.readthedocs.io/en/latest/datasets.html?highlight=--track-abundance#genbank-lca-dataset)**  
~87,000 microbial genomes (includes viral and fungal) from NCBI GenBank, also gives an idea of what is your sequencing data.  

**[PLSDB Mash Sketch & BLAST](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)**  
Includes meta data/annotations, Mash sketches, and BLAST database files of all plasmids stored in PLSDB.  

### Species Specific

**[PubMLST.org MLST Schemas](https://pubmlst.org/datasets/)**  
Multi-locus sequence typing (MLST) allelic profiles and seqeunces for a many different bacterial species (and even a few eukaryotes!).  

**Clustered RefSeq Proteins**  
For the given bacterial species, completed RefSeq genomes are downloaded and then the proteins are clustered and formatted for usage with Prokka.  

**Minmer Sketch of RefSeq Genomes**  
Using the completed genomes downloaded for clustering proteins a Mash sketch and Sourmash signature is created for these genomes. These sketches can then be used for automatic selection of reference genomes for variant calling.

**Optional User Populated Folders**  
A few folders for things such as calling variants, insertion sequences and primers are created that the user can manually populate. 

## Setting Up
Included in Bactopia is the `bactopia datasets` command to automate the process of downloading and/or building these datasets.

### Quick Start
``` bash
bactopia datasets datasets
```

This will set up Ariba datasets (`card` and `vfdb_core`), RefSeq Mash sketch, GenBank Sourmash Signatures, and PLSDB in the newly created `datasets` folder.


### A Single Bacterial Species
``` bash
bactopia datasets datasets --species "Haemophilus influenzae" --include_genus
```

### Multiple Bacterial Species
You can also set up datasets for multiple bacterial species at a time. There are two options to do so.

#### Comma-Separated 
At runtime, you can separate the the different species
``` bash
bactopia datasets datasets --species "Haemophilus influenzae,Staphylococcus aureus" --include_genus
```
#### Text File

In order to do so, you will need to create a text file where each line is the name of a species to set up.

For example, you could create a `species.txt` file and include the following species in it.
``` bash
Haemophilus influenzae
Staphylococcus aureus
Mycobacterium tuberculosis
```

The new command becomes:

``` bash
bactopia datasets datasets --species species.txt --include_genus
```

This will setup the MLST schema (if available) and a protein cluster FASTA file for each species in `species.txt`. 

## Usage
``` 
usage: bactopia datasets [-h] [--ariba STR] [--species STR]
                              [--skip_prokka] [--include_genus]
                              [--identity FLOAT] [--overlap FLOAT]
                              [--max_memory INT] [--fast_cluster]
                              [--skip_minmer] [--skip_plsdb] [--cpus INT]
                              [--clear_cache] [--force] [--force_ariba]
                              [--force_mlst] [--force_prokka]
                              [--force_minmer] [--force_plsdb]
                              [--keep_files] [--list_datasets] [--depends]
                              [--version] [--verbose] [--silent]
                              OUTPUT_DIRECTORY

bactopia datasets - Setup public datasets for Bactopia

positional arguments:
  OUTPUT_DIRECTORY  Directory to write output.

optional arguments:
  -h, --help        show this help message and exit

Ariba Reference Datasets:
  --ariba STR       Setup Ariba datasets for a given reference or a list of
                    references in a text file. (Default: card,vfdb_core)

Bacterial Species:
  --species STR     Download available (cg)MLST schemas and completed genomes
                    for a given species or a list of species in a text file.

Custom Prokka Protein FASTA:
  --skip_prokka     Skip creation of a Prokka formatted fasta for each species
  --include_genus   Include all genus members in the Prokka proteins FASTA
  --identity FLOAT  CD-HIT (-c) sequence identity threshold. (Default: 0.9)
  --overlap FLOAT   CD-HIT (-s) length difference cutoff. (Default: 0.8)
  --max_memory INT  CD-HIT (-M) memory limit (in MB). (Default: unlimited
  --fast_cluster    Use CD-HIT's (-g 0) fast clustering algorithm, instead of
                    the accurate but slow algorithm.

Minmer Datasets:
  --skip_minmer     Skip download of pre-computed minmer datasets (mash,
                    sourmash)

PLSDB (Plasmid) BLAST/Sketch:
  --skip_plsdb      Skip download of pre-computed PLSDB datbases (blast, mash)

Helpful Options:
  --cpus INT        Number of cpus to use. (Default: 1)
  --clear_cache     Remove any existing cache.
  --force           Forcibly overwrite existing datasets.
  --force_ariba     Forcibly overwrite existing Ariba datasets.
  --force_mlst      Forcibly overwrite existing MLST datasets.
  --force_prokka    Forcibly overwrite existing Prokka datasets.
  --force_minmer    Forcibly overwrite existing minmer datasets.
  --force_plsdb     Forcibly overwrite existing PLSDB datasets.
  --keep_files      Keep all downloaded and intermediate files.
  --list_datasets   List Ariba reference datasets and (cg)MLST schemas
                    available for setup.
  --depends         Verify dependencies are installed.

Adjust Verbosity:
  --version         show program's version number and exit
  --verbose         Print debug related text.
  --silent          Only critical errors will be printed.

example usage:
  bactopia datasets outdir
  bactopia datasets outdir --ariba 'card'
  bactopia datasets outdir --species 'Staphylococcus aureus' --include_genus
```

### Useful Parameters
#### --clear_cache
To determine which MLST schemas are available, PubMLST.org is queryied. To prevent a query every run, a list of available schemas is cached to `$HOME/.bactopia/datasets.json`. The cache expires after 15 days, but in case a new species has been made available `--clear_cache` will force a requery of PubMLST.org.

#### --cpus
Increasing `--cpus` (it defaults to 1) is useful for speeding up the download and clustering steps.

#### --force*
If a dataset exists, it will only be overwritten if one of the `--force` parameters are used.
 
#### --include_genus
Completed RefSeq genomes are downloaded for a given species to be used for protein clustering. `--include_genus` will also download completed RefSeq genomes for each genus member.

#### --keep_files
Many intermediate files are downloaded/created (e.g. completed genomes) and deleted during the building process, use `--keep_files` to retain these files.

#### Tweaking CD-HIT
There are parameters (`--identity`, `--overlap`, `--max_memory`, and `--fast_cluster`) to tweak CD-HIT if you find it necessary. Please keep in mind, the only goal of the protein clustering step is to help speed up Prokka, by providing a decent set of proteins to annotate against first.

# Available Public Datasets
Below is a list of public datasets (alphabetical) that could have potentially 
been included during the *[Build Datasets](docs/datasets.md)* step.

### Ariba Reference Datasets
* __[ARG-ANNOT](http://en.mediterranee-infection.com/article.php?laref=283%26titre=arg-annot)__  
_Gupta, S. K. et al. [ARG-ANNOT, a new bioinformatic tool to 
discover antibiotic resistance genes in bacterial genomes.](http://www.ncbi.nlm.nih.gov/pubmed/24145532) 
Antimicrob. Agents Chemother. 58, 212–220 (2014)._  

* __[CARD](https://card.mcmaster.ca/)__  
_McArthur, A. G. et al. [The comprehensive antibiotic resistance database.](http://www.ncbi.nlm.nih.gov/pubmed/23650175) 
Antimicrob. Agents Chemother. 57, 3348–3357 (2013)._  

* __[MEGARes](https://megares.meglab.org/)__  
_Lakin, S. M. et al. [MEGARes: an antimicrobial resistance database for high 
throughput sequencing](http://www.ncbi.nlm.nih.gov/pubmed/27899569). 
Nucleic Acids Res. 45, D574–D580 (2017)._  

* __[NCBI Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047)__  
_Feldgarden, M. et al. [Validating the NCBI AMRFinder Tool and Resistance Gene Database Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of NARMS Isolates](https://doi.org/10.1128/AAC.00483-19). Antimicrob. Agents Chemother. (2019)_ 

* __[plasmidfinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/)__  
_Carattoli, A. et al. [In silico detection and typing of plasmids using 
PlasmidFinder and plasmid multilocus sequence typing.](http://www.ncbi.nlm.nih.gov/pubmed/24777092) 
Antimicrob. Agents Chemother. 58, 3895–3903 (2014)._

* __[resfinder](https://cge.cbs.dtu.dk//services/ResFinder/)__  
_Zankari, E. et al. [Identification of acquired antimicrobial resistance genes.](http://www.ncbi.nlm.nih.gov/pubmed/22782487) 
J. Antimicrob. Chemother. 67, 2640–2644 (2012)._  

* __[SRST2](https://github.com/katholt/srst2)__  
_Inouye, M. et al. [SRST2: Rapid genomic surveillance for public health and 
hospital microbiology labs.](http://www.ncbi.nlm.nih.gov/pubmed/25422674) 
Genome Med. 6, 90 (2014)._  

* __[VFDB](http://www.mgc.ac.cn/VFs/)__  
_Chen, L., Zheng, D., Liu, B., Yang, J. & Jin, Q. [VFDB 2016: hierarchical 
and refined dataset for big data analysis--10 years on.](http://www.ncbi.nlm.nih.gov/pubmed/26578559) 
Nucleic Acids Res. 44, D694–7 (2016)._  

* __[VirulenceFinder](https://cge.cbs.dtu.dk/services/VirulenceFinder/)__  
_Joensen, K. G. et al. [Real-time whole-genome sequencing for routine typing, 
surveillance, and outbreak detection of verotoxigenic Escherichia coli.](http://www.ncbi.nlm.nih.gov/pubmed/24574290) 
J. Clin. Microbiol. 52, 1501–1510 (2014)._  

### Minmer Datasets
* __[Mash Refseq (release 88) Sketch](https://mash.readthedocs.io/en/latest/data.html)__  
_Ondov, B. D. et al. [Mash Screen: High-throughput sequence containment 
estimation for genome discovery.](https://doi.org/10.1101/557314) bioRxiv 557314 (2019)._  

* __[Sourmash Genbank LCA Signature](https://sourmash.readthedocs.io/en/latest/databases.html)__  
_Titus Brown, C. & Irber, L. [sourmash: a library for MinHash sketching of DNA.](http://joss.theoj.org/papers/10.21105/joss.00027) 
JOSS 1, 27 (2016)._  

### Everything Else
* __[NCBI RefSeq Database](https://www.ncbi.nlm.nih.gov/refseq/)__  
_O’Leary, N. A. et al. [Reference sequence (RefSeq) database at NCBI: current status, 
taxonomic expansion, and functional annotation](http://dx.doi.org/10.1093/nar/gkv1189). 
Nucleic Acids Res. 44, D733–45 (2016)._  

* __[PLSDB - A plasmid database](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)__  
_Galata, V., Fehlmann, T., Backes, C. & Keller, A. [PLSDB: a resource of complete bacterial 
plasmids](http://dx.doi.org/10.1093/nar/gky1050). 
Nucleic Acids Res. 47, D195–D202 (2019)._  

* __[PubMLST.org](https://pubmlst.org/)__  
_Jolley, K. A., Bray, J. E. & Maiden, M. C. J. [Open-access bacterial population genomics: BIGSdb 
software, the PubMLST.org website and their applications](http://dx.doi.org/10.12688/wellcomeopenres.14826.1). 
Wellcome Open Res 3, 124 (2018)._  

# Software Included In Bactopia
Below is a list of software (alphabetical) used (directly and indirectly) by Bactopia. A link to the software page as well as the citation (if available) have been included.

* __[AMRFinderPlus](https://github.com/ncbi/amr)__  
Find acquired antimicrobial resistance genes and some point mutations in protein or assembled nucleotide sequences.  
_Feldgarden, M. et al. [Validating the NCBI AMRFinder Tool and Resistance Gene Database Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of NARMS Isolates](https://doi.org/10.1128/AAC.00483-19). Antimicrob. Agents Chemother. (2019)_  

* __[Aragorn](http://130.235.244.92/ARAGORN/Downloads/)__  
Finds transfer RNA features (tRNA)  
_Laslett D. and B. Canback, [ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences.](https://doi.org/10.1093/nar/gkh152) Nucleic Acids Res. 32(1):11-6. (2004)_  

* __[Ariba](https://github.com/sanger-pathogens/ariba)__  
Antimicrobial Resistance Identification By Assembly  
_Hunt, M. et al. [ARIBA: rapid antimicrobial resistance genotyping directly from 
sequencing reads](http://dx.doi.org/10.1099/mgen.0.000131). 
Microb Genom 3, e000131 (2017)._  

* __[Assembly-Scan](https://github.com/rpetit3/assembly-scan)__  
Generate basic stats for an assembly.  
_Petit III, R. A. [assembly-scan: generate basic stats for an 
assembly](https://github.com/rpetit3/assembly-scan)._  

* __[Barrnap](https://github.com/tseemann/barrnap)__  
Bacterial ribosomal RNA predictor  
_Seemann, T. [Barrnap: Bacterial ribosomal RNA predictor](https://github.com/tseemann/barrnap)_  

* __[BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)__  
BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data.  
_Bushnell, B. [BBMap short read aligner, and other bioinformatic tools.](http://sourceforge.net/projects/bbmap/)_  

* __[BCFtools](https://github.com/samtools/bcftools)__  
Utilities for variant calling and manipulating VCFs and BCFs.  
_Danecek, P. et al. [BCFtools - Utilities for variant calling and manipulating VCFs and BCFs.](http://github.com/samtools/bcftools)_  

* __[Bedtools](https://github.com/arq5x/bedtools2)__  
A powerful toolset for genome arithmetic.  
_Quinlan, A. R. & Hall, I. M. [BEDTools: a flexible suite of utilities for 
comparing genomic features](http://dx.doi.org/10.1093/bioinformatics/btq033). 
Bioinformatics 26, 841–842 (2010)._  

* __[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)__  
Basic Local Alignment Search Tool  
_Camacho, C. et al. [BLAST+: architecture and applications](http://dx.doi.org/10.1186/1471-2105-10-421). 
BMC Bioinformatics 10, 421 (2009)._  

* __[BWA](https://github.com/lh3/bwa/)__  
Burrow-Wheeler Aligner for short-read alignment  
_Li, H. [Aligning sequence reads, clone sequences and assembly contigs with 
BWA-MEM](http://arxiv.org/abs/1303.3997). arXiv [q-bio.GN] (2013)._  

* __[CD-Hit](https://github.com/weizhongli/cdhit)__  
Accelerated for clustering the next-generation sequencing data  
_Li, W. & Godzik, A. [Cd-hit: a fast program for clustering and comparing large sets of protein 
or nucleotide sequences](http://dx.doi.org/10.1093/bioinformatics/btl158). 
Bioinformatics 22, 1658–1659 (2006)._  
_Fu, L., Niu, B., Zhu, Z., Wu, S. & Li, W. [CD-HIT: accelerated for clustering the next-generation 
sequencing data](http://dx.doi.org/10.1093/bioinformatics/bts565). 
Bioinformatics 28, 3150–3152 (2012)._  

* __[FastQC](https://github.com/s-andrews/FastQC)__  
A quality control analysis tool for high throughput sequencing data.  
_Andrews, S. [FastQC: a quality control tool for high throughput sequence data.](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)._  

* __[Fastq-Scan](https://github.com/rpetit3/fastq-scan)__  
Output FASTQ summary statistics in JSON format  
_Petit III, R. A. [fastq-scan: generate summary statistics of input FASTQ sequences.](https://github.com/rpetit3/fastq-scan)_  

* __[FLASH](https://ccb.jhu.edu/software/FLASH/)__  
A fast and accurate tool to merge paired-end reads.  
_Magoč, T., and S. L. Salzberg, [FLASH: fast length adjustment of short reads to improve genome assemblies.](https://doi.org/10.1093/bioinformatics/btr507) Bioinformatics 27.21 (2011): 2957-2963._  

* __[freebayes](https://github.com/ekg/freebayes)__  
Bayesian haplotype-based genetic polymorphism discovery and genotyping  
_Garrison E., and G. Marth, [Haplotype-based variant detection from short-read sequencing.](https://arxiv.org/abs/1207.3907) arXiv preprint arXiv:1207.3907 [q-bio.GN] (2012)_  

* __[GNU Parallel](https://www.gnu.org/software/parallel/)__  
A shell tool for executing jobs in parallel  
_Tange, O. [GNU Parallel](https://doi.org/10.5281/zenodo.1146014) 2018, March 2018_  

* __[HMMER](http://hmmer.org/)__  
Biosequence analysis using profile hidden Markov models  
_Finn R. D. et al. [HMMER web server: interactive sequence similarity searching.](https://doi.org/10.1093/nar/gkr367) Nucleic Acids Res. ;39:W29-37. (2011)_  

* __[Infernal](http://eddylab.org/infernal/)__  
Searches DNA sequence databases for RNA structure and sequence similarities  
_Nawrocki, E. P., and S. R. Eddy, [Infernal 1.1: 100-fold faster RNA homology searches.](https://doi.org/10.1093/bioinformatics/btt509) Bioinformatics, 29(22), 2933-2935. (2013)_  

* __[ISMapper](https://github.com/jhawkey/IS_mapper)__  
IS mapping software  
_Hawkey, J. et al. [ISMapper: identifying transposase insertion sites in 
bacterial genomes from short read sequence data](http://dx.doi.org/10.1186/s12864-015-1860-2). 
BMC Genomics 16, 667 (2015)._  

* __[Lighter](https://github.com/mourisl/Lighter)__  
Fast and memory-efficient sequencing error corrector  
_Song, L., Florea, L. and B. Langmead, [Lighter: Fast and Memory-efficient Sequencing Error Correction without Counting](http://genomebiology.com/2014/15/11/509/). Genome Biol. 2014 Nov 15;15(11):509._

* __[Mash](https://github.com/marbl/Mash)__  
Fast genome and metagenome distance estimation using MinHash  
_Ondov, B. D. et al. [Mash: fast genome and metagenome distance 
estimation using MinHash](http://dx.doi.org/10.1186/s13059-016-0997-x). 
Genome Biol. 17, 132 (2016)._  
_Ondov, B. D. et al. [Mash Screen: High-throughput sequence 
containment estimation for genome discovery](http://dx.doi.org/10.1101/557314). 
bioRxiv 557314 (2019)._  

* __[McCortex](https://github.com/mcveanlab/mccortex)__  
De novo genome assembly and multisample variant calling  
_Turner, I., Garimella, K. V., Iqbal, Z. and G. McVean, [Integrating long-range 
connectivity information into de Bruijn graphs.](http://dx.doi.org/10.1093/bioinformatics/bty157) 
Bioinformatics 34, 2556–2565 (2018)._  

* __[MEGAHIT](https://github.com/voutcn/megahit)__  
Ultra-fast and memory-efficient (meta-)genome assembler  
_Li, D., et al. [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.](https://doi.org/10.1093/bioinformatics/btv033) Bioinformatics 31.10 (2015): 1674-1676._  

* __[MinCED](https://github.com/ctSkennerton/minced)__  
Mining CRISPRs in Environmental Datasets  
_Skennerton, C. [MinCED: Mining CRISPRs in Environmental Datasets](https://github.com/ctSkennerton/minced)_  

* __[Minimap2](https://github.com/lh3/minimap2)__
A versatile pairwise aligner for genomic and spliced nucleotide sequences  
_Li, H. [Minimap2: pairwise alignment for nucleotide sequences.](https://doi.org/10.1093/bioinformatics/bty191) Bioinformatics, 34:3094-3100. (2018)_  

* __[NCBI Genome Download](https://github.com/kblin/ncbi-genome-download)__  
Scripts to download genomes from the NCBI FTP servers  
_Blin, K. [NCBI Genome Download: Scripts to download genomes from the NCBI FTP 
servers](https://github.com/kblin/ncbi-genome-download)_  

* __[Nextflow](https://github.com/nextflow-io/nextflow)__  
A DSL for data-driven computational pipelines.  
_Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P.P., Palumbo, E., Notredame, C., 2017. [Nextflow enables reproducible computational workflows.](https://www.nature.com/articles/nbt.3820.pdf?origin=ppub) Nat. Biotechnol. 35, 316–319._

* __[Pigz](https://zlib.net/pigz/)__  
A parallel implementation of gzip for modern multi-processor, multi-core machines.  
_Adler, M. [pigz: A parallel implementation of gzip for modern multi-processor, multi-core machines.](https://zlib.net/pigz/) Jet Propulsion Laboratory (2015)._  

* __[Pilon](https://github.com/broadinstitute/pilon/)__  
An automated genome assembly improvement and variant detection tool  
_Walker, B. J., et al. [Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement.](https://doi.org/10.1371/journal.pone.0112963) PloS one 9.11 (2014): e112963._  

* __[Prodigal](https://github.com/hyattpd/Prodigal)__  
Fast, reliable protein-coding gene prediction for prokaryotic genomes.  
_Hyatt, D., et al. [Prodigal: prokaryotic gene recognition and translation initiation site identification.](https://doi.org/10.1186/1471-2105-11-119) BMC Bioinformatics 11.1 (2010): 119._  

* __[Prokka](https://github.com/tseemann/prokka)__  
Rapid prokaryotic genome annotation  
_Seemann, T. [Prokka: rapid prokaryotic genome annotation](http://dx.doi.org/10.1093/bioinformatics/btu153). 
Bioinformatics 30, 2068–2069 (2014)._  

* __[RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/)__  
Consistent and rapid annotation of ribosomal RNA genes  
_Lagesen, K., et al. [RNAmmer: consistent annotation of rRNA genes in genomic sequences.](https://dx.doi.org/10.1093%2Fnar%2Fgkm160) Nucleic Acids Res 35.9: 3100-3108. (2007)_  

* __[samclip](https://github.com/tseemann/samclip)__  
Filter SAM file for soft and hard clipped alignments  
_Seemann, T. [Samclip: Filter SAM file for soft and hard clipped alignments](https://github.com/tseemann/samclip)_  

* __[Samtools](https://github.com/samtools/samtools)__  
Tools for manipulating next-generation sequencing data  
_Li, H. et al. [The Sequence Alignment/Map format and SAMtools](http://dx.doi.org/10.1093/bioinformatics/btp352). 
Bioinformatics 25, 2078–2079 (2009)._

* __[Seqtk](https://github.com/lh3/seqtk)__  
A fast and lightweight tool for processing sequences in the FASTA or FASTQ format.  
_Li, H. [Toolkit for processing sequences in FASTA/Q formats](https://github.com/lh3/seqtk)_  

* __[Shovill](https://github.com/tseemann/shovill)__  
Faster assembly of Illumina reads  
_Seemann, T. [Shovill: De novo assembly pipeline for Illumina paired reads](https://github.com/tseemann/shovill)_  

* __[SignalP](http://www.cbs.dtu.dk/services/SignalP-4.0/)__  
Finds signal peptide features in CDS  
_Petersen, T. N., et al. [SignalP 4.0: discriminating signal peptides from transmembrane regions.](https://doi.org/10.1038/nmeth.1701) Nature methods 8.10: 785.(2011)_  

* __[SKESA](https://github.com/ncbi/SKESA)__  
Strategic Kmer Extension for Scrupulous Assemblies  
_Souvorov, A., Agarwala, R. and D. J. Lipman. [SKESA: strategic k-mer extension for scrupulous assemblies.](https://doi.org/10.1186/s13059-018-1540-z) Genome Biology 19:153 (2018).__  

* __[Snippy](https://github.com/tseemann/snippy)__  
Rapid haploid variant calling and core genome alignment  
_Seemann, T. [Snippy: fast bacterial variant calling from NGS reads](https://github.com/tseemann/snippy)
(2015)_  

* __[SnpEff](http://snpeff.sourceforge.net/)__  
Genomic variant annotations and functional effect prediction toolbox.  
_Cingolani, P., et al. [A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.](https://doi.org/10.4161/fly.19695) Fly, 6(2), 80-92 (2012)_  

* __[SNP-sites](https://github.com/sanger-pathogens/snp-sites)__  
Rapidly extracts SNPs from a multi-FASTA alignment.  
_Page, A. J., et al. [SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments.](https://dx.doi.org/10.1099%2Fmgen.0.000056) Microbial Genomics 2.4 (2016)._  

* __[Sourmash](https://github.com/dib-lab/sourmash)__  
Compute and compare MinHash signatures for DNA data sets.  
_Titus Brown, C. and L. Irber [sourmash: a library for MinHash sketching 
of DNA](http://dx.doi.org/10.21105/joss.00027). JOSS 1, 27 (2016)._  

* __[SPAdes](https://github.com/ablab/spades)__  
An assembly toolkit containing various assembly pipelines.  
_Bankevich, A., et al. [SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing.](https://doi.org/10.1089/cmb.2012.0021) Journal of computational biology 19.5 (2012): 455-477._  

* __[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)__  
A flexible read trimming tool for Illumina NGS data  
_Bolger, A. M., Lohse, M., and B. Usadel. [Trimmomatic: a flexible trimmer for Illumina sequence data.](https://doi.org/10.1093/bioinformatics/btu170) Bioinformatics 30.15 (2014): 2114-2120._  

* __[VCF-Annotator](https://github.com/rpetit3/vcf-annotator)__  
Add biological annotations to variants in a VCF file.  
_Petit III, R. A. [VCF-Annotator: Add biological annotations to variants 
in a VCF file.](https://github.com/rpetit3/vcf-annotator)._  

* __[Vcflib]()__  
a simple C++ library for parsing and manipulating VCF files  
_Garrison, E. [Vcflib: A C++ library for parsing and manipulating VCF files](https://github.com/vcflib/vcflib)_  

* __[Velvet](https://github.com/dzerbino/velvet)__  
Short read de novo assembler using de Bruijn graphs  
_Zerbino, D. R., and E. Birney. [Velvet: algorithms for de novo short read assembly using de Bruijn graphs.](http://www.genome.org/cgi/doi/10.1101/gr.074492.107) Genome research 18.5 (2008): 821-829._  

* __[vt](https://github.com/atks/vt)__  
A tool set for short variant discovery in genetic sequence data.  
_Tan, A., Abecasis, G. R., and H. M. Kang, [Unified representation of genetic variants.](https://doi.org/10.1093/bioinformatics/btv112) Bioinformatics, 31(13), 2202-2204. (2015)_  

## Please Cite Datasets and Tools
If you have used Bactopia in your work, please be sure to cite any datasets or tools you may have used. A citation link for each dataset/tool has been made available. If a citation needs to updated please let me know!

A BibTeX file of each citation is also available at [Bactopia Datasets and Software BibTeX](docs/bactopia-datasets-software.bib)

# Acknowledgements
Bactopia is truly a case of *"standing upon the shoulders of giants"*. As seen above, nearly every component of Bactopia was created by others and made freely available to the public.

I would like to personally extend my many thanks and gratitude to the authors
of these software packages and public datasets. If you've made it this far, I 
owe you a beer 🍻 (or coffee ☕!) if we ever encounter one another in person. 
Really, thank you very much!

