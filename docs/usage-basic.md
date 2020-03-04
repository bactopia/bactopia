# Basic Usage For Bactopia
Bactopia is a wrapper around many different tools. Each of these tools may (or may not) have there own configurable parameters for you to tweak. In order to facilitate getting started with Bactopia, this section has been limited to discussion of only a few parameters. However, if you are interested in the full list of configurable parameters in Bactopia, please check out the [Complete Usage](usage-complete.md) section.

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

    --genome_size INT       Expected genome size (bp) for all samples, a value of '0'
                                will disable read error correction and read subsampling.
                                Special values (requires --species):
                                    'min': uses minimum completed genome size of species
                                    'median': uses median completed genome size of species
                                    'mean': uses mean completed genome size of species
                                    'max': uses max completed genome size of species
                                Default: Mash estimate

    --outdir DIR            Directory to write results to
                                Default: .

Nextflow Queue Parameters:
    At execution, Nextflow creates a queue and the number of slots in the queue is determined by the total number 
    of cores on the system. When a task is submitted to the queue, the total number of slots it occupies is 
    determined by the value set by "--cpus". 

    This can have a significant effect on the efficiency of the Nextflow's queue system. If "--cpus" is set to a 
    value that is equal to the number of cores availabe, in most cases only a single task will be able to run 
    because its occupying all available slots.

    When in doubt, "--cpus 4" is a safe bet, it is also the default value if you don't use "--cpus".

    --max_time INT          The maximum number of minutes a task should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single task.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single task.
                                Default: 4

Nextflow Related Parameters:
    --infodir DIR           Directory to write Nextflow summary files to
                                Default: ./bactopia-info

    --condadir DIR          Directory to Nextflow should use for Conda environments
                                Default: Bactopia's Nextflow directory

    --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows 
                                you to create profiles specific to your environment (e.g. SGE,
                                AWS, SLURM, etc...). This config file is loaded last and will 
                                overwrite existing variables if set.
                                Default: Bactopia's default configs

    --nfdir                 Print directory Nextflow has pulled Bactopia to

    --overwrite             Nextflow will overwrite existing output files.
                                Default: false

    --conatainerPath        Path to Singularity containers to be used by the 'slurm'
                                profile.
                                Default: /opt/bactopia/singularity

    --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                will wait before execution.
                                Default: 5 seconds

    -resume                 Nextflow will attempt to resume a previous run. Please notice it is 
                                only a single '-'

Useful Parameters:
    --available_datasets    Print a list of available datasets found based
                                on location given by "--datasets"

    --example_fastqs        Print example of expected input for FASTQs file

    --check_fastqs          Verify "--fastqs" produces the expected inputs

    --compress              Compress (gzip) select outputs (e.g. annotation, variant calls)
                                to reduce overall storage footprint.

    --keep_all_files        Keeps all analysis files created. By default, intermediate
                                files are removed. This will not affect the ability
                                to resume Nextflow runs, and only occurs at the end
                                of the process.


    --dry_run               Mimics workflow execution, to help determine if conda environments
                                or container images are properly set up.

    --version               Print workflow version information

    --help                  Show this message and exit

    --help_all              Show a complete list of adjustable parameters
```

## FASTQ Inputs
Bactopia has multiple approaches to specify your input sequences. You can make use of your local FASTQs or download FASTQs from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena). Which approach really depends on what you need to achieve! The following sections describe methods to process single samples, multiple samples, downloading samples from the ENA.

### Local
#### Single Sample
When you only need to process a single sample at a time, Bactopia allows that! You only have to the sample name (`--sample`) and the whether the read set is paired-end (`--R1` and `--R2`) or a single-end (`--SE`). 

!!! info "Use --R1, --R2 for Paired-End FASTQs"
    `bactopia --sample my-sample --R1 /path/to/my-sample_R1.fastq.gz --R2 /path/to/my-sample_R2.fastq.gz`

!!! info "Use --SE for Single-End FASTQs"
    `bactopia --sample my-sample --SE /path/to/my-sample.fastq.gz`

#### Multiple Samples
For multiple samples, you must create a file with information about the inputs, a *file of filenames* (FOFN). This file specifies sample names and location of FASTQs to be processed. Using this information, paired-end or single-end information can be extracted as well as naming output files.

While this is an additional step for you, the user, it helps to avoid potential pattern matching errors. 

Most importantly, by taking this approach, you can process hundreds of samples in a single command. There is also the added benefit of knowing which FASTQs were analysed and their location at a later time!

!!! info "Use --fastqs for Multiple Samples"
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
`bactopia prepare` has been included to help aid (hopefully!) the process of creating a FOFN for your samples. This script will attempt to find FASTQ files in a given directory and output the expected FOFN format. It will also output any potential issues associated with the pattern matching.

!!! error "Verify accuracy of FOFN"
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
bactopia prepare  tests/dummy-fastqs/
sample  r1      r2
SRR00000        /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000.fastq.gz /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_1.fastq.gz       /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_2.fastq.gz
ERROR: "SRR00000" has more than two different FASTQ files, please check.
```

After tweaking the `--pattern` parameter a little bit. The error is corrected and sample *SRR00000* is properly recognized as a paired-end read set.

```
bactopia prepare  tests/dummy-fastqs/ --pattern *_[12].fastq.gz
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

### ENA & SRA
There are a lot of publicly avilable sequences available from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) (ENA) and the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA). There's a good chance you might want to include some of those sequences in your analysis! If that sounds like you, Bactopia has that built in for you! You can give a single *Experiment* accession (`--accession`) or a file where each line is a single *Experiment* accession (`--accessions`). Bactopia will then query ENA to determine *Run* accession(s) associated with the given Experiment accession and proceed download the corresponding FASTQ files from either the SRA (default) or ENA (`--use_ena`). After the download is completed, it will be processed through Bactopia.

!!! info "Use --accession for a Single Experiment Accession"
    SRA: `bactopia --accession SRX476958`
    ENA: `bactopia --accession SRX476958 --use_ena`

!!! info "Use --accessions for Multiple Experiment Accessions"
    SRA: `bactopia --accessions my-accessions.txt`
    ENA: `bactopia --accessions my-accessions.txt --use_ena`

#### Generating Accession List
`bactopia search` has been made to help assist in generating a list of Experiment accessions to be procesed by Bactopia (via `--accessions`). Users can provide a Taxon ID (e.g. 1280), a binary name (e.g. Staphylococcus aureus), or Study accessions (e.g. PRJNA480016). This value is then queried against ENA's [Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)), and a list of all Experiment accessions associated with the query is returned.

##### Usage
```
usage: bactopia search [-h] [--exact_taxon] [--outdir OUTPUT_DIRECTORY]
                       [--prefix PREFIX] [--limit INT] [--version]
                       STR

bactopia search - Search ENA for associated WGS samples

positional arguments:
  STR                   Taxon ID or Study accession

optional arguments:
  -h, --help            show this help message and exit
  --exact_taxon         Exclude Taxon ID descendents.
  --outdir OUTPUT_DIRECTORY
                        Directory to write output. (Default: .)
  --prefix PREFIX       Prefix to use for output file names. (Default: ena)
  --limit INT           Maximum number of results to return. (Default:
                        1000000)
  --version             show program's version number and exit

example usage:
  bactopia search PRJNA480016 --limit 20
  bactopia search 1280 --exact_taxon --limit 20'
  bactopia search "staphylococcus aureus" --limit 20
```

##### Example
```
bactopia search PRJNA480016 --limit 5
```

When completed three files are produced:

1. `ena-accessions.txt` - Contains a list of Experiment accessions to be processed.
   ```
   SRX4563686
   SRX4563689
   SRX4563687
   SRX4563690
   SRX4563688
   ```

!!! info "Input for Bactopia"
    This file can be used in conjunction with the `--accessions` parameter for Bactopia processing.


2. `ena-results.txt` - Contains the full results of the API query. This includes multiples fields (sample_accession, tax_id, sample_alias, center_name, etc...)

3. `ena-summary.txt` - Contains a small summary of the completed request
    ```
    QUERY: (study_accession=PRJNA480016 OR secondary_study_accession=PRJNA480016)
    LIMIT: 5
    RESULTS: 5 (./ena-results.txt)
    ILLUMINA ACCESSIONS: 5 (./ena-accessions.txt)
    ```


## `--cpus`
At execution, Nextflow creates a queue and the number of slots in the queue is determined by the total number of cores on the system. So if you have a 24-core system, that means Nextflow will have a queue with 24-slots available. This feature kind of makes `--cpus` a little misleading. Typically when you give `--cpus` you are saying *"use this amount of cpus"*. But that is not the case for Nextflow and Bactopia. When you use `--cpus` what you are actually saying is *"for any particular task, use this amount of slots"*. Commands within a task processors will use the amount specified by `--cpus`.

!!! error "`--cpus` can have a significant effect on the efficiency of Bactopia"
    So for example if you have a system with 24-cores.

    This command, `bactopia ... --cpus 24`, says *for any particular task, use 24 slots*. Nextflow will give tasks in Bactopia 24 slots out of 24 available (24-core machine). In other words the queue can one have one task running at once because each task occupies 24 slots.

    On the other hand, `bactopia ... --cpus 4` says *for any particular task, use 4 slots*. Now, for Nextflow will give each task 4 slots out of 24 slots. Which means 6 tasks can be running at once. This can lead to much better efficiency because less jobs are stuck waiting in line. 

    There are some tasks in Bactopia that will only ever use a single slot because they are single-core tasks. But for example the `annotation` step will always use the number of slots specified by `--cpus`. If the `--cpus` is too high, the `annotation` will get bogged down, which causes tasks dependent on `annotation` to also get bogged down.

!!! info "When in doubt `--cpus 4` is a safe value."
    This is also the default value for Bactopia.


## `--genome_size`
Throughout the Bactopia workflow a genome size is used for various tasks. By default, a genome size is estimated using Mash. However, users can provide their own value for genome size, use values based on [Species Specific Datasets](/datasets/#species-specific), or completely disable it.

| Value | Result |
|-------|--------|
| *empty*  | Mash is used to estimate the genome size |
| integer | Uses the genome size (e.g. `--genome_size 2800000`) provided by the user |
| 0 | Read error correct and read subsampling will be disabled. |
| min | Requires `--species`, the minimum completed genome size for a species is used |
| median | Requires `--species`, the median completed genome size for a species is used | 
| mean |  Requires `--species`, the mean completed genome size for a species is used | 
| max | Requires `--species`, the maximum completed genome size for a species is used | 

!!! error "Mash may not be the most accurate estimate"
    Mash is very convenient to quickly estimate a genome size, but it may not be the most accurate in all cases and will differ between samples. It is recommended that when possible a known genome size or one based off completed genomes should be used. 

## `--nfconfig`
A key feature of Nextflow is you can provide your own config files. What this boils down to you can easily set Bactopia to run on your environment. With `--nfconfig` you can tell Bactopia to import your config file. 

`--nfconfig` has been set up so that it is the last config file to be loaded by Nextflow. This means that if your config file contains variables (e.g. params or profiles) already set they will be overwritten by your values.

[Nextflow goes into great details on how to create configuration files.](https://www.nextflow.io/docs/latest/config.html) Please check the following links for adjustsments you be interested in making.

| Scope   | Description |
|---------|-------------|
| [env](https://www.nextflow.io/docs/latest/config.html#scope-env)     | Set any environment variables that might be required |
| [params](https://www.nextflow.io/docs/latest/config.html#scope-params)  | Change the default values of command line arguments  |
| [process](https://www.nextflow.io/docs/latest/config.html#scope-process) | Adjust perprocess configurations such as containers, conda envs, or resource usage |
| [profile](https://www.nextflow.io/docs/latest/config.html#config-profiles) | Create predefined profiles for your [Executor](https://www.nextflow.io/docs/latest/operator.html#filtering-operators) |

There are [many other scopes](https://www.nextflow.io/docs/latest/config.html#config-scopes) that you might be interested in checking out.

You are most like going to want to create a custom profile. By doing so you can specify it at runtime (`-profile myProfile`) and Nextflow will be excuted based on that profile. Often times your custom profile will include information on the executor (queues, allocations, apths, etc...).

If you need help please [reach out](https://github.com/bactopia/bactopia/issues/new/choose)!

*If you're using the standard profile (did not specify -profile 'xyz') this might not be necessary.*

## `-resume`
Bactopia relies on [Nextflow's Resume Feature](https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume) to resume runs. You can tell Bactopia to resume by adding `-resume` to your command line. When `-resume` is used, Nextflow will review the cache and determine if the previous run is resumable. If the previous run is not resumable, execution will
start at the beginning.

## `--keep_all_files`
In some processes, Bactopia will delete large intermediate files (e.g. multiple uncompressed FASTQs) **only** after a process successfully completes. Since this a per-process function, it does not affect Nextflow's ability to resume (`-resume`)a workflow. You can deactivate this feature using `--keep_all_files`. Please, keep in mind the *work* directory is already large, this will make it 2-3 times larger.




