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
A script named `prepare-fofn` has been included to help aid (hopefully!) the process of creating a FOFN for your samples. This script will attempt to find FASTQ files in a given directory and output the expected FOFN format. It will also output any potential issues associated with the pattern matching.

!!! error "Verify accuracy of FOFN"
    This is currently an experimental function. There are likely bugs to be ironed out. Please be sure to give the resulting FOFN a quick look over.

###### Usage
```
prepare-fofn [-h] [-e STR] [-s STR] [--pattern STR] [--version] STR

prepare-fofn - Read a directory and prepare a FOFN of FASTQs

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
prepare-fofn  tests/dummy-fastqs/
sample  r1      r2
SRR00000        /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000.fastq.gz /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_1.fastq.gz       /home/rpetit/projects/bactopia/bactopia/tests/dummy-fastqs/SRR00000_2.fastq.gz
ERROR: "SRR00000" has more than two different FASTQ files, please check.
```

After tweaking the `--pattern` parameter a little bit. The error is corrected and sample *SRR00000* is properly recognized as a paired-end read set.

```
prepare-fofn  tests/dummy-fastqs/ --pattern *_[12].fastq.gz
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

!!! info "Use --accession for a Single Experiment Accession"
    `bactopia --accession SRX476958`

!!! info "Use --accessions for Multiple Experiment Accessions"
    `bactopia --accessions my-accessions.txt`

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
