# Basic Usage For Bactopia
Bactopia is a wrapper around many different tools. Each of these tools may (or may not) have there own configurable parameters for you to tweak. In order to facilitate getting started with Bactopia, this section has been limited to discussion of only a few parameters. However, if you are interested in the full list of configurable parameters in Bactopia, please check out the [Complete Usage](usage-complete.md).

## Usage
```
bactopia v0.0.1

Required Parameters:
    --fastqs STR            An input file containing the sample name and
                                absolute paths to FASTQs to process

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

    --keep_cache            Keeps 'work' and '.nextflow' logs.
                                Default: Delete on successful completion

    --version               Print workflow version information
    --help                  Show this message and exit
    --help_all              Show a complete list of adjustable parameters
```

## `--fastqs`
Bactopia was developed to take in a file with information about the inputs, a *file of filenames* (FOFN). Essentially what this file does is specify the sample's name and location of FASTQs to be processed. With this information, paired-end/single-end information can be extracted as well as naming output files.

While this is an additional step for you, the user, it avoids potenial bugs associated with pattern matching FASTQs to extract a name and determine the SE vs PE status. 

Most importantly, by taking this approach, Nextflow can easily queue the analysis of a single FASTQ or hundreds of FASTQs in a single command. There is also the added benefit of knowing which FASTQs were analysed and their location at a later time!

### Expected Format
You can use the `--example_fastqs` to get an example of the expected structure for the input FASTQs FOFN.

```
illumina-cleanup --example_fastqs
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

All lines after the header line, contain unique sample names and location(s) to associated FASTQ file(s). For the purpose of the FOFN, absolute paths are preferred as they are more informative in the future. Relative paths may be ok, but have not been thoroughly tested. 

In the example above, two samples would be cleaned up. Sample `test001` has two FASTQs and would be processed as pair-end reads. While sample `test002` only has a single FASTQ and would be processed as single-end reads.

### Validating
When a FOFN is given, the first thing *illumina-cleanup* does is verify the path all FASTQ files is valid. If all paths are valid, each sample will then be processed, otherwise a list of samples with errors will be output to STDERR. 

If you would like to only validate your FOFN (and not run the full pipeline), you can use the `--check_fastqs` parameter.

#### Without Errors
```
illumina-cleanup --check_fastqs --fastqs example-data/good-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [astonishing_colden] - revision: 96c6a1a7ae
Printing what would have been processed. Each line consists of an array of
three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]

Found:
[test001, false, [/home/rpetit/illumina-cleanup/test/fastqs/test_R1.fastq.gz, /home/rpetit/illumina-cleanup/test/fastqs/test_R2.fastq.gz]]
[test002, true, [/home/rpetit/illumina-cleanup/test/fastqs/test.fastq.gz]]
```
Each sample has passed validation and is put into a three element array:

1. sample - the name for this sample
2. is_single_end - the reads are single-end (true) or paired-end (false)
3. fastq_array - the fastqs associated with the sample

This array is then automatically queued up for proccessing by Nextflow.

#### With errors
```
illumina-cleanup --check_fastqs --fastqs example-data/bad-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [kickass_mestorf] - revision: 222a5ad8b1
LINE 4:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test002_R2.fastq.gz exists, and try again
Sample name "test002" is not unique, please revise sample names
The header line (line 1) does not follow expected structure.
Verify sample names are unique and/or FASTQ paths are correct
See "--example_fastqs" for an example
Exiting
```

In the above example, there are mulitple errors. Lines 4 and 5 (`LINE 4:ERROR` or `LINE 5:ERROR`) suggest that based on the given paths the FASTQs do not exist. The sample name `test002` has been used multiple times, and must be corrected. There is also an issue with the header line that must be looked into.


## `--max_cpus` & `--cpus`
By default when Nextflow executes, it uses all available cpus to queue up processes. As you might imagine, if you are on a single server with multiple users, this approach of using all cpus might annoy other users! (Whoops sorry!) To circumvent this feature, two parmeters have been included `--max_cpus` and `--cpus`.

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


## `--keep_cache`
By default, *Bactopia* will remove all intermediate files **after successfully completing all steps**. It does this by removing the `work` and `.nextflow` directories created by Nextflow. An important consequence of this action is you can no longer resume the workflow. But, again, clearing the cache only occurs after all samples have been successfully processed. If there is an error during execution, the cache will remain intact and the process can still be resumed (`-resume`). 

If you would like to retain the cache files, you can use `--keep_cache` to do so. This is useful to use when you are still working out which data you want to include in your datasets. For example, if you are still deciding on which reference genome(s) to call variants against.

Although, please keep in mind there is storage overhead for maintaining the cache. It will contain multiple intermediate files (e.g. uncompressed FASTQs, BAMs, etc...) for each sample that was processed. In other words, it can get pretty large!

