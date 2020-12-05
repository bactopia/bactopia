# Basic Usage For Bactopia
Bactopia is a wrapper around many different tools. Each of these tools may (or may not) have there own configurable parameters for you to tweak. In order to facilitate getting started with Bactopia, this section has been limited to discussion of only a few parameters. However, if you are interested in the full list of configurable parameters in Bactopia, please check out the [Complete Usage](usage-complete.md) section.

## Usage
```
Required Parameters:
    ### For Procesessing Multiple Samples
    --fastqs STR            An input file containing the sample name and
                                absolute paths to FASTQ/FASTAs to process

    ### For Processing A Single Sample
    --R1 STR                First set of reads for paired end in compressed (gzip)
                                FASTQ format

    --R2 STR                Second set of reads for paired end in compressed (gzip)
                                FASTQ format

    --SE STR                Single end set of reads in compressed (gzip) FASTQ format

    --hybrid                The SE should be treated as long reads for hybrid assembly.

    --sample STR            The name of the input sequences

    ### For Downloading from SRA/ENA or NCBI Assembly
    **Note: Assemblies will have error free Illumina reads simulated for processing.**
    --accessions            An input file containing ENA/SRA Experiment accessions or
                                NCBI Assembly accessions to be processed

    --accession             A single ENA/SRA Experiment accession or NCBI Assembly accession
                                to be processed

    ### For Processing an Assembly
    **Note: The assembly will have error free Illumina reads simulated for processing.**
    --assembly STR          A assembled genome in compressed FASTA format.

    --reassemble            The simulated reads will be used to create a new assembly.
                                Default: Use the original assembly, do not reassemble

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

    --max_time INT          The maximum number of minutes a single task should run before being halted.
                                Default: 240 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single task.
                                Default: 64 Gb

    --cpus INT              Number of processors made available to a single task.
                                Default: 4

    -qs                     Nextflow queue size. This parameter is very useful to limit the total number of
                                processors used on desktops, laptops or shared resources.
                                Default: Nextflow defaults to the total number of processors on your system.


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

    --publish_mode          Set Nextflow's method for publishing output files. Allowed methods are:
                                'copy' (default)    Copies the output files into the published directory.

                                'copyNoFollow' Copies the output files into the published directory
                                               without following symlinks ie. copies the links themselves.

                                'link'    Creates a hard link in the published directory for each
                                          process output file.

                                'rellink' Creates a relative symbolic link in the published directory
                                          for each process output file.

                                'symlink' Creates an absolute symbolic link in the published directory
                                          for each process output file.

                                Default: copy

    --force                 Nextflow will overwrite existing output files.
                                Default: false

    -resume                 Nextflow will attempt to resume a previous run. Please notice it is
                                only a single '-'

    --cleanup_workdir       After Bactopia is successfully executed, the work firectory will be deleted.
                                Warning: by doing this you lose the ability to resume workflows.

Useful Parameters:
    --skip_logs             Logs for each process per sample will not be kept.

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

    --version               Print workflow version information

    --help                  Show this message and exit

    --help_all              Show a complete list of adjustable parameters
```

## Inputs
Bactopia has multiple approaches to specify your input sequences. Bactopia can process Illumina FASTQs, assemblies, and long reads for hybrid assembly.

Illumina FASTQs a can be your local FASTQs or a Experiment accession to download associated FASTQs from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena). 

Likewise assemblies can be local, or a GenBank/RefSeq accession to download from NCBI Assembly. Input assemblies will have Illumina reads simulated so that the complete Bactopia pipeline run. By default, the assembly will not be reassembled.

Which approach really depends on what you need to achieve! The following sections describe methods to process single samples, multiple samples, downloading samples from the ENA.

### Local
#### Single Sample
When you only need to process a single sample at a time, Bactopia allows that! You only have to the sample name (`--sample`) and the whether the read set is paired-end (`--R1` and `--R2`), single-end (`--SE`), Illumina paired-end + long reads (`--hybrid`), or an assembly (`--assembly`).

!!! info "Use --R1, --R2 for Paired-End FASTQs"
    `bactopia --sample my-sample --R1 /path/to/my-sample_R1.fastq.gz --R2 /path/to/my-sample_R2.fastq.gz`

!!! info "Use --SE for Single-End FASTQs"
    `bactopia --sample my-sample --SE /path/to/my-sample.fastq.gz`

!!! info "Use --R1, --R2, --SE, and --hybrid for Paired-End FASTQs with Long Reads"
    At the assembly step, Unicycler will be used to create a hybrid assembly using the paired-end reads and the long reads.
    ```
    bactopia --sample my-sample 
             --R1 /path/to/my-sample_R1.fastq.gz \
             --R2 /path/to/my-sample_R2.fastq.gz \
             --SE /path/to/my-sample.fastq.gz \
             --hybrid
    ```

!!! info "Use --assembly for an assembled FASTA"
    Assemblies will have 2x250bp Illumina reads simulated without insertions or deletions in the sequence and a minimum PHRED score of Q33. By default, the input assembly will be used for all downstream analyses (e.g. annotation) which use an assembly. If the `--reassemble` parameter is given, then the a assembly will be created from the simulated reads.
    ```
    bactopia --sample my-sample --assembly /path/to/my-sample.fna.gz 

    ```

#### Multiple Samples
For multiple samples, you must create a file with information about the inputs, a *file of filenames* (FOFN). This file specifies sample names and location of FASTQs/FASTAs to be processed. Using this information, paired-end, single-end, hybrid or assembly information can be extracted as well as naming output files.

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

sample	runtype	r1	r2	extra
SA103113	assembly			/example/SA103113.fna.gz
SA110685	hybrid	/example/SA110685_R1.fastq.gz	/SA110685_R2.fastq.gz	/example/SA110685.fastq.gz
SA123186	paired-end	/example/SA123186_R1.fastq.gz	/example/SA123186_R2.fastq.gz
SA123456	single-end	/example/SA12345.fastq.gz
```

The expected structure is a **tab-delimited** table with three columns:

1. `sample`: A unique prefix, or unique name, to be used for naming output files
2. `runtype`: Informs Bactopia what type of input the sample is
3. `r1`: If paired-end, the first pair of reads, else the single-end reads
4. `r2`: If paired-end, the second pair of reads
5. `extra`: Either the assembly or long reads associated with a sample.

These five columns are used as the header for the file. In other words, all input FOFNs require their first line to be:
```
sample	runtype	r1	r2	extra
```

All lines after the header line, contain unique sample names and location(s) to associated FASTQ file(s). Absolute paths should be used to prevent any *file not found* errors due to the relative path changing.

In the example above, four samples would be processed by Bactopia. 

1. `SA103113` would have simulated reads crreated from the assembly
2. `SA110685` would have a hybrid assembly created using the paired-end reads and long-reads
3. `SA123186` would be processed as paired-end reads
4. `SA123456` would be processed as single-end reads

##### Generating A FOFN
`bactopia prepare` has been included to help aid (hopefully!) the process of creating a FOFN for your samples. This script will attempt to find FASTQ files in a given directory and output the expected FOFN format. It will also output any potential issues associated with the pattern matching.

!!! error "Verify accuracy of FOFN"
    This is currently an experimental function. There are likely bugs to be ironed out. Please be sure to give the resulting FOFN a quick look over.

###### Usage
```
usage: bactopia prepare [-h] [-f STR] [-a STR] [--fastq_seperator STR]
                        [--fastq_pattern STR] [--assembly_pattern STR]
                        [--long_reads] [--version]
                        STR

bactopia prepare - Read a directory and prepare a FOFN of FASTQs/FASTAs

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
  --pe1_pattern STR     Designates difference first set of paired-end reads.
                        Default: ([Aa]|[Rr]1) (R1, r1, 1, A, a)
  --pe2_pattern STR     Designates difference second set of paired-end reads.
                        Default: ([Bb]|[Rr]2) (R2, r2, 2, AB b)
  --assembly_pattern STR
                        Glob pattern to match assembly FASTAs. Default:
                        *.fna.gz
  -r, --recursive       Directories will be traversed recursively
  --long_reads          Single-end reads should be treated as long reads
  --merge               Flag samples with multiple read sets to be merged by
                        Bactopia
  --version             show program's version number and exit
```

##### Validating FOFN
When a FOFN is given, the first thing Bactopia does is verify all FASTQ files are found. If everything checks out, each sample will then be processed, otherwise a list of samples with errors will be output to STDERR. 

If you would like to only validate your FOFN (and not run the full pipeline), you can use the `--check_fastqs` parameter.

###### Without Errors
```
N E X T F L O W  ~  version 20.01.0
Launching `/home/rpetit3/repos/bactopia/main.nf` [gigantic_meitner] - revision: 6a0fbfbd9c
Printing what would have been processed. Each line consists of an array of
five elements: [SAMPLE_NAME, RUNTYPE, IS_SINGLE_END, [FASTQ_1, FASTQ_2], EXTRA]

Found:

[SA103113, assembly, false, [null, null], /example/SA103113.fna.gz]
[SA110685, hybrid, false, [/example/SA110685_R1.fastq.gz, /example/SA110685_R2.fastq.gz], /example/SA110685.fastq.gz]
[SA123186, paired-end, false, [/example/SA123186_R1.fastq.gz, /example/SA123186_R2.fastq.gz], null]
[SA12345, single-end, true, [/example/SA12345.fastq.gz], null]
```
Each sample has passed validation and is put into a three element array:

1. sample - the name for this sample
2. is_single_end - the reads are single-end (true) or paired-end (false)
3. fastq_array - the fastqs associated with the sample

This array is then automatically queued up for proccessing by Nextflow.

###### With errors
```
N E X T F L O W  ~  version 20.01.0
Launching `/home/rpetit3/repos/bactopia/main.nf` [special_ampere] - revision: 6a0fbfbd9c
LINE 4:ERROR: Please verify /example-bad/SA123186_R1.fastq.gz exists, and try again
LINE 4:ERROR: Please verify /example-bad/SA123186_R2.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /example-bad/SA12345.fastq.gz exists, and try again
Sample name "SA123186" is not unique, please revise sample names
Verify sample names are unique and/or FASTA/FASTQ paths are correct
See "--example_fastqs" for an example
Exiting
```

In the above example, there are multiple errors. Lines 4 and 5 (`LINE 4:ERROR` or `LINE 5:ERROR`) suggest that based on the given paths the FASTQs do not exist. The sample name `SA123186` has been used multiple times, and must be corrected.

### ENA & SRA
There are a lot of publicly avilable sequences available from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) (ENA) and the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA). There's a good chance you might want to include some of those sequences in your analysis! If that sounds like you, Bactopia has that built in for you! You can give a single *Experiment* accession (`--accession`) or a file where each line is a single *Experiment* accession (`--accessions`). Bactopia will then query ENA to determine *Run* accession(s) associated with the given Experiment accession and proceed download the corresponding FASTQ files from either the SRA (default) or ENA (`--use_ena`). 

After the download is completed, it will be processed through Bactopia.

!!! info "Use --accession for a Single Experiment Accession"
    SRA: `bactopia --accession SRX476958`  
    ENA: `bactopia --accession SRX476958 --use_ena`

!!! info "Use --accessions for Multiple Experiment Accessions"
    SRA: `bactopia --accessions my-accessions.txt`  
    ENA: `bactopia --accessions my-accessions.txt --use_ena`

!!! info "What happens when an Experiment has multiple Runs?"
    In cases where a single Experiment might have multiple Run accessions associated with it, the FASTQ files from each Run are merged into a single set of sequences.

#### Generating Accession List
`bactopia search` has been made to help assist in generating a list of Experiment accessions to be procesed by Bactopia (via `--accessions`). Users can provide a Taxon ID (e.g. 1280), a binary name (e.g. Staphylococcus aureus), a Study accession (e.g. PRJNA480016), a BioSample accession (e.g. SAMN01737350), or a Run accession (e.g. SRR578340). This value is then queried against ENA's [Data Warehouse API](https://www.ebi.ac.uk/ena/browse/search-rest)), and a list of all Experiment accessions associated with the query is returned.

##### Usage
```
usage: bactopia search [-h] [--exact_taxon] [--outdir OUTPUT_DIRECTORY]
                       [--prefix PREFIX] [--limit INT] [--version]
                       STR

bactopia search - Search ENA for associated WGS samples

positional arguments:
  STR                   Taxon ID or Study, BioSample, or Run accession

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
  bactopia search SAMN01737350
  bactopia search SRR578340
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

## `--cleanup_workdir`
After you run Bactopia, you will notice a directory called `work`. This directory is where Nextflow runs all the processes and stores the intermediate files. After a process completes successfully, the appropriate results are pulled out and placed in the sample's result folder. The `work` directory can grow very large very quickly! Please keep this in mind when using Bactopia. To help prevent the build up of the `work` directory you can use `--cleanup_workdir` to delete intermediate files after a successful execution of Bactopia.

!!! info "Bactopia and Bactopia Tools use separate `work` directories"
    Inside the `work` directory there will be separate subfolders that correspond to a Bactopia run or a specific Bactopia Tool run. This allows you to more easily identify which are ok to delete. The `work` directory is always ok to delete after a successful run.

## `--cpus`
At execution, Nextflow creates a queue and the number of slots in the queue is determined by the total number of cores on the system. So if you have a 24-core system, that means Nextflow will have a queue with 24-slots available. This feature kind of makes `--cpus` a little misleading. Typically when you give `--cpus` you are saying *"use this amount of cpus"*. But that is not the case for Nextflow and Bactopia. When you use `--cpus` what you are actually saying is *"for any particular task, use this amount of slots"*. Commands within a task processors will use the amount specified by `--cpus`.

!!! error "`--cpus` can have a significant effect on the efficiency of Bactopia"
    So for example if you have a system with 24-cores.

    This command, `bactopia ... --cpus 24`, says *for any particular task, use 24 slots*. Nextflow will give tasks in Bactopia 24 slots out of 24 available (24-core machine). In other words the queue can one have one task running at once because each task occupies 24 slots.

    On the other hand, `bactopia ... --cpus 4` says *for any particular task, use 4 slots*. Now, for Nextflow will give each task 4 slots out of 24 slots. Which means 6 tasks can be running at once. This can lead to much better efficiency because less jobs are stuck waiting in line. 

    There are some tasks in Bactopia that will only ever use a single slot because they are single-core tasks. But for example the `annotation` step will always use the number of slots specified by `--cpus`. If the `--cpus` is too high, the `annotation` will get bogged down, which causes tasks dependent on `annotation` to also get bogged down.

!!! info "When in doubt `--cpus 4` is a safe value."
    This is also the default value for Bactopia.

## `-qs`
The `-qs` parameter is short for *queue size*. As described above for `--cpus`, the default value for `-qs` is set to the total number of cores on the system. This parameter allows you to adjust the maximum number of cores Nextflow can use at any given moment.

!!! error "`-qs` allows you to play nicely on shared resources"
    From the example above, if you have a system with 24-cores. The default queue size if 24 slots.

    `bactopia ... --cpus 4` says *for any particular task, use a maximum of 4 slots*. Nextflow will give each task 4 slots out of 24 slots. But there might be other people also using the server.

    `bactopia ... --cpus 4 -qs 12` says *for any particular task, use a maximum of 4 slots, but don't use more than 12 slots*. Nextflow will give each task 4 slots out of 12 slots. Now instead of using all the cores on the server, the maximum that can be used in 12.

!!! info "`-qs` might need adjusting for job schedulers."
    The default value for `-qs` is set to 100 when using a job scheduler (e.g. SLURM, AWS Batch). There may be times when you need adjust this to meet your needs. For example, if using AWS Batch you might want to increase the value to have more jobs processed at once (e.g. 100 vs 500).


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
