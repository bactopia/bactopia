# Runtime Parameters

Bactopia includes numerous (100+) configurable parameters. Basically for each step of the pipeline, you can modify the default parameters of a specific tool.

## Required
The required parameters depends on how many samples are to be proccessed. You can learn more about which approach to take at [Specifying Input FASTQs](usage-basic.md#fastq-inputs).
```
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
```

## Dataset
If you followed the steps in [Build Datasets](datasets.md), you can use the following parameters to point Bactopia to you datasets.

```
    --datasets DIR          The path to available datasets that have
                                already been set up

    --species STR           Determines which species-specific dataset to
                                use for the input sequencing
```

## Optional
These optional parameters, while not required, will be quite useful to tweak.

```
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

    --max_time INT          The maximum number of minutes a task should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single task.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single task.
                                Default: 4

    -qs                     Nextflow queue size. This parameter is very useful to limit the total number of
                                processors used on desktops, laptops or shared resources.
                                Default: Nextflow defaults to the total number of processors on your system.
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

## `-qs`
The `-qs` parameter is short for *queue size*. As described above for `--cpus`, the default value for `-qs` is set to the total number of cores on the system. This parameter allows you to adjust the maximum number of cores Nextflow can use at any given moment.

!!! error "`-qs` allows you to play nicely on shared resources"
    From the example above, if you have a system with 24-cores. The default queue size if 24 slots.

    `bactopia ... --cpus 4` says *for any particular task, use a maximum of 4 slots*. Nextflow will give each task 4 slots out of 24 slots. But there might be other people also using the server.

    `bactopia ... --cpus 4 -qs 12` says *for any particular task, use a maximum of 4 slots, but don't use more than 12 slots*. Nextflow will give each task 4 slots out of 12 slots. Now instead of using all the cores on the server, the maximum that can be used in 12.

!!! info "`-qs` might need adjusting for job schedulers."
    The default value for `-qs` is set to 100 when using a job scheduler (e.g. SLURM, AWS Batch). There may be times when you need adjust this to meet your needs. For example, if using AWS Batch you might want to increase the value to have more jobs processed at once (e.g. 100 vs 500).


### `--genome_size`
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


## Helpers
The following parameters are useful to test input parameters.
```
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

## `--cleanup_workdir`
After you run Bactopia, you will notice a directory called `work`. This directory is where Nextflow runs all the processes and stores the intermediate files. After a process completes successfully, the appropriate results are pulled out and placed in the sample's result folder. The `work` directory can grow very large very quickly! Please keep this in mind when using Bactopia. To help prevent the build up of the `work` directory you can use `--cleanup_workdir` to delete intermediate files after a successful execution of Bactopia.

!!! info "Bactopia and Bactopia Tools use separate `work` directories"
    Inside the `work` directory there will be separate subfolders that correspond to a Bactopia run or a specific Bactopia Tool run. This allows you to more easily identify which are ok to delete. The `work` directory is always ok to delete after a successful run.

### `--nfconfig`
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

### `-resume`
Bactopia relies on [Nextflow's Resume Feature](https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume) to resume runs. You can tell Bactopia to resume by adding `-resume` to your command line. When `-resume` is used, Nextflow will review the cache and determine if the previous run is resumable. If the previous run is not resumable, execution will
start at the beginning.


### `--keep_all_files`
In some processes, Bactopia will delete large intermediate files (e.g. multiple uncompressed FASTQs) **only** after a process successfully completes. Since this a per-process function, it does not affect Nextflow's ability to resume (`-resume`)a workflow. You can deactivate this feature using `--keep_all_files`. Please, keep in mind the *work* directory is already large, this will make it 2-3 times larger.

## Additional Parameters
The remaining parameters are associated with specific programs. In the following sections, these parameters are grouped by which Nextflow process they are applicable to.

The description and default values for these parameters were taken from the program to which they apply.

It is important to note, not all of the available parameters for each and every program are available in Bactopia. If there is a parameter that was overlooked and should probably be included, please make a suggestion!

### ENA Download Parameters
```
ENA Download Parameters:
    --max_retry INT         Maximum times to retry downloads
                                Default: 10

    --use_ena               Download FASTQs from ENA with Aspera Connect.
                                Default: Download from SRA

    --ftp_only              If "--use_ena" is enabled, FTP will be used to
                                download FASTQs from ENA.

    --aspera_speed STR      Speed at which Aspera Connect will download.
                                Default: 100M

    --no_cache              Don't cache the assembly summary file from ncbi-genome-download
```

### FASTQ Minimum Requirements Parameters
```
FASTQ Minimum Requirements Parameters:
    --min_basepairs INT     The minimum amount of input sequenced basepairs required
                                to continue downstream analyses.
                                Default: 2241820

    --min_reads INT         The minimum amount of input sequenced reads required
                                to continue downstream analyses.
                                Default: 7472

    --min_proportion FLOAT  The minimum proportion of basepairs for paired-end reads to continue
                                downstream analyses. Example: If set to 0.75 the R1 and R2 must
                                have > 75% proportion of reads (e.g. R1 100bp, R2 75bp, not
                                R1 100bp, R2 50bp)
                                Default: 0.5

    --skip_fastq_check      The input FASTQs will not be check to verify they meet the
                                minimum requirements to be processed. This parameter
                                is useful if you are confident your sequences will
                                pass the minimum requirements.

```

### Estimate Genome Size Parameters
```
Estimate Genome Size Parameters:
    Only applied if the genome size is estimated.

    --min_genome_size INT   The minimum estimated genome size allowed for the input sequence
                            to continue downstream analyses.
                            Default: 100000

    --max_genome_size INT   The maximum estimated genome size allowed for the input sequence
                            to continue downstream analyses.
                            Default: 18040666

```

### QC Reads Parameters
```
QC Reads Parameters:
    --skip_qc               The QC step qill be skipped and it will be assumed the inputs
                                sequences have already been QCed.

    --skip_error_correction FLASH error correction of paired-end reads will be skipped.

    --qc_ram INT            Try to keep RAM usage below this many GB
                                Default: 3 GB

    --adapters FASTA        Illumina adapters to remove
                                Default: BBmap adapters

    --adapter_k INT         Kmer length used for finding adapters. Adapters
                                shorter than k will not be found
                                Default: 23

    --phix FASTA            phiX174 reference genome to remove
                                Default: NC_001422

    --phix_k INT            Kmer length used for finding phiX174.
                                Contaminants shorter than k will not be
                                found
                                Default: 31

    --ktrim STR             Trim reads to remove bases matching reference
                                kmers
                                Values:
                                    f (do not trim)
                                    r (trim to the right, Default)
                                    l (trim to the left)

    --mink INT              Look for shorter kmers at read tips down to this
                                length, when k-trimming or masking. 0 means
                                disabled. Enabling this will disable
                                maskmiddle
                                Default: 11

    --hdist INT             Maximum Hamming distance for ref kmers (subs only)
                                Memory use is proportional to (3*K)^hdist
                                Default: 1

    --tpe BOOL              When kmer right-trimming, trim both reads to the
                                minimum length of either
                                Values:
                                    f (do not equally trim)
                                    t (equally trim to the right, Default)

    --tbo BOOL              Trim adapters based on where paired reads overlap
                                Values:
                                    f (do not trim by overlap)
                                    t (trim by overlap, Default)

    --qtrim STR             Trim read ends to remove bases with quality
                                below trimq. Performed AFTER looking for
                                kmers
                                Values:
                                    rl (trim both ends, Default)
                                    f (neither end)
                                    r (right end only)
                                    l (left end only)
                                    w (sliding window)

    --trimq FLOAT           Regions with average quality BELOW this will be
                                trimmed if qtrim is set to something other
                                than f
                                Default: 6

    --maq INT               Reads with average quality (after trimming)
                                below this will be discarded
                                Default: 10

    --minlength INT         Reads shorter than this after trimming will be
                                discarded. Pairs will be discarded if both
                                are shorter
                                Default: 35

    --ftm INT               If positive, right-trim length to be equal to
                                zero, modulo this number
                                Default: 5

    --tossjunk              Discard reads with invalid characters as bases
                                Values:
                                    f (keep all reads)
                                    t (toss reads with ambiguous bases, Default)

    --qout STR              Output quality offset
                                Values:
                                    33 (PHRED33 offset quality scores, Default)
                                    64 (PHRED64 offset quality scores)
                                    auto (keeps the current input offset)

    --xmx STR               This will be passed to Java to set memory usage
                                Examples:
                                    '8g' will specify 8 gigs of RAM (Default)
                                    '20g' will specify 20 gigs of RAM
                                    '200m' will specify 200 megs of RAM

    --maxcor INT            Max number of corrections within a 20bp window
                                Default: 1

    --sampleseed INT        Set to a positive number to use as the rng seed
                                for sampling
                                Default: 42
```

### Assembly Parameters
```
Assembly Parameters:
    Standard Assembly:
    --shovill_ram INT       Try to keep RAM usage below this many GB
                                Default: 8 GB

    --assembler STR         Assembler: megahit velvet skesa spades
                                Default: skesa

    --min_contig_len INT    Minimum contig length <0=AUTO>
                                Default: 500

    --min_contig_cov INT    Minimum contig coverage <0=AUTO>
                                Default: 2

    --contig_namefmt STR    Format of contig FASTA IDs in 'printf' style
                                Default: contig%05d

    --shovill_opts STR      Extra assembler options in quotes eg.
                                spades: "--untrusted-contigs locus.fna" ...
                                Default: ''

    --shovill_kmers STR     K-mers to use <blank=AUTO>
                                Default: ''

    --trim                  Enable adaptor trimming

    --nostitch              Disable read stitching

    --nocorr                Disable post-assembly correction

    Hybrid Assembly:
    --unicycler_ram INT       Try to keep RAM usage below this many GB
                                Default: 32 GB

    --unicycler_mode STR    Bridging mode used by Unicycler, choices are:
                                conservative = smaller contigs, lowest
                                               misassembly rate
                                normal = moderate contig size and
                                         misassembly rate (Default)
                                bold = longest contigs, higher misassembly
                                       rate

    --min_polish_size INT   Contigs shorter than this value (bp) will not be
                                polished using Pilon
                                Default: 10000

    --min_component_size INT
                            Graph components smaller than this size (bp) will
                                be removed from the final graph
                                Default: null

    --min_dead_end_size INT
                            Graph dead ends smaller than this size (bp) will
                                be removed from the final graph
                                Default: 1000

    --no_miniasm            Skip miniasm+Racon bridging
                                Default: Produce long-read bridges

    --no_rotate             Do not rotate completed replicons to start at a
                                standard gene

    --no_pilon              Do not use Pilon to polish the final assembly
```

### Assembly Quality Control Parameters
```
Assembly Quality Control Parameters:
    --checkm_unique INT     Minimum number of unique phylogenetic markers required
                                to use lineage-specific marker set.
                                Default: 10

    --checkm_multi INT      Maximum number of multi-copy phylogenetic markers before
                                defaulting to domain-level marker set.
                                Default: 10

    --aai_strain FLOAT      AAI threshold used to identify strain heterogeneity
                                Default: 0.9

    --checkm_length FLOAT   Percent overlap between target and query
                                Default: 0.7

    --full_tree             Use the full tree (requires ~40GB of memory) for determining
                                lineage of each bin.
                                Default: Use reduced tree (<16gb memory)

    --skip_pseudogene_correction
                            Skip identification and filtering of pseudogene

    --ignore_thresholds     Ignore model-specific score thresholds

    --checkm_ali            Generate HMMER alignment file for each bin

    --checkm_nt             Generate nucleotide gene sequences for each bin

    --force_domain          Use domain-level sets for all bins

    --no_refinement         Do not perform lineage-specific marker set refinement

    --individual_markers    Treat marker as independent (i.e., ignore co-located
                                set structure.

    --skip_adj_correction   Do not exclude adjacent marker genes when estimating
                                contamination

    --contig_thresholds STR Comma-separated list of contig length thresholds
                                Default: 0,1000,10000,100000,250000,1000000

    --plots_format STR      Save plots in specified format.
                                Supported formats: emf, eps, pdf, png, ps, raw,
                                                    rgba, svg, svgz
                                Default: pdf
```

### Count 31mers Parameters
```
Count 31mers Parameters:
    --cortex_ram INT        Try to keep RAM usage below this many GB
                                Default: 8 GB

    --keep_singletons       Keep all counted 31-mers
                                Default: Filter out singletons
```

### Annotation Parameters
```
Annotation Parameters:
    --compliant             Force Genbank/ENA/DDJB compliance: --genes --mincontiglen 500 --centre 'Bactopia'
                                Default: false

    --centre STR            Sequencing centre ID
                                Default: 'Bactopia'

    --addmrna               Add 'mRNA' features for each 'CDS' feature

    --rawproduct            Do not clean up /product annotation

    --cdsrnaolap            Allow [tr]RNA to overlap CDS

    --prokka_evalue STR     Similarity e-value cut-off
                                Default: 1e-09

    --prokka_coverage INT   Minimum coverage on query protein
                                 Default: 80

    --nogenes               Do not add 'gene' features for each 'CDS' feature

    --norrna                Don't run rRNA search

    --notrna                Don't run tRNA search

    --rnammer               Prefer RNAmmer over Barrnap for rRNA prediction

    --rfam                  Enable searching for ncRNAs with Infernal+Rfam

    --skip_prodigal_tf      If a Prodigal training file was found, it will not be used
```

### Minmer Sketch Parameters
```
Minmer Sketch Parameters:
    --mash_sketch INT       Sketch size. Each sketch will have at most this
                                many non-redundant min-hashes.
                                Default: 10000

    --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                input k-mers
                                Default: 10000
```

### Minmer Query Parameters
```
Minmer Query Parameters:
    --minmer_ram INT        Try to keep RAM usage below this many GB
                                Default: 5 GB

    --screen_w              Winner-takes-all strategy for identity estimates.
                                After counting hashes for each query, hashes
                                that appear in multiple queries will be
                                removed from all except the one with the best
                                identity (ties broken by larger query), and
                                other identities will be reduced. This
                                removes output redundancy, providing a rough
                                compositional outline.
                                Default: True

    --screen_i FLOAT        Minimum identity to report. Inclusive unless set
                                to zero, in which case only identities greater
                                than zero (i.e. with at least one shared hash)
                                will be reported. Set to -1 to output
                                everything.
                                Default: 0.8
```

### Ariba Parameters
```
Ariba Parameters:
    --nucmer_min_id INT     Minimum alignment identity (delta-filter -i)
                                Default: 90

    --nucmer_min_len INT    Minimum alignment length (delta-filter -i)
                                Default: 20

    --nucmer_breaklen INT   Value to use for -breaklen when running nucmer
                                Default: 200

    --assembly_cov INT      Target read coverage when sampling reads for
                                assembly
                                Default: 50

    --min_scaff_depth INT   Minimum number of read pairs needed as evidence
                                for scaffold link between two contigs
                                Default: 10

    --spades_options STR    Extra options to pass to Spades assembler
                                Default: null

    --assembled_threshold FLOAT (between 0 and 1)
                            If proportion of gene assembled (regardless of
                                into how many contigs) is at least this
                                value then the flag gene_assembled is set
                                Default: 0.95

    --gene_nt_extend INT    Max number of nucleotides to extend ends of gene
                                matches to look for start/stop codons
                                Default: 30

    --unique_threshold FLOAT (between 0 and 1)
                            If proportion of bases in gene assembled more
                                than once is <= this value, then the flag
                                unique_contig is set
                                Default: 0.03

    --ariba_no_clean        Do not clean up intermediate files created by
                                Ariba. By default, the local assemblies are
                                deleted.
```

### Call Variant Parameters
```
Call Variant Parameters:
    --snippy_ram INT        Try and keep RAM under this many GB
                                Default: 4 GB

    --mapqual INT           Minimum read mapping quality to consider
                                Default: 60

    --basequal INT          Minimum base quality to consider
                                Default: 13

    --mincov INT            Minimum site depth to for calling alleles
                                Default: 10

    --minfrac FLOAT         Minimum proportion for variant evidence (0=AUTO)
                                Default: 0

    --minqual INT           Minimum QUALITY in VCF column 6
                                Default: 100

    --maxsoft INT           Maximum soft clipping to allow
                                Default: 10

    --bwaopt STR            Extra BWA MEM options, eg. -x pacbio
                                Default: ''

    --fbopt STR             Extra Freebayes options,
                                eg. --theta 1E-6 --read-snp-limit 2
                                Default: ''
```

### Nearest Neighbor Reference Genomes
```
Nearest Neighbor Reference Genomes:
    --max_references INT    Maximum number of nearest neighbor reference genomes to
                                download for variant calling.
                                Default: 1

    --random_tie_break      On references with matching distances, randomly select one.
                                Default: Pick earliest accession number

    --disable_auto_variants Disable automatic selection of reference genome based on
                                Mash distances.
```

### BLAST Parameters
```
BLAST Parameters:
    --perc_identity INT     Percent identity
                                Default: 50

    --qcov_hsp_perc INT     Percent query coverage per hsp
                                Default: 50

    --max_target_seqs INT   Maximum number of aligned sequences to
                                keep
                                Default: 2000
```

### Mapping Parameters
```
Mapping Parameters:
    --keep_unmapped_reads   Keep unmapped reads, this does not affect variant
                                calling.

    --bwa_mem_opts STR      Extra BWA MEM options
                                Default: ''

    --bwa_aln_opts STR      Extra BWA ALN options
                                Default: ''

    --bwa_samse_opts STR    Extra BWA SAMSE options
                                Default: ''

    --bwa_sampe_opts STR    Extra BWA SAMPE options
                                Default: ''

    --bwa_n INT             Maximum number of alignments to output in the XA
                                tag for reads paired properly. If a read has
                                more than INT hits, the XA tag will not be
                                written.
                                Default: 9999
```

### Antimicrobial Resistance Parameters
```
Antimicrobial Resistance Parameters:
    --update_amr            Force amrfinder to update its database.

    --amr_ident_min         Minimum identity for nucleotide hit (0..1). -1
                                means use a curated threshold if it exists and
                                0.9 otherwise
                                Default: -1

    --amr_coverage_min      Minimum coverage of the reference protein (0..1)
                                Default: 0.5

    --amr_organism          Taxonomy group: Campylobacter, Escherichia, Klebsiella
                                Salmonella, Staphylococcus, Vibrio
                                Default: ''

    --amr_translation_table NCBI genetic code for translated BLAST
                                Default: 11

    --amr_plus              Add the plus genes to the report

    --amr_report_common     Suppress proteins common to a taxonomy group
```
