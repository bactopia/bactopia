/*
  Bactopia Usage
*/

def print_usage() {
    usage_text = params.help_all ? full_help() : basic_help()
    log.info"""
    ${PROGRAM_NAME} v${VERSION}
    ${basic_help()}
    ${params.help_all ? full_help() : ""}
    """.stripIndent()

    if (params.conda_help) {
        // Cleanup up the directory
        // This is only meant to be used with tests for conda build
        file("./work/").deleteDir()
        file("./.nextflow/").deleteDir()
        def files = new FileNameByRegexFinder().getFileNames('./', '.nextflow.log*')
        files.each { new File(it).delete()}
    }
    exit 0
}

def basic_help() {
    genome_size = params.genome_size ? params.genome_size : "Mash Estimate"
    return """
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
        --assembly STR          A assembled genome in compressed (gzip) FASTA format.

        --reassemble            The simulated reads will be used to create a new assembly.
                                    Default: Use the original assembly, do not reassemble

    Dataset Parameters:
        --datasets DIR          The path to available datasets that have
                                    already been set up

        --species STR           Determines which species-specific dataset to
                                    use for the input sequencing

    Optional Parameters:
        --coverage INT          Reduce samples to a given coverage
                                    Default: ${params.coverage}x

        --genome_size INT       Expected genome size (bp) for all samples, a value of '0'
                                    will disable read error correction and read subsampling.
                                    Special values (requires --species):
                                        'min': uses minimum completed genome size of species
                                        'median': uses median completed genome size of species
                                        'mean': uses mean completed genome size of species
                                        'max': uses max completed genome size of species
                                    Default: Mash estimate

        --outdir DIR            Directory to write results to
                                    Default: ${params.outdir}

    Nextflow Queue Parameters:
        At execution, Nextflow creates a queue and the number of slots in the queue is determined by the total number
        of cores on the system. When a task is submitted to the queue, the total number of slots it occupies is
        determined by the value set by "--cpus".

        This can have a significant effect on the efficiency of the Nextflow's queue system. If "--cpus" is set to a
        value that is equal to the number of cores availabe, in most cases only a single task will be able to run
        because its occupying all available slots.

        When in doubt, "--cpus 4" is a safe bet, it is also the default value if you don't use "--cpus".

        --max_retry INT         Maximum times to retry a process before allowing it to fail.
                                    Default: ${params.max_retry}

        --min_time INT          The minimum number of minutes a single task should run before being halted.
                                    Default: ${params.min_time} minutes

        --max_time INT          The maximum number of minutes a single task should run before being halted.
                                    Default: ${params.max_time} minutes

        --max_memory INT        The maximum amount of memory (Gb) allowed to a single task.
                                    Default: ${params.max_memory} Gb

        --cpus INT              Number of processors made available to a single task.
                                    Default: ${params.cpus}

        -qs INT                 Nextflow queue size. This parameter is very useful to limit the total number of
                                    processors used on desktops, laptops or shared resources.
                                    Default: Nextflow defaults to the total number of processors on your system.


    Nextflow Related Parameters:
        --infodir DIR           Directory to write Nextflow summary files to
                                    Default: ${params.infodir}

        --condadir DIR          Directory to Nextflow should use for Conda environments
                                    Default: Bactopia's Nextflow directory

        --registry STR          Docker registry to pull containers from.
                                    Available options: dockerhub, quay, or github
                                    Default: dockerhub

        --singularity_cache STR Directory where remote Singularity images are stored. If using a cluster, it must
                                    be accessible from all compute nodes.
                                    Default: NXF_SINGULARITY_CACHEDIR evironment variable, otherwise ${params.singularity_cache}

        --queue STR             The name of the queue(s) to be used by a job scheduler (e.g. AWS Batch or SLURM).
                                    If using multiple queues, please seperate queues by a comma without spaces.
                                    Default: ${params.queue}

        --disable_scratch       All intermediate files created on worker nodes of will be transferred to the head node.
                                    Default: Only result files are transferred back

        --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows
                                    you to create profiles specific to your environment (e.g. SGE,
                                    AWS, SLURM, etc...). This config file is loaded last and will
                                    overwrite existing variables if set.
                                    Default: Bactopia's default configs

        --nfdir                 Print directory Nextflow has pulled Bactopia to

        --overwrite             Nextflow will overwrite existing output files.
                                    Default: ${params.overwrite}

        --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                    will wait before execution.
                                    Default: ${params.sleep_time} seconds

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

                                    Default: ${params.publish_mode}

        --force                 Nextflow will overwrite existing output files.
                                    Default: ${params.force}

        -resume                 Nextflow will attempt to resume a previous run. Please notice it is
                                    only a single '-'

        --cleanup_workdir       After Bactopia is successfully executed, the work directory will be deleted.
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
    """
}


def full_help() {
    return """
    Additional Parameters:
    The description of the following parameters were taken from the program for
    which they apply to.

    Many of the default values were also taken from the program for which they
    apply to.

    AWS Batch Profile  (-profile awsbatch) Parameters:
        --aws_region STR        AWS Region to be used by Nextflow
                                    Default: ${params.aws_region}

        --aws_volumes STR       Volumes to be mounted from the EC2 instance to the Docker container
                                    Default: ${params.aws_volumes}

        --aws_cli_path STR       Path to the AWS CLI for Nextflow to use.
                                    Default: ${params.aws_cli_path}

        --aws_upload_storage_class STR
                                The S3 storage slass to use for storing files on S3
                                    Default: ${params.aws_upload_storage_class}

        --aws_max_parallel_transfers INT
                                The number of parallele transfers between EC2 and S3
                                    Default: ${params.aws_max_parallel_transfers}

        --aws_delay_between_attempts INT
                                The duration of sleep (in seconds) between each transfer between EC2 and S3
                                    Default: ${params.aws_delay_between_attempts}

        --aws_max_transfer_attempts INT
                                The maximum number of times to retry transferring a file between EC2 and S3
                                    Default: ${params.aws_max_transfer_attempts}

        --aws_max_retry INT     The maximum number of times to retry a process on AWS Batch
                                    Default: ${params.aws_max_retry}

        --aws_ecr_registry STR  The ECR registry containing Bactopia related containers.
                                    Default: Use the registry given by --registry

    ENA/SRA Download Parameters:
        --max_downloads INT     Maximum number of FASTQs to download at once.
                                    Warning: Setting this value too high can lead to NCBI temporarily
                                             blocking your IP addess. 3-5 is reasonable, >10 is likely
                                             to be excessive.
                                    Default: ${params.max_downloads}

        --use_ena               Download FASTQs from ENA with Aspera Connect.
                                    Default: Download from SRA

        --no_cache              Don't cache the assembly summary file from ncbi-genome-download

    FASTQ Minimum Requirements Parameters:
        --min_basepairs INT     The minimum amount of input sequenced basepairs required
                                    to continue downstream analyses.
                                    Default: ${params.min_basepairs}

        --min_coverage INT      The minimum coverage of input sequences required
                                    to continue downstream analyses.
                                    Default: ${params.min_coverage}

        --min_reads INT         The minimum amount of input sequenced reads required
                                    to continue downstream analyses.
                                    Default: ${params.min_reads}

        --min_proportion FLOAT  The minimum proportion of basepairs for paired-end reads to continue
                                    downstream analyses. Example: If set to 0.75 the R1 and R2 must
                                    have > 75% proportion of reads (e.g. R1 100bp, R2 75bp, not
                                    R1 100bp, R2 50bp)
                                    Default: ${params.min_proportion}

        --skip_fastq_check      The input FASTQs will not be check to verify they meet the
                                    minimum requirements to be processed. This parameter
                                    is useful if you are confident your sequences will
                                    pass the minimum requirements.

    Estimate Genome Size Parameters:
        Only applied if the genome size is estimated.

        --min_genome_size INT   The minimum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.min_genome_size}

        --max_genome_size INT   The maximum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.max_genome_size}

    QC Reads Parameters:
        --skip_qc               The QC step qill be skipped and it will be assumed the inputs
                                    sequences have already been QCed.

        --skip_error_correction FLASH error correction of paired-end reads will be skipped.

        --qc_ram INT            Try to keep RAM usage below this many GB
                                    Default: ${params.qc_ram} GB

        --adapters FASTA        Illumina adapters to remove
                                    Default: BBmap adapters

        --adapter_k INT         Kmer length used for finding adapters. Adapters
                                    shorter than k will not be found
                                    Default: ${params.adapter_k}

        --phix FASTA            phiX174 reference genome to remove
                                    Default: NC_001422

        --phix_k INT            Kmer length used for finding phiX174.
                                    Contaminants shorter than k will not be
                                    found
                                    Default: ${params.phix_k}

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
                                    Default: ${params.mink}

        --hdist INT             Maximum Hamming distance for ref kmers (subs only)
                                    Memory use is proportional to (3*K)^hdist
                                    Default: ${params.hdist}

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
                                    Default: ${params.trimq}

        --maq INT               Reads with average quality (after trimming)
                                    below this will be discarded
                                    Default: ${params.maq}

        --minlength INT         Reads shorter than this after trimming will be
                                    discarded. Pairs will be discarded if both
                                    are shorter
                                    Default: ${params.minlength}

        --ftm INT               If positive, right-trim length to be equal to
                                    zero, modulo this number
                                    Default: ${params.ftm}

        --tossjunk              Discard reads with invalid characters as bases
                                    Values:
                                        f (keep all reads)
                                        t (toss reads with ambiguous bases, Default)

        --qout STR              Output quality offset
                                    Values:
                                        33 (PHRED33 offset quality scores, Default)
                                        64 (PHRED64 offset quality scores)
                                        auto (keeps the current input offset)

        --maxcor INT            Max number of corrections within a 20bp window
                                    Default: ${params.maxcor}

        --sampleseed INT        Set to a positive number to use as the rng seed
                                    for sampling
                                    Default: ${params.sampleseed}


    Assembly Parameters:
        Standard Assembly:
        --shovill_ram INT       Try to keep RAM usage below this many GB
                                    Default: ${params.shovill_ram} GB

        --assembler STR         Assembler: megahit velvet skesa spades unicycler
                                    Default: ${params.assembler}

        --min_contig_len INT    Minimum contig length <0=AUTO>
                                    Default: ${params.min_contig_len}

        --min_contig_cov INT    Minimum contig coverage <0=AUTO>
                                    Default: ${params.min_contig_cov}

        --contig_namefmt STR    Format of contig FASTA IDs in 'printf' style
                                    Default: "SAMPLE_NAME_%05d"

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
                                    Default: ${params.unicycler_ram} GB

        --unicycler_mode STR    Bridging mode used by Unicycler, choices are:
                                    conservative = smaller contigs, lowest
                                                   misassembly rate
                                    normal = moderate contig size and
                                             misassembly rate (Default)
                                    bold = longest contigs, higher misassembly
                                           rate

        --min_polish_size INT   Contigs shorter than this value (bp) will not be
                                    polished using Pilon
                                    Default: ${params.min_polish_size}

        --min_component_size INT
                                Graph components smaller than this size (bp) will
                                    be removed from the final graph
                                    Default: ${params.min_component_size}

        --min_dead_end_size INT
                                Graph dead ends smaller than this size (bp) will
                                    be removed from the final graph
                                    Default: ${params.min_dead_end_size}

        --no_miniasm            Skip miniasm+Racon bridging
                                    Default: Produce long-read bridges

        --no_rotate             Do not rotate completed replicons to start at a
                                    standard gene

        --no_pilon              Do not use Pilon to polish the final assembly

    Assembly Quality Control Parameters:
        --skip_checkm           CheckM analysis will be skipped. This is useful for systems
                                    with less than 8GB of memory.

        --checkm_unique INT     Minimum number of unique phylogenetic markers required
                                    to use lineage-specific marker set.
                                    Default: ${params.checkm_unique}

        --checkm_multi INT      Maximum number of multi-copy phylogenetic markers before
                                    defaulting to domain-level marker set.
                                    Default: ${params.checkm_multi}

        --aai_strain FLOAT      AAI threshold used to identify strain heterogeneity
                                    Default: ${params.aai_strain}

        --checkm_length FLOAT   Percent overlap between target and query
                                    Default: ${params.checkm_length}

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
                                    Default: ${params.contig_thresholds}

        --plots_format STR      Save plots in specified format.
                                    Supported formats: emf, eps, pdf, png, ps, raw,
                                                        rgba, svg, svgz
                                    Default: ${params.plots_format}

    Count 31mers Parameters:
        --cortex_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.cortex_ram} GB

        --keep_singletons       Keep all counted 31-mers
                                    Default: Filter out singletons

    Annotation Parameters:
        --compliant             Force Genbank/ENA/DDJB compliance: --genes --mincontiglen ${params.min_contig_len} --centre '${params.centre}'
                                    Default: ${params.compliant}

        --centre STR            Sequencing centre ID
                                    Default: '${params.centre}'

        --addmrna               Add 'mRNA' features for each 'CDS' feature

        --rawproduct            Do not clean up /product annotation

        --cdsrnaolap            Allow [tr]RNA to overlap CDS

        --prokka_evalue STR     Similarity e-value cut-off
                                    Default: ${params.prokka_evalue}

        --prokka_coverage INT   Minimum coverage on query protein
                                     Default: ${params.prokka_coverage}

        --nogenes               Do not add 'gene' features for each 'CDS' feature

        --norrna                Don't run rRNA search

        --notrna                Don't run tRNA search

        --rnammer               Prefer RNAmmer over Barrnap for rRNA prediction

        --rfam                  Enable searching for ncRNAs with Infernal+Rfam

        --skip_prodigal_tf      If a Prodigal training file was found, it will not be used

    Minmer Sketch Parameters:
        --mash_sketch INT       Sketch size. Each sketch will have at most this
                                    many non-redundant min-hashes.
                                    Default: ${params.mash_sketch}

        --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                    input k-mers
                                    Default: ${params.sourmash_scale}


    Minmer Query Parameters:
        --minmer_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.minmer_ram} GB

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
                                    Default: ${params.screen_i}

    Ariba Parameters:
        --nucmer_min_id INT     Minimum alignment identity (delta-filter -i)
                                    Default: ${params.nucmer_min_id}

        --nucmer_min_len INT    Minimum alignment length (delta-filter -i)
                                    Default: ${params.nucmer_min_len}

        --nucmer_breaklen INT   Value to use for -breaklen when running nucmer
                                    Default: ${params.nucmer_breaklen}

        --assembly_cov INT      Target read coverage when sampling reads for
                                    assembly
                                    Default: ${params.assembly_cov}

        --min_scaff_depth INT   Minimum number of read pairs needed as evidence
                                    for scaffold link between two contigs
                                    Default: ${params.min_scaff_depth}

        --spades_options STR    Extra options to pass to Spades assembler
                                    Default: ${params.spades_options}

        --assembled_threshold FLOAT (between 0 and 1)
                                If proportion of gene assembled (regardless of
                                    into how many contigs) is at least this
                                    value then the flag gene_assembled is set
                                    Default: ${params.assembled_threshold}

        --gene_nt_extend INT    Max number of nucleotides to extend ends of gene
                                    matches to look for start/stop codons
                                    Default: ${params.gene_nt_extend}

        --unique_threshold FLOAT (between 0 and 1)
                                If proportion of bases in gene assembled more
                                    than once is <= this value, then the flag
                                    unique_contig is set
                                    Default: ${params.unique_threshold}

        --ariba_no_clean        Do not clean up intermediate files created by
                                    Ariba. By default, the local assemblies are
                                    deleted.

    Call Variant Parameters:
        --snippy_ram INT        Try and keep RAM under this many GB
                                    Default: ${params.snippy_ram} GB

        --mapqual INT           Minimum read mapping quality to consider
                                    Default: ${params.mapqual}

        --basequal INT          Minimum base quality to consider
                                    Default: ${params.basequal}

        --mincov INT            Minimum site depth to for calling alleles
                                    Default: ${params.mincov}

        --minfrac FLOAT         Minimum proportion for variant evidence (0=AUTO)
                                    Default: ${params.minfrac}

        --minqual INT           Minimum QUALITY in VCF column 6
                                    Default: ${params.minqual}

        --maxsoft INT           Maximum soft clipping to allow
                                    Default: ${params.maxsoft}

        --bwaopt STR            Extra BWA MEM options, eg. -x pacbio
                                    Default: ''

        --fbopt STR             Extra Freebayes options,
                                    eg. --theta 1E-6 --read-snp-limit 2
                                    Default: ''

    Nearest Neighbor Reference Genomes:
        --max_references INT    Maximum number of nearest neighbor reference genomes to
                                    download for variant calling.
                                    Default: 1

        --random_tie_break      On references with matching distances, randomly select one.
                                    Default: Pick earliest accession number

        --disable_auto_variants Disable automatic selection of reference genome based on
                                    Mash distances. This will not skip reference genomes
                                    you may have provided in the optional folders.

    BLAST Parameters:
        --perc_identity INT     Percent identity
                                    Default: ${params.perc_identity}

        --qcov_hsp_perc INT     Percent query coverage per hsp
                                    Default: ${params.qcov_hsp_perc}

        --max_target_seqs INT   Maximum number of aligned sequences to
                                    keep
                                    Default: ${params.max_target_seqs}

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
                                    Default: ${params.bwa_n}

    Antimicrobial Resistance Parameters:
        --skip_amr              AMRFinder+ analysis will be skipped. This is useful
                                    if the AMRFinder+ software and database versions are
                                    no longer compatible.

        --amr_ident_min         Minimum identity for nucleotide hit (0..1). -1
                                    means use a curated threshold if it exists and
                                    0.9 otherwise
                                    Default: ${params.amr_ident_min}

        --amr_coverage_min      Minimum coverage of the reference protein (0..1)
                                    Default: ${params.amr_coverage_min}

        --amr_organism          Taxonomy group: Campylobacter, Escherichia, Klebsiella
                                    Salmonella, Staphylococcus, Vibrio
                                    Default: ''

        --amr_translation_table NCBI genetic code for translated BLAST
                                    Default: ${params.amr_translation_table}

        --amr_plus              Add the plus genes to the report

        --amr_report_common     Suppress proteins common to a taxonomy group

    """
}
