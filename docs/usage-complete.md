# Runtime Parameters

Bactopia includes numerous (100+) configurable parameters. Basically for each step of the pipeline, you can modify the default parameters of a specific tool.

## Required
The required parameters depends on how many samples are to be proccessed. You can learn more about which approach to take at [Specifying Input FASTQs](usage-basic.md#fastq-inputs).
```
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
```

## Dataset
If you followed the steps in [Build Datasets](datasets.md), you can use the following parameters to point Bactopia to you datasets.

```
    --datasets DIR          The path to available public datasets that have
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

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process. 
                                Default: 4
```

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
The following parameters are useful to test to test input parameters.
```
    --infodir DIR           Directory to write Nextflow summary files to
                                Default: ./bactopia-info

    --condadir DIR          Directory to Nextflow should use for Conda environments
                                Default: Bactopia's Nextflow directory

    --nfdir                 Print directory Nextflow has pulled Bactopia to

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

### `--keep_all_files`
In some processes, Bactopia will delete large intermediate files (e.g. multiple uncompressed FASTQs) **only** after a process successfully completes. Since this a per-process function, it does not affect Nextflow's ability to resume (`-resume`)a workflow. You can deactivate this feature using `--keep_all_files`. Please, keep in mind the *work* directory is already large, this will make it 2-3 times larger.

## Program Specific
The remaining parameters are associated with specific programs. In the following sections, these parameters are grouped by which Nextflow process they are applicable to.

The description and default values for these parameters were taken from the program to which they apply.

It is important to note, not all of the available parameters for each and every program are available in Bactopia. If there is a parameter that was overlooked and should probably be included, please make a suggestion!

### Annotation
```
    --centre STR            Sequencing centre ID
                                Default: ''

    --addgenes              Add 'gene' features for each 'CDS' feature

    --addmrna               Add 'mRNA' features for each 'CDS' feature

    --rawproduct            Do not clean up /product annotation

    --cdsrnaolap            Allow [tr]RNA to overlap CDS

    --prokka_evalue STR     Similarity e-value cut-off
                                Default: 1e-09

    --prokka_coverage INT   Minimum coverage on query protein
                                 Default: 80

    --norrna                Don't run rRNA search

    --notrna                Don't run tRNA search

    --rnammer               Prefer RNAmmer over Barrnap for rRNA prediction
```

### Antimicrobial Resistance
```
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

### Ariba
```
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


### Assembly
```
    --shovill_ram INT       Try to keep RAM usage below this many GB
                                Default: 6 GB

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
```


### BLAST
```
    --perc_identity INT     Percent identity
                                Default: 50

    --qcov_hsp_perc INT     Percent query coverage per hsp
                                Default: 50

    --max_target_seqs INT   Maximum number of aligned sequences to
                                keep
                                Default: 2000

    --outfmt STR            BLAST alignment view options
                                Default: '6 qseqid qlen qstart qend sseqid slen sstart send length evalue bitscore pident nident mismatch gaps qcovs qcovhsp'
```


### Counting 31mers
```
    --cortex_ram INT        Try to keep RAM usage below this many GB
                                Default: 8 GB

    --keep_singletons       Keep all counted 31-mers
                                Default: Filter out singletons
```

### Download FASTQ
```
    --aspera_speed STR      Speed at which Aspera Connect will download.
                                Default: 100M

    --max_retry INT         Maximum times to retry downloads
                                Default: 10

    --ftp_only              Only use FTP to download FASTQs from ENA
```

### Download Reference Genome
```
    --max_references INT    Maximum number of nearest neighbor reference genomes to
                                download for variant calling.
                                Default: 1

    --random_tie_break      On references with matching distances, randomly select one.
                                Default: Pick earliest accession number

    --disable_auto_variants Disable automatic selection of reference genome based on
                                Mash distances.
```

### Insertion Mapping
```
    --min_clip INT          Minimum size for softclipped region to be
                                extracted from initial mapping
                                Default: 10

    --max_clip INT          Maximum size for softclipped regions to be
                                included
                                Default: 30

    --cutoff INT            Minimum depth for mapped region to be kept in
                                bed file
                                Default: 6

    --novel_gap_size INT    Distance in base pairs between left and right
                                flanks to be called a novel hit
                                Default: 15

    --min_range FLOAT       Minimum percent size of the gap to be called a
                                known hit
                                Default: 0.9

    --max_range FLOAT       Maximum percent size of the gap to be called a
                                known hit
                                Default: 1.1

    --merging INT           Value for merging left and right hits in bed
                                files together to simply calculation of
                                closest and intersecting regions
                                Default: 100

    --ismap_all             Switch on all alignment reporting for bwa

    --ismap_minqual INT     Mapping quality score for bwa
                                Default: 30
```


### Mapping
```
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


### Minmer Query
```
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


### Minmer Sketch
```
    --minmer_ram INT        Try to keep RAM usage below this many GB
                                Default: 5 GB

    --mash_sketch INT       Sketch size. Each sketch will have at most this
                                many non-redundant min-hashes.
                                Default: 10000

    --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                input k-mers
                                Default: 10000
```


### Quality Control 
```
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
                                Default: 20

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


### Variant Calling
```
    --snippy_ram INT        Try and keep RAM under this many GB
                                Default: 8 GB

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
