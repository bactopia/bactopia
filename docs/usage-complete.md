# Runtime Parameters

Bactopia includes numerous configurable parameters. Currently 103 to be exact! Basically for each step of the pipeline, you can modify the default parameters of a specific tool.

## Required
The required parameters depends on how many samples are to be proccessed. You can learn more about which approach to take at [Specifying Input FASTQs](usage-basic.md#specifying-input-fastqs).
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
These optional parameters, while not required, will be quite useful (especially `--max_cpus` and `--cpus`!) to tweak.

```
    --coverage INT          Reduce samples to a given coverage
                                Default: 100x

    --genome_size INT       Expected genome size (bp) for all samples
                                Default: Mash Estimate

    --outdir DIR            Directory to write results to
                                Default .
```


### `--max_cpus` vs `--cpus`
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

By default `--max_cpus` is set to 1 and if `--cpus` is set to a value greater than `--max_cpus` it will be set equal to `--max_cpus`. This appoach errs on the side of caution, by not occupying all cpus on the server without the users consent!


## Helpers
The following parameters are useful to test to test input parameters.
```
    --available_datasets    Print a list of available datasets found based
                                on location given by "--datasets"

    --example_fastqs        Print example of expected input for FASTQs file

    --check_fastqs          Verify "--fastqs" produces the expected inputs

    --clean_cache           Removes 'work' and '.nextflow' logs. Caution, if used,
                                the Nextflow run cannot be resumed.

    --version               Print workflow version information
    --help                  Show this message and exit
    --help_all              Show a complete list of adjustable parameters
```

### `--clean_cache`
By default, *Bactopia* will keep all intermediate files even after successfully completing all steps. While the cache is maintained Bactopia is resumable using the `-resume` parameter. This does however introduce a potentential storage overhead. The cache will contain multiple intermediate files (e.g. uncompressed FASTQs, BAMs, etc...) for each sample that was processed. In other words, it can get pretty large!

If you would like to clean up the cache you can use `--clean_cache`. This will remove the cache **only** after a successful execution (e.g. everything thing finished without errors). This is accomplished by removing the `work` directory created by Nextflow. As you might have guessed, by using removing the cache, Bactopia will no longer be resumeable.

At the end of the day, you can always decide to not use `--clean_cache` and manually remove the `work` directory when you feel it is safe!

## Program Specific
The remaining parameters are associated with specific programs. In the following sectionsm, these parameters are grouped by which Nextflow process they are applicable to.

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
```


### Assembly
```
    --shovill_ram INT       Try to keep RAM usage below this many GB
                                Default: 32

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
    --keep_singletons       Keep all counted 31-mers
                                Default: Filter out singletons
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
    --mash_sketch INT       Sketch size. Each sketch will have at most this
                                many non-redundant min-hashes.
                                Default: 10000

    --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                input k-mers
                                Default: 10000
```


### Quality Control 
```
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
                                Default: 8

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
