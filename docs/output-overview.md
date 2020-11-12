# Overview of Bactopia Output
After a successful run, Bactopia will have produced numerous output files. Just how many output files depends on the input datasets used (e.g. none, general datasets, species specific datasets, user populated datasets).

Here is the complete directory structure that is possible (using all available dataset options) with Bactopia. 

```
${SAMPLE_NAME}/
├── annotation
├── antimicrobial_resistance
├── ariba
├── assembly
├── blast
├── kmers
├── logs
├── mapping
├── minmers
├── mlst
├── quality-control
├── variants
└── ${SAMPLE_NAME}-genome-size.txt
```

For each type of analysis in Bactopia, a separate directory is created to hold the results. All samples processed by Bactopia will have this directory structure. The only difference is the usage of `${SAMPLE_NAME}` as a prefix for naming some output files.

## Directories
The following sections include a list of expected outputs as well as a brief description of each output file.

There are instances where additional files (e.g. `--keep_all_files` and `--ariba_noclean`) may be encountered. These files aren't described below, just the defaults. Also, using `--compress` will add a *gz* extension, but the original extension is maintained and its description still applies. 


!!! info "Developer Descriptions Take Priority"
    If a developer described their software's outputs, their description was used with a link back to the software's documentation (major thanks for taking the time to do that!). In some cases there may have been slight formatting modifications made. In any case, if descriptions are not original credit will be properly given to the source.

### `annotation`
The `annotation` directory will contain the outputs from [Prokka](https://github.com/tseemann/prokka) annotation. These outputs include FASTA (proteins and genes), GFF3, GenBank, and many more. By default the included Prokka databases are used for annotation. However, if a [Species Specific Dataset](/datasets/#species-specific) was a created the RefSeq clustered proteins are used first for annotation.

File descriptions were directly taken from [Prokka's Output Files](https://github.com/tseemann/prokka#output-files) section and slight modifications were made to the order of rows.

```
${SAMPLE_NAME}/
└── annotation
    ├── ${SAMPLE_NAME}.err
    ├── ${SAMPLE_NAME}.faa
    ├── ${SAMPLE_NAME}.ffn
    ├── ${SAMPLE_NAME}.fna
    ├── ${SAMPLE_NAME}.fsa
    ├── ${SAMPLE_NAME}.gbk
    ├── ${SAMPLE_NAME}.gff
    ├── ${SAMPLE_NAME}.log
    ├── ${SAMPLE_NAME}.sqn
    ├── ${SAMPLE_NAME}.tbl
    ├── ${SAMPLE_NAME}.tsv
    └── ${SAMPLE_NAME}.txt
```

| Extension | Description |
|-----------|-------------|
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. |
| .gbk | This is a standard GenBank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-GenBank, with one record for each sequence. |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .log | Contains all the output that Prokka produced during its run. This is a record of what settings you used. |
| .sqn | An ASN1 format "Sequin" file for submission to GenBank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .tbl | Feature Table file, used by "tbl2asn" to create the .sqn file. |
| .tsv | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product |
| .txt | Statistics relating to the annotated features found. |

### `antimicrobial_resistance`
The `antimicrobial_resistance` directory will contain the output from [NCBI's AMRFinderPlus](https://github.com/ncbi/amr). The results of AMRFinderPlus using genes as and input, and proteins as an input are available. More information about the output format is available from the [AMRFinderPlus Wiki](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format).

```
${SAMPLE_NAME}/
└── antimicrobial_resistance/
    ├── ${SAMPLE_NAME}-gene-report.txt
    └── ${SAMPLE_NAME}-protein-report.txt
```

| Extension | Description |
|----------|-------------|
| -gene-report.txt | Results of using gene sequences as an input |
| -protein-report.txt | Results of using protein sequences as an input |

### `ariba`
The `ariba` directory will contain the results of any [Ariba](https://github.com/sanger-pathogens/ariba/) analysis (excluding MLST). Only the [Ariba databases](/datasets/#general) created during the dataset setup are used for analysis. For each Ariba database (e.g. `card` or `vfdb`), a separate folder with the name of the database is included in the `ariba` folder.

The file descriptions below were modified from Ariba's wiki entries for [`run`](https://github.com/sanger-pathogens/ariba/wiki/Task:-run) and [`summary`](https://github.com/sanger-pathogens/ariba/wiki/Task:-summary).

```
${SAMPLE_NAME}/
└── ariba
    └── ARIBA_DATABASE_NAME
        ├── assembled_genes.fa.gz
        ├── assembled_seqs.fa.gz
        ├── assemblies.fa.gz
        ├── debug.report.tsv
        ├── log.clusters.gz
        ├── report.tsv
        ├── summary.csv
        └── version_info.txt
```

| Filename | Description |
|----------|-------------|
| assembled_genes.fa.gz | A gzipped FASTA file of only assembled gene sequences (with extensions). |
| assembled_seqs.fa.gz | A gzipped FASTA of the assembled sequences (genes and non-coding). |
| assemblies.fa.gz | A gzipped FASTA file of the assemblies (complete, unedited, contigs). |
| debug.report.tsv | The complete list of clusters, including those that did not pass filtering. |
| log.clusters.gz | Detailed logging for the progress of each cluster. |
| report.tsv | A detailed [report file](https://github.com/sanger-pathogens/ariba/wiki/Task:-run#report-file) of clusters which passed filtering. |
| summary.csv | A more condensed [summary](https://github.com/sanger-pathogens/ariba/wiki/Task:-summary) of the report.tsv | 
| version_info.txt | Information on the versions of ARIBA and its dependencies at runtime. |

### `assembly`
The `assembly` folder contains the results of the sample's assembly. 

#### standard
The standard assembly is managed by [Shovill](https://github.com/tseemann/shovill) and by default [SKESA](https://github.com/ncbi/SKESA) is used for assembly. Alternative assemblers include [SPAdes](https://github.com/ablab/spades), [MEGAHIT](https://github.com/voutcn/megahit), and [Velvet](https://github.com/dzerbino/velvet). Depending on the choice of assembler, additional output files (e.g. assembly graphs) may be given.

Files descriptions with some modifications were directly taken from [Shovill's Output Files](https://github.com/tseemann/shovill#output-files) section as well as the [FLASH usage](https://sourceforge.net/p/flashpage/code/ci/master/tree/flash.c#l114).
```
${SAMPLE_NAME}/
└── assembly
    ├── cointigs.fa
    ├── flash.hist
    ├── flash.histogram
    ├── shovill.corrections
    ├── shovill.log
    ├── ${SAMPLE_NAME}.fna
    └── ${SAMPLE_NAME}.fna.json
```

| Filename | Description |
|----------|-------------|
| contigs.fa | Final assembly without renamed headers. |
| flash.hist | Numeric histogram of merged read lengths. |
| flash.histogram | Visual histogram of merged read lengths |
| shovill.log | Full log file for bug reporting |
| shovill.corrections | List of post-assembly corrections |
| ${SAMPLE_NAME}.fna | The final assembly, with renamed header to include sample name |
| ${SAMPLE_NAME}.fna.json | Summary statistics of the assembly |

!!! info "FASTA inputs are not reassembled by default"
    In the case where an assembly is given as an input, the only files that will be available are `${SAMPLE_NAME}.fna` (the original unmodified assembly) and `${SAMPLE_NAME}.fna.json`. If `--reassemble` is also given, then all the files seen above will be available.

#### hybrid
If long reads are available to supplement input paired-end Illumina reads, a hybrid assembly can be created using [Unicycler](https://github.com/rrwick/Unicycler).

Files descriptions with some modifications were directly taken from [Unicycler's Output Files](https://github.com/rrwick/Unicycler#output-files).
```
${SAMPLE_NAME}/
└── assembly
    ├── 001_best_spades_graph.gfa
    ├── 002_overlaps_removed.gfa
    ├── 003_long_read_assembly.gfa
    ├── 004_bridges_applied.gfa
    ├── 005_final_clean.gfa
    ├── 006_polished.gfa
    ├── 007_rotated.gfa
    ├── assembly.fasta
    ├── assembly.gfa
    ├── ${SAMPLE_NAME}.fna
    ├── ${SAMPLE_NAME}.fna.json
    └── unicycler.log
```

| Filename | Description |
|----------|-------------|
| 001_best_spades_graph.gfa | The best SPAdes short-read assembly graph, with a bit of graph clean-up |
| 002_overlaps_removed.gfa | Overlap-free version of the SPAdes graph, with some more graph clean-up |
| 003_long_read_assembly.gfa | Assembly graph after long read assembly |
| 004_bridges_applied.gfa | Bridges applied, before any cleaning or merging |
| 005_final_clean.gfa | Assembly graph after more redundant contigs removed |
| 006_polished.gfa | Assembly graph after a round of Pilon polishing |
| 007_rotated.gfa | Assembly graph after ircular replicons rotated and/or flipped to a start position |
| assembly.fasta | The final assembly in FASTA format (same contigs names as in assembly.gfa) |
| assembly.gfa | The final assembly in [GFA v1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) graph format |
| ${SAMPLE_NAME}.fna | The final assembly with renamed header to include sample name |
| ${SAMPLE_NAME}.fna.json | Summary statistics of the assembly |
| unicycler.log | Unicycler log file (same info as stdout) |

#### quality reports
Each assembly will have its biological and technical quality assessed with [CheckM](https://github.com/Ecogenomics/CheckM) and [QUAST](http://quast.sourceforge.net/). This assessment is done no matter the input type (paired, single, hybrid, or assembly).

Files descriptions with some modifications were directly taken from [CheckM's Usage](https://github.com/Ecogenomics/CheckM/blob/master/bin/checkm)  and [QUAST's Output Files](http://quast.sourceforge.net/docs/manual.html#sec3).
```
assembly/
├── checkm
│   ├── bins/
│   ├── checkm-genes.aln
│   ├── checkm.log
│   ├── checkm-results.txt
│   ├── lineage.ms
│   └── storage/
└── quast
    ├── basic_stats/
    ├── icarus.html
    ├── icarus_viewers
    │   └── contig_size_viewer.html
    ├── predicted_genes
    │   ├── GCF_003431365_glimmer_genes.gff.gz
    │   └── GCF_003431365_glimmer.stderr
    ├── quast.log
    ├── report.{html|pdf|tex|tsv|txt}
    ├── transposed_report.tex
    ├── transposed_report.tsv
    └── transposed_report.txt
```

_CheckM Outputs_  

| Filename | Description |
|----------|-------------|
| bins/ | A folder with inputs (e.g. proteins) for processing by CheckM  |
| checkm-genes.aln | Alignment of multi-copy genes and their AAI identity |
| checkm.log | The output log of CheckM |
| checkm-results.txt | Final results of [CheckM's lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) |
| lineage.ms | Output file describing marker set for each bin |
| storage/ | A folder with intermediate results from CheckM processing |


_QUAST Outputs_  

| Filename | Description |
|----------|-------------|
| basic_stats | A folder with plots of assembly metrics (e.g. GC content, NGx, Nx) |
| icarus.html | Icarus main menu with links to interactive viewers. |
| icarus_viewers/ | Additional reports for Icarus  |
| predicted_genes/ | Predicted gene output |
| quast.log | Detailed information about the QUAST run |
| report.{html\|pdf }| Assessement summary including all tables and plots |
| report.{tex\|tsv\|txt} | Assessment summary in multiple different formats |
| transposed_report.{tex\|tsv\|txt} | Transposed version of the assessment summary |


### `blast`
The `blast` directory contains the [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) results and a BLAST database of the sample assembly.

Each of the [User Populated BLAST Sequences](/datasets/#blast) (gene, primer, or protein) is BLASTed against the sample assembly. Also if [setup](/datasets/#general), annotated genes are BLASTed against the [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/) BLAST database.

By default, results are returned in tabular format.

```
${SAMPLE_NAME}/
└── blast
    ├── blastdb
    │   ├── ${SAMPLE_NAME}.nhr
    │   ├── ${SAMPLE_NAME}.nin
    │   └── ${SAMPLE_NAME}.nsq
    ├── genes
    │   └── ${INPUT_NAME}.txt
    ├── primers
    │   └── ${INPUT_NAME}.txt
    ├── proteins
    │   └── ${INPUT_NAME}.txt
    └── ${SAMPLE_NAME}-plsdb.txt
```

| Extension | Description |
|----------|-------------|
| .nhr | Sample assembly BLAST database header file |
| .nin | Sample assembly BLAST database index file |
| .nsq | Sample assembly BLAST database sequence file |
| -plsdb.txt| The BLAST results against the PLSDB BALST database assembly. |
| .txt | The BLAST results of user input sequence(s) against the sample assembly. |

### `genome-size`
For every sample `${SAMPLE_NAME}-genome-size.txt` file is created. This file contains the genome size that was used for analysis. Genome size is used throughout Bactopia for various tasks including error correction, subsampling, and assembly.

By default, the genome size is estimated with Mash, but users have the option to provide their own value or completely disable genome size related features. Learn more about changing the genome size at [Basic Usage - Genome Size](/usage-basic/#-genome_size)

### `kmers`
The `kmers` directory contains [McCortex](https://github.com/mcveanlab/mccortex) 31-mer counts of the cleaned up FASTQ sequences. 

```
${SAMPLE_NAME}/
└── kmers
    └── ${SAMPLE_NAME}.ctx
```

| Extension | Description |
|-----------|-------------|
| .ctx      | A Cortex graph file of the 31-mer counts |

### `logs`
The `logs` folder will contain useful files for debugging or reviewing what was executed. For each process (e.g. annotation or assembly) the STDOUT and STDERR is log, as well as the time of execution and program versions. These outputs are completely optional, and can be disabled using `--skip_logs` at runtime.

```
${SAMPLE_NAME}/
└── logs/
    ├── ${PROCESS_NAME}
    │   ├── ${PROCESS_NAME}.err
    │   ├── ${PROCESS_NAME}.out
    │   ├── ${PROCESS_NAME}.sh
    │   ├── ${PROCESS_NAME}.trace
    │   ├── ${PROCESS_NAME}.versions
    │   ├── ${PROGRAM}.err
    │   └── ${PROGRAM}.out
    └── bactopia.versions

```

| Filename | Description |
|-----------|-------------|
| ${PROCESS_NAME}.err      | Any STDERR captured by the process. |
| ${PROCESS_NAME}.out      | Any STDERR captured by the process.|
| ${PROCESS_NAME}.sh       | The shell script that process used. |
| ${PROCESS_NAME}.trace    | Compute resource usage by the process (this will not always be available) |
| ${PROCESS_NAME}.versions | Date and program versions used by the process |
| ${PROGRAM}.err           | STDERR captured for a specific program. |
| ${PROGRAM}.err           | STDOUT captured for a specific program. |
| bactopia.versions        | Date and Bactopia/Nextflow versions used |

#### Example `versions`
Here is an example of the `bactopia.versions` file. 
```
# Timestamp
2020-11-11T11:12:31-05:00
# Bactopia Version
bactopia 1.4.11
# Nextflow Version
nextflow 20.07.1
```

All the `.versions` files will follow this format. The first line is always `# Timestamp` followed by the output of `date`. Then each line beginning with `#` will represent a new program and its version. 


### `mapping`
The `mapping-sequences` directory contains [BWA](https://github.com/lh3/bwa/) (bwa-mem) mapping results for each of the [User Populated Mapping Sequences](/datasets/#mapping). 

```
${SAMPLE_NAME}/
└── mapping
    └── ${MAPPING_INPUT}
        ├── ${MAPPING_INPUT}.bam
        └── ${MAPPING_INPUT}.coverage.txt
```

| Extension | Description |
|-----------|-------------|
| .bam      | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. |
| .coverage.txt | The per-base coverage of the mapping results |

### `minmers`
The `minmers` directory contains [Mash](https://github.com/marbl/Mash) and [Sourmash](https://github.com/dib-lab/sourmash) sketches of the cleaned FASTQs. If setup, it also contains the results of queries against [RefSeq, GenBank and PLSDB](/datasets/#general)

```
${SAMPLE_NAME}/
└── minmers
    ├── ${SAMPLE_NAME}-genbank-k21.txt
    ├── ${SAMPLE_NAME}-genbank-k31.txt
    ├── ${SAMPLE_NAME}-genbank-k51.txt
    ├── ${SAMPLE_NAME}-k21.msh
    ├── ${SAMPLE_NAME}-k31.msh
    ├── ${SAMPLE_NAME}-plsdb-k21.txt
    ├── ${SAMPLE_NAME}-refseq-k21.txt
    └── ${SAMPLE_NAME}.sig
```

| Extension | Description |
|-----------|-------------|
| -genbank-k(21\|31\|51).txt | Sourmash LCA Gather results against [Sourmash GenBank Signature](https://sourmash.readthedocs.io/en/latest/databases.html) (k=21,31,51)  |
| -k(21\|31).msh | A Mash sketch (k=21,31) of the sample |
| -plsdb-k21.txt | Mash Screen results against [PLSDB Mash Sketch](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/) |
| -refseq-k21.txt | Mash Screen results against [Mash Refseq Sketch](https://mash.readthedocs.io/en/latest/data.html) |
| .sig | A Sourmash signature (k=21,31,51) of the sample |


### `mlst`
If a [Species Specific Dataset](/datasets/#species-specific) has been set up, the `mlst` directory will contain [Ariba](https://github.com/sanger-pathogens/ariba) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) results for a [PubMLST.org](https://pubmlst.org/) schema. For most organisms there is only one MLST schema available, and it will be labeled as `default`. In cases where a organism has multiple schemas available they will be named following pubMLST's naming.

```
${SAMPLE_NAME}/
└── mlst
    └── ${SCHEMA}
       ├── ariba
       │   ├── assembled_genes.fa.gz
       │   ├── assembled_seqs.fa.gz
       │   ├── assemblies.fa.gz
       │   ├── debug.report.tsv
       │   ├── log.clusters.gz
       │   ├── mlst_report.details.tsv
       │   ├── mlst_report.tsv
       │   ├── report.tsv
       │   └── version_info.txt
       └── blast
           └── ${SAMPLE_NAME}-blast.json
```

| Filename | Description |
|----------|-------------|
| assembled_genes.fa.gz | A gzipped FASTA file of only assembled gene sequences (with extensions). |
| assembled_seqs.fa.gz | A gzipped FASTA of the assembled sequences (genes and non-coding). |
| assemblies.fa.gz | A gzipped FASTA file of the assemblies (complete, unedited, contigs). |
| debug.report.tsv | The complete list of clusters, including those that did not pass filtering. |
| log.clusters.gz | Detailed logging for the progress of each cluster. |
| mlst_report.details.tsv | A more detailed summary of the allele calls.  |
| mlst_report.tsv | A summary of the allele calls and identified sequence type.  |
| report.tsv | A detailed [report file](https://github.com/sanger-pathogens/ariba/wiki/Task:-run#report-file) of clusters which passed filtering. |
| summary.csv | A more condensed [summary](https://github.com/sanger-pathogens/ariba/wiki/Task:-summary) of the report.tsv | 
| version_info.txt | Information on the versions of ARIBA and its dependencies at runtime. |
| -blast.json| Allele calls and identified sequence type based on BLAST |

### `quality-control`
The `quality-control` directory contains the cleaned up FASTQs ([BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) and [Lighter](https://github.com/mourisl/Lighter)) and summary statitics ([FastQC](https://github.com/s-andrews/FastQC) and [Fastq-Scan](https://github.com/rpetit3/fastq-scan)) before and after cleanup.

```
${SAMPLE_NAME}/
└── quality-control
    ├── logs
    │   ├── bbduk-adapter.log
    │   └── bbduk-phix.log
    ├── ${SAMPLE_NAME}(|_R1|_R2).fastq.gz
    └── summary-(original|final)
        ├── ${SAMPLE_NAME}(|_R1|_R2)-(original|final)_fastqc.html
        ├── ${SAMPLE_NAME}(|_R1|_R2)-(original|final)_fastqc.zip
        └── ${SAMPLE_NAME}(|_R1|_R2)-(original|final).json
```

| Extension | Description |
|-----------|-------------|
| -adapter.log | A description of how many reads were filtered during the adapter removal step |
| -phix.log | A description of how many reads were filtered during the PhiX removal step |
| .fastq.gz | The cleaned up FASTQ(s), `_R1` and `_R2` for paired-end reads, and an empty string for single-end reads. |
| _fastqc.html | The FastQC html report of the original and final FASTQ(s) |
| _fastqc.zip | The zipped FastQC full report of the original and final FASTQ(s) |
| .json | Summary statistics (e.g. qualities and read lengths) of the original and final FASTQ(s) |


### `variants`
The `variants` directory contains the results of [Snippy](https://github.com/tseemann/snippy) variant calls against one or more reference genomes. There are two subdirectories `auto` and `user`. 

The `auto` directory includes variants calls against automatically selected reference genome(s) based on nearest Mash distance to RefSeq completed genomes. This process only happens if a [Species Specific Dataset](/datasets/#species-specific) was a created. By default, only a single reference genome (nearest) is selected. This feature can be disabled (`--disable_auto_variants`) or the number of genomes changed (`--max_references INT`).

The `user` directory contains variant calls against for each of the [User Populated Reference Genomes](/datasets/#reference-genomes).

The following description of files was directly taken from [Snippy's Output Files](https://github.com/tseemann/snippy#output-files) section. Slight modifications were made to the order of rows.

```
${SAMPLE_NAME}/
└── variants
    └── (auto|user)
        └── ${REFERENCE_NAME}
            ├── ${SAMPLE_NAME}.aligned.fa
            ├── ${SAMPLE_NAME}.annotated.vcf
            ├── ${SAMPLE_NAME}.bam
            ├── ${SAMPLE_NAME}.bam.bai
            ├── ${SAMPLE_NAME}.bed
            ├── ${SAMPLE_NAME}.consensus.fa
            ├── ${SAMPLE_NAME}.consensus.subs.fa
            ├── ${SAMPLE_NAME}.consensus.subs.masked.fa
            ├── ${SAMPLE_NAME}.coverage.txt
            ├── ${SAMPLE_NAME}.csv
            ├── ${SAMPLE_NAME}.filt.vcf
            ├── ${SAMPLE_NAME}.gff
            ├── ${SAMPLE_NAME}.html
            ├── ${SAMPLE_NAME}.log
            ├── ${SAMPLE_NAME}.raw.vcf
            ├── ${SAMPLE_NAME}.subs.vcf
            ├── ${SAMPLE_NAME}.tab
            ├── ${SAMPLE_NAME}.txt
            └── ${SAMPLE_NAME}.vcf
```

| Extension | Description |
|----------|--------------|
| .aligned.fa | A version of the reference but with `-` at position with `depth=0` and `N` for `0 < depth < --mincov` (**does not have variants**) |
| .annotated.vcf | The final variant calls with additional annotations from Reference genome's GenBank file |
| .bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Includes unmapped, multimapping reads. Excludes duplicates. |
| .bam.bai | Index for the .bam file |
| .bed | The variants in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format |
| .consensus.fa | A version of the reference genome with *all* variants instantiated |
| .consensus.subs.fa | A version of the reference genome with *only substitution* variants instantiated |
| .consensus.subs.masked.fa | A version of the reference genome with *only substitution* variants instantiated and low-coverage regions masked |
| .coverage.txt | The per-base coverage of each position in the reference genome |
| .csv | A [comma-separated](http://en.wikipedia.org/wiki/Comma-separated_values) version of the .tab file |
| .filt.vcf | The filtered variant calls from Freebayes |
| .gff | The variants in [GFF3](http://www.sequenceontology.org/gff3.shtml) format |
| .html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file |
| .log | A log file with the commands run and their outputs |
| .raw.vcf | The unfiltered variant calls from Freebayes |
| .subs.vcf | *Only substitution* variants from the final annotated variants |
| .tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants |
| .txt | A summary of the Snippy run. |
| .vcf | The final annotated variants in [VCF](http://en.wikipedia.org/wiki/Variant_Call_Format) format |
