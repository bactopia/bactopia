# Overview of Bactopia Output
After a successful run Bactopia will have produced a lot of output files. Just how many output files depends on the input datasets used (e.g. none, general datasets, species specific datasets, user populated datasets).

Here is the complete directory structure that is possible (using all available dataset options) with Bactopia. 

```
${SAMPLE_NAME}/
├── annotation
├── antimicrobial_resistance
├── ariba
├── assembly
├── blast
├── insertion-sequences
├── kmers
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
    If a developer described their tool's outputs, their description was used with a link back to the tool's documentation (major thanks for taking the time to do that!). In some cases there may have been slight formatting modifications made. In any case, if descriptions are not original credit will be properly given to the source.

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
The `assembly` folder contains the results of the sample's assembly. Assembly is managed by [Shovill](https://github.com/tseemann/shovill) and by default [SKESA](https://github.com/ncbi/SKESA) is used for assembly. Alternative assemblers include [SPAdes](https://github.com/ablab/spades), [MEGAHIT](https://github.com/voutcn/megahit), and [Velvet](https://github.com/dzerbino/velvet). Depending on the choice of assembler, additional output files (e.g. assembly graphs) may be given.

Files descriptions with some modifications were directly taken from [Shovill's Output Files](https://github.com/tseemann/shovill#output-files) section as well as the [FLASH usage](https://sourceforge.net/p/flashpage/code/ci/master/tree/flash.c#l114).


```
${SAMPLE_NAME}/
└── assembly
    ├── flash.hist
    ├── flash.histogram
    ├── shovill.corrections
    ├── shovill.log
    ├── ${SAMPLE_NAME}.fna
    └── ${SAMPLE_NAME}.fna.json
```

| Filename | Description |
|----------|-------------|
| flash.hist | Numeric histogram of merged read lengths. |
| flash.histogram | Visual histogram of merged read lengths. |
| shovill.log | Full log file for bug reporting |
| shovill.corrections | List of post-assembly corrections |
| ${SAMPLE_NAME}.fna | The final assembly you should use |
| ${SAMPLE_NAME}.fna.json | Summary statistics of the assembly |

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

### `insertion-sequences`
The `insertion-sequences` directory contains [ISMapper](https://github.com/jhawkey/IS_mapper) results for each of the [User Populated Insertion Sequences](/datasets/#insertion-sequences). 
```
${SAMPLE_NAME}/
└── insertion-sequences
    ├── ${INSERTION_NAME}
    │   ├── ${SAMPLE_NAME}_${INSERTION_NAME}_(left|right)_final.fastq
    │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_finalcov.bed
    │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_merged.sorted.bed
    │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}.sorted.bam
    │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}.sorted.bam.bai
    │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_unpaired.bed
    │   ├── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_closest.bed
    │   ├── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_intersect.bed
    │   ├── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_table.txt
    └── ${SAMPLE_NAME}-${INSERTION_NAME}.log
```

| Extension | Description |
|-----------|-------------|
| _final.fastq | Sequences (FASTQ format) that mapped to the flanking regions of the IS query |
| _finalcov.bed | Contains information about the coverage of the IS query |
| _merged.sorted.bed | Merged overlapping regions that passed coverage cutoffs |
| .sorted.bam | Reads mapped to the IS query. |
| .sorted.bam.bai | An index of the sorted BAM file. |
| _unpaired.bed | All unpaired mappings to the IS query  |
| _closest.bed | Merged regions that are close but do not overlap |
| _intersect.bed | An intersection of merged regions from the left and right flanks. |
| _table.txt | A [detailed description](https://github.com/jhawkey/IS_mapper#single-isolate-output) of the IS query results. |
| .log | Information logged during the execution of ISMapper |

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
If a [Species Specific Dataset](/datasets/#species-specific) has been set up, the `mlst` directory will contain [Ariba](https://github.com/sanger-pathogens/ariba) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) results for a [PubMLST.org](https://pubmlst.org/) schema.

```
${SAMPLE_NAME}/
└── mlst
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
