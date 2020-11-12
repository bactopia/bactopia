# Bactopia Tools - *ismapper*
The `ismapper` tool uses [ISMapper](https://github.com/jhawkey/IS_mapper) to search for insertion sites in your samples.


## Example
The following command will run `ismapper` on each available sample.
```
bactopia tools ismapper --bactopia ~/bactopia-tutorial/bactopia \
    --reference ~/bactopia-tutorial/bactopia-datasets/ismapper/reference.gbk \
    --insertions ~/bactopia-tutorial/bactopia-datasets/ismapper/insertions.fasta \
    --cpus 4
```

## Output Overview
Below is the default output structure for the `ismapper` tool. Where possible the 
file descriptions below were modified from a tools description.
```
bactopia-tools/
└── ismapper/
    └── ${PREFIX}/
        ├── bactopia-info
        │   ├── ismapper-report.html
        │   ├── ismapper-timeline.html
        │   └── ismapper-trace.txt
        └── ${SAMPLE_NAME}
        │  ├── ${INSERTION_NAME}
        │  │   ├── ${SAMPLE_NAME}_${INSERTION_NAME}_(left|right)_final.fastq
        │  │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_finalcov.bed
        │  │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_merged.sorted.bed
        │  │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}.sorted.bam
        │  │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}.sorted.bam.bai
        │  │   ├── ${SAMPLE_NAME}_(left|right)_${SAMPLE_NAME}_${CONTIG_NUMBER}_unpaired.bed
        │  │   ├── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_closest.bed
        │  │   ├── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_intersect.bed
        │  │   └── ${SAMPLE_NAME}__${SAMPLE_NAME}_${CONTIG_NUMBER}_table.txt
        │  └── ${SAMPLE_NAME}-${INSERTION_NAME}.log
        └── genbank
            └── ${REFERENCE_NAME}.gbk
```

Below is a description of ISMapper outputs.

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


### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| TOOL_NAME-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| TOOL_NAME-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| TOOL_NAME-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### genbank
| Dxtension | Description |
|----------|-------------|
| .gbk | The reference used for analysis |

## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

    --insertions STR        Multifasta file with insertion sequence(s) to be mapped to.

    *Note: --reference or --accession is required.*
    --reference STR         Reference genome for typing against in GenBank format.

    --accession STR         The Assembly accession (e.g. GC(A|F)*.*) of the reference to
                                download from RefSeq.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: TOOL_NAME

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

ISMapper Parameters:
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

Nextflow Related Parameters:
    --condadir DIR          Directory to Nextflow should use for Conda environments
                                Default: Bactopia's Nextflow directory

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

    --conatainerPath        Path to Singularity containers to be used by the 'slurm'
                                profile.
                                Default: /opt/bactopia/singularity

    --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                will wait before execution.
                                Default: 5 seconds

    --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows
                                you to create profiles specific to your environment (e.g. SGE,
                                AWS, SLURM, etc...). This config file is loaded last and will
                                overwrite existing variables if set.
                                Default: Bactopia's default configs

    -resume                 Nextflow will attempt to resume a previous run. Please notice it is
                                only a single '-'

Useful Parameters:
    --version               Print workflow version information
    --help                  Show this message and exit
```

## Conda Environment
Below is the command used to create the Conda environment.
```
conda create -n ismapper -c conda-forge -c bioconda ismapper
```

## References
* __[Bedtools](https://github.com/arq5x/bedtools2)__  
A powerful toolset for genome arithmetic.  
_Quinlan, A. R. & Hall, I. M. [BEDTools: a flexible suite of utilities for 
comparing genomic features](http://dx.doi.org/10.1093/bioinformatics/btq033). 
Bioinformatics 26, 841–842 (2010)._  

* __[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)__  
Basic Local Alignment Search Tool  
_Camacho, C. et al. [BLAST+: architecture and applications](http://dx.doi.org/10.1186/1471-2105-10-421). 
BMC Bioinformatics 10, 421 (2009)._  

* __[BWA](https://github.com/lh3/bwa/)__  
Burrow-Wheeler Aligner for short-read alignment  
_Li, H. [Aligning sequence reads, clone sequences and assembly contigs with 
BWA-MEM](http://arxiv.org/abs/1303.3997). arXiv [q-bio.GN] (2013)._  

* __[ISMapper](https://github.com/jhawkey/IS_mapper)__  
IS mapping software  
_Hawkey, J. et al. [ISMapper: identifying transposase insertion sites in 
bacterial genomes from short read sequence data](http://dx.doi.org/10.1186/s12864-015-1860-2). 
BMC Genomics 16, 667 (2015)._  

* __[Samtools](https://github.com/samtools/samtools)__  
Tools for manipulating next-generation sequencing data  
_Li, H. et al. [The Sequence Alignment/Map format and SAMtools](http://dx.doi.org/10.1093/bioinformatics/btp352). 
Bioinformatics 25, 2078–2079 (2009)._
