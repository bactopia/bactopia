# Bactopia Tools - *fastani*
The `fastani` tool uses [FastANI](https://github.com/ParBLiSS/FastANI) to calcualte 
the average nucleotide identity (ANI) between your samples. Although, sometimes you
might be more interested in calculating the ANI of your samples against a reference
genome. Fortunately, using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download), 
the `fastani` tool allows you  specify either a specific NCBI Assembly RefSeq 
accession (`--accession`) or a species name (`--species`) for which to download 
all RefSeq genomes for.

## Example
```
bactopia tools fastani \
    --bactopia ~/bactopia-tutorial/bactopia \
    --phyloflash ~/bactopia-tutorial/bactopia-datasets/16s/138 \
    --exclude ~/bactopia-tutorial/bactopia-tools/summary/bactopia-exclude.txt \
    --accession "GCF_900475245.1" \

awk '{if ($3 > 95){print $0}}' ~/bactopia-tutorial/bactopia-tools/fastani/fastani.tsv | \
    grep -v "GCF_900475245" > ~/bactopia-tutorial/GCF_900475245-include.txt
```

Above is a good example of subsetting your samples. In the example, all samples would
have had their ANI to GCF_900475245 calculated. Then with awk, all samples that had 
greater than 95% ANI to GCF_900475245 were output to a text file. This text file could 
then be used with the `--include` parameter for other Bactopia Tools.

## Output Overview
Below is the default output structure for the `fastani` tool. Where possible the 
file descriptions below were modified from a tools description.

```
bactopia-tools/
└── fastani/
    └── ${PREFIX}
        ├── bactopia-info
        │   ├── fastani-report.html
        │   ├── fastani-timeline.html
        │   └── fastani-trace.txt
        ├── fastani.tsv
        ├── references
        │   └── ${REFERENCE}.tsv
        └── refseq
            └── fasta
                └── ${REFERENCE}.fna
```

| Filename | Description |
|-----------|-------------|
| fastani.tsv | All the FastANI results (_references/*.tsv_) merged into a single file.  |


### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| fastani-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| fastani-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| fastani-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### references
| Filename | Description |
|----------|-------------|
| ${REFERENCE}.tsv | FastANI results of all samples against a reference genome |

#### refseq/fasta
| Filename | Description |
|----------|-------------|
| ${REFERENCE}.fna | FASTA formated genome downloaded from NCBI Assembly database. |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: fastani

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

RefSeq Assemblies Related Parameters:
    This is a completely optional step and is meant to supplement your dataset with
    high-quality completed genomes.

    --species STR           The name of the species to download RefSeq assemblies for. This
                                is a completely optional step and is meant to supplement
                                your dataset with high-quality completed genomes.

    --accession STR         A NCBI Assembly database RefSeq accession to be downloaded and included
                                in the pan-genome analysis.

    --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                number of available genomes exceeds the given limit, a 
                                random subset will be selected.
                                Default: Download all available genomes

    --refseq_only           Pairwise ANI's will only be calulated against download RefSeq genomes.

FastANI Related Parameters:
    --kmer INT              kmer size <= 16
                                Default: null

    --fragLen INT           fragment length
                                Default: 3000

    --minFraction FLOAT     Minimum fraction of genome that must be shared for trusting ANI.
                                If reference and query genome size differ, smaller one among
                                the two is considered.
                                Default: 0.2

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
Useful Parameters:
    --version               Print workflow version information
    --help                  Show this message and exit
```

## Conda Environment
Below is the command used to create the Conda environment.
```
conda create -y -n bactopia-fastani -c conda-forge -c bioconda \
    fastani \
    ncbi-genome-download \
    rename 
```

## References
* __[FastANI](https://github.com/ParBLiSS/FastANI)__  
_Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S. 
[High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries.](http://dx.doi.org/10.1038/s41467-018-07641-9)
 Nat. Commun. 9, 5114 (2018)_  

* __[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)__  
_Blin, K. [ncbi-genome-download: Scripts to download genomes from the NCBI FTP 
servers](https://github.com/kblin/ncbi-genome-download)_  
