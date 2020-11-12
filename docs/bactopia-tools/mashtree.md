# Bactopia Tools - *mashtree*
The `mashtree` tool allows you to create a tree of your samples using 
[Mashtree](https://github.com/lskatz/mashtree). 

Often times, you may also want to see how your samples compare to completed genomes.
This is possible with the `mashtree` tool. If you use the `--species` parameter,
all completed genomes available from RefSeq will be downloaded with 
[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) and 
included in your tree.

## Example
The following command will run Mashtree on a set of Bactopia samples, as well as 
all completed *Bacillus cereus* genomes from RefSeq.
```
bactopia tools mashree \
    --bactopia ~/bactopia-tutorial/bactopia \
    --species "Bacillus cereus" \
    --cpus 4
```

## Output Overview
Below is the default output structure for the `mashtree` tool. Where possible the 
file descriptions below were modified from a tools description.
```
bactopia-tools/
└── mashtree
    └── ${PREFIX}
        ├── bactopia-info
        │   ├── mashtree-report.html
        │   ├── mashtree-timeline.html
        │   └── mashtree-trace.txt
        ├── refseq
        │   └── fasta
        │       └── ${SAMPLE_NAME}.fna
        ├── ${PREFIX}-mashtree.dnd
        └── ${PREFIX}-matrix.txt
```

| Filename | Description |
|-----------|-------------|
| ${PREFIX}-mashtree.dnd| A newick formatted tree based on Mash distances |
| ${PREFIX}-matrix.txt | Pair-wise Mash distance for each sample |

### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| mashtree-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| mashtree-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| mashtree-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### refseq
| Extension | Description |
|----------|-------------|
| .fna | FASTA formated genome downloaded from NCBI Assembly database. |


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
                                Default: mashtree

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

    --species STR           The name of the species to download RefSeq assemblies for.

    --accession STR         The Assembly accession (e.g. GCF*.*) download from RefSeq.

    --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                number of available genomes exceeds the given limit, a
                                random subset will be selected.
                                Default: Download all available genomes

User Procided Reference:
    --reference STR         A reference genome to calculate

Mashtree Related Parameters
    --trunclength INT       How many characters to keep in a filename
                                Default: 250

    --sortorder STR         For neighbor-joining, the sort order can make a difference.
                                Options include:  ABC (alphabetical), random, input-order
                                Default: ABC

    --genomesize INT        Genome size of the input samples.
                                Default: 5000000

    --mindepth INT          If mindepth is zero, then it will be chosen in a smart but slower
                                method, to discard lower-abundance kmers.
                                Default: 5

    --kmerlength INT        Hashes will be based on strings of this many nucleotides.
                                Default: 21

    --sketchsize INT        Each sketch will have at most this
                                many non-redundant min-hashes.
                                Default: 10000

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
Below is the command that was used to create the Conda environment.
```
conda create -n bactopia-mashtree -c conda-forge -c bioconda \
    mashtree \
    ncbi-genome-download \
    rename
```

## References

* __[Mash](https://github.com/marbl/Mash)__  
_Ondov, B. D. et al. [Mash: fast genome and metagenome distance 
estimation using MinHash](http://dx.doi.org/10.1186/s13059-016-0997-x). 
Genome Biol. 17, 132 (2016)._  

* __[Mashtree](https://github.com/lskatz/mashtree)__  
_Katz, L. S., Griswold, T., Morrison, S., Caravas, J., Zhang, S., den Bakker, H.C., Deng, X., and Carleton, H. A. [Mashtree: a rapid comparison of whole genome sequence files.](https://doi.org/10.21105/joss.01762) Journal of Open Source Software, 4(44), 1762, (2019)_  

* __[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)__  
_Blin, K. [ncbi-genome-download: Scripts to download genomes from the NCBI FTP 
servers](https://github.com/kblin/ncbi-genome-download)_  
