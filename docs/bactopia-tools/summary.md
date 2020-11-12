# Bactopia Tools - *summary*
The `summary` tool allows you to quickly aggregate the results of your Bactopia
analysis. For each sample the sequence stats (before and after QC), assembly stats, 
and the annotation stats are put into a single tab-delimited file.

For each sample, the `summary` assigns a rank of *Gold*, *Silver*, *Bronze*,
or *Fail*. The rank is determined by sequence quality and assembly quality. Below
is the default cutoffs for each rank.

| Rank | Coverage | Mean Per-Read Quality | Mean Read Length | Total Contigs |
|----------|:-------------:|:---:|:---:|:--:|
| Gold | 100x | Q30 | 95bp | 100 |
| Silver | 50x | Q20 | 75bp | 200 |
| Bronze | 20x | Q12 | 49bp | 500 |
| Fail | <20x | <@12 | <49bp | >500 |

Samples that fail to meet all the cutoffs for at least a *Bronze* rank are added
to a *exclude* file. This turns out to be a useful feature beacuse all other 
Bactopia Tools can read this file and automatically
exclude the samples marked as *Fail* from downstream analysis.

## Example
```
bactopia tools summary --bactopia ~/bactopia-tutorial/bactopia
```

## Output Overview
```
bactopia-tools/
└── summary/
    ├── amrfinder
    │   ├── amrfinder-(gene|protein)-detailed-summary.txt
    │   └── amrfinder-(gene|protein)-summary.txt
    ├── ariba
    │   ├── ariba-(card|vfdb|etc...)-detailed-summary.txt
    │   └── ariba-(card|vfdb|etc...)-summary.txt
    ├── bactopia-exclude.txt
    ├── bactopia-info
    │   ├── summary-report.html
    │   ├── summary-timeline.html
    │   └── summary-trace.txt
    ├── bactopia-results.txt
    └── bactopia-summary.txt
```

| Filename | Description |
|----------|-------------|
| bactopia-exclude.txt | A list of samples and the reason they failed quality cutoffs |
| bactopia-results.txt | A tab-delimited file containing sequence, assembly and annotation stats for all samples |
| bactopia-summary.txt | Brief breakdown of ranks and qc-failures |

### Directory Description

#### amrfinder
| Filename | Description |
|----------|-------------|
| amrfinder-(gene\|protein)-detailed-summary.txt | Detailed information about each hit against a specific antimicrobial resistance |
| amrfinder-(gene\|protein)-summary.txt | A presence/absence matrix for hits against a specific antimicrobial resistance |

#### ariba
| Filename | Description |
|----------|-------------|
| ariba-(card\|vfdb\|etc...)-detailed-summary.txt | Detailed information about each hit against a reference Ariba dataset |
| ariba-(card\|vfdb\|etc...)-summary.txt | A presence/absence matrix for hits against a reference Ariba dataset |

#### bactopia-info
| Filename | Description |
|----------|-------------|
| summary-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| summary-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| summary-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

Bactopia Summary Parameters:
    --gold_coverage FLOAT   Minimum amount of coverage required for Gold status
                                Default: 100

    --gold_quality INT      Minimum per-read mean quality score required for Gold
                                status
                                Default: 30

    --gold_read_length INT  Minimum mean read length required for Gold status
                                Default: 95

    --gold_contigs INT      Maximum contig count required for Gold status
                                Default: 100

    --silver_coverage FLOAT Minimum amount of coverage required for Silver status
                                Default: 50

    --silver_quality INT    Minimum per-read mean quality score required for
                                Silver status
                                Default: 20

    --silver_read_length INT
                            Minimum mean read length required for Silver status
                                Default: 75

    --silver_contigs INT    Maximum contig count required for Silver status
                                Default: 200

    --min_coverage FLOAT    Minimum amount of coverage required to pass
                                Default: 20

    --min_quality INT       Minimum per-read mean quality score required to pass
                                Default: 12

    --min_read_length INT   Minimum mean read length required to pass
                                Default: 49

    --max_contigs INT       Maximum contig count required to pass
                                Default: 500

    --min_genome_size INT   Minimum assembled genome size.
                                Default: null

    --max_genome_size INT   Maximum assembled genome size.
                                Default: null

Ariba Summary Parameters:
    --all_hits              Include all hits (matches and partials) in the summary
                                Default: Only report hits that are a match

AMRFinder+ Summary Parameters:
    --subclass              Group the report by subclass (ex. Streptomycin).
                                Default: Group by class (ex. Aminoglycoside)

Optional Parameters:
    --prefix STR            Prefix to use for final output files
                                Default: bactopia

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 4

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
    --verbose               Increase the verbosity of processes.
    --version               Print workflow version information
    --help                  Show this message and exit
```
