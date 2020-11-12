# Bactopia Tools - *gtdb*
The `gtdb` tool uses [GTDB-Tk's](https://github.com/Ecogenomics/GTDBTk) classify 
workflow to assign taxonomic classifications to your set of samples. This is 
done through the use of the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/). 
If you are unsure of your sequences, `gtdb` is useful tool to help determine
the taxonomy of your samples.

## Example
The following command will use `gtdb` to assign a taxonomic classification on all 
samples except those listed in the *exclude* file.
```
bactopia tools gtdb \
    --bactopia ~/bactopia-tutorial/bactopia \
    --gtdb ~/bactopia-tutorial/bactopia-datasets/gtdb/db \
    --exclude ~/bactopia-tutorial/bactopia-tools/summary/bactopia-exclude.txt \
    --cpus 4
```

## Output Overview
Below is the default output structure for the `gtdb` tool. Where possible the 
file descriptions below were modified from descriptions at the 
[GTDB-Tk Classify Workflow](https://github.com/Ecogenomics/GTDBTk#classify-workflow) page.

```
bactopia-tools/
└── gtdb
    └── ${PREFIX}
        ├── bactopia-info
        │   ├── gtdb-report.html
        │   ├── gtdb-timeline.html
        │   └── gtdb-trace.txt
        └── classify
            ├── align
            │   ├── gtdb.(ar122|bac120).filtered.tsv
            │   ├── gtdb.(ar122|bac120).msa.fasta
            │   ├── gtdb.(ar122|bac120).user_msa.fasta
            │   └── intermediate_results
            │       └── gtdb.(ar122|bac120).marker_info.tsv
            ├── classify
            │   ├── gtdb.(ar122|bac120).classify.tree
            │   ├── gtdb.(ar122|bac120).summary.tsv
            │   └── intermediate_results
            │       ├── gtdb.(ar122|bac120).classification_pplacer.tsv
            │       ├── gtdb.(ar122|bac120).red_dictionary.tsv
            │       └── pplacer
            │           ├── pplacer.(ar122|bac120).json
            │           └── pplacer.(ar122|bac120).out
            ├── gtdb.(ar122|bac120).classify.tree
            ├── gtdb.(ar122|bac120).filtered.tsv
            ├── gtdb.(ar122|bac120).markers_summary.tsv
            ├── gtdb.(ar122|bac120).msa.fasta
            ├── gtdb.(ar122|bac120).summary.tsv
            ├── gtdb.(ar122|bac120).user_msa.fasta
            ├── gtdbtk.log
            ├── gtdbtk.warnings.log
            ├── gtdb.translation_table_summary.tsv
            └── identify
                ├── gtdb.ar122.markers_summary.tsv
                ├── gtdb.bac120.markers_summary.tsv
                ├── gtdb.translation_table_summary.tsv
                └── intermediate_results
                    └── marker_genes
                        └── ${SAMPLE_NAME}
                            ├── prodigal_translation_table.tsv
                            ├── prodigal_translation_table.tsv.sha256
                            ├── ${SAMPLE_NAME}_pfam_tophit.tsv
                            ├── ${SAMPLE_NAME}_pfam_tophit.tsv.sha256
                            ├── ${SAMPLE_NAME}_pfam.tsv
                            ├── ${SAMPLE_NAME}_pfam.tsv.sha256
                            ├── ${SAMPLE_NAME}_protein.faa
                            ├── ${SAMPLE_NAME}_protein.faa.sha256
                            ├── ${SAMPLE_NAME}_protein.fna
                            ├── ${SAMPLE_NAME}_protein.fna.sha256
                            ├── ${SAMPLE_NAME}_protein.gff
                            ├── ${SAMPLE_NAME}_protein.gff.sha256
                            ├── ${SAMPLE_NAME}_tigrfam.out
                            ├── ${SAMPLE_NAME}_tigrfam.out.sha256
                            ├── ${SAMPLE_NAME}_tigrfam_tophit.tsv
                            ├── ${SAMPLE_NAME}_tigrfam_tophit.tsv.sha256
                            ├── ${SAMPLE_NAME}_tigrfam.tsv
                            └── ${SAMPLE_NAME}_tigrfam.tsv.sha256
```

| Filename | Description |
|-----------|-------------|
| gtdb.(ar122\|bac120).classify.tree | Reference tree in Newick format containing query genomes placed with pplacer |
| gtdb.(ar122\|bac120).filtered.tsv | List of genomes with an insufficient number of amino acids in MSA |
| gtdb.(ar122\|bac120).markers_summary.tsv | Markers used in generation of the concatenated MSA and the order in which they were applied |
| gtdb.(ar122\|bac120).msa.fasta | FASTA file containing MSA of submitted and reference genomes |
| gtdb.(ar122\|bac120).summary.tsv | A summary of classifications provided by GTDB-Tk, see [classification summary](https://github.com/Ecogenomics/GTDBTk#classification-summary-file) for more details |
| gtdb.(ar122\|bac120).user_msa.fasta | FASTA file containing MSA of the submitted genomes |
| gtdbtk.log | A log of the run |
| gtdbtk.warnings.log | A log of any warnings produced by the run |
| gtdb.translation_table_summary.tsv | Summary of the tranlastlation table used for each genome |

### Directory Description

#### align
| Filename | Description |
|----------|-------------|
| gtdb.(ar122\|bac120).filtered.tsv | List of genomes with an insufficient number of amino acids in MSA |
| gtdb.(ar122\|bac120).msa.fasta | FASTA file containing MSA of submitted and reference genomes |
| gtdb.(ar122\|bac120).user_msa.fasta | FASTA file containing MSA of the submitted genomes |
| gtdb.(ar122\|bac120).marker_info.tsv | Markers used in generation of the concatenated MSA and the order in which they were applied |

#### bactopia-info
| Filename | Description |
|----------|-------------|
| gtdb-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| gtdb-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| gtdb-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### classify
| Filename | Description |
|----------|-------------|
| gtdb.(ar122\|bac120).classify.tree | Reference tree in Newick format containing query genomes placed with pplacer |
| gtdb.(ar122\|bac120).summary.tsv | Classification of query genomes based on their placement in the reference tree, relative evolutionary divergence, and ANI to reference genomes |
| gtdb.(ar122\|bac120).classification_pplacer.tsv | Classification of query genomes based only on their placement in the reference tree |
| gtdb.(ar122\|bac120).red_dictionary.tsv | Median RED values for taxonomic ranks |
| pplacer.(ar122\|bac120).json | Output information generated by pplacer in JSON format |
| pplacer.(ar122\|bac120).out | Output information generated by pplacer |

#### identify
| Filename | Description |
|----------|-------------|
| gtdb.ar122.markers_summary.tsv | Summary of unique, duplicated, and missing markers within the 122 archaeal marker set for each submitted genome |
| gtdb.bac120.markers_summary.tsv | Summary of unique, duplicated, and missing markers within the 120 bacterial marker set for each submitted genome |
| gtdb.translation_table_summary.tsv | The predicted translation table used for gene calling for each genome |
| identify/intermediate_results/marker_genes/ | Contains individual genome results for gene calling using Prodigal and gene identification based on TIGRFAM and Pfam HMMs |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

    --gtdb STR              Location of a GTDB database. If a database is not found, you must
                                use '--download_gtdb'.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: gtdb

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

GTDB-Tk Related Parameters:
    --download_gtdb         Download the latest GTDB database, even it exists.

    --min_perc_aa INT       Filter genomes with an insufficient percentage of AA in the MSA
                                Default: 10

    --recalculate_red       Recalculate RED values based on the reference tree and all added
                                user genomes

    --force_gtdb            Force GTDB to continue processing if an error occurrs on a single
                                genome

    --debug                 Create intermediate files for debugging purposes

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
conda create -y -n bactopia-gtdb -c conda-forge -c bioconda gtdbtk
```

## References
* __[FastANI](https://github.com/ParBLiSS/FastANI)__  
_Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S. 
[High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries.](http://dx.doi.org/10.1038/s41467-018-07641-9)
 Nat. Commun. 9, 5114 (2018)_  

* __[FastTree 2](http://www.microbesonline.org/fasttree)__  
_Price, M. N., Dehal, P. S. & Arkin, A. P. 
[FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments.](https://dx.doi.org/10.1371%2Fjournal.pone.0009490)
 PLoS One 5, e9490 (2010)_  

* __[Genome Taxonomy Database](https://gtdb.ecogenomic.org/)__  
_Parks, D. H. et al. 
[A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life.](https://doi.org/10.1038/nbt.4229)
 Nat. Biotechnol. 36, 996–1004 (2018) _  
_Parks, D. H. et al. 
[Selection of representative genomes for 24,706 bacterial and archaeal species clusters provide a complete genome-based taxonomy.](https://doi.org/10.1101/771964)
 bioRxiv 771964 (2019)_  

* __[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)__  
_Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. 
[GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database.](https://doi.org/10.1093/bioinformatics/btz848)
 Bioinformatics (2019)_  

* __[HMMER3](http://hmmer.org/)__  
_Eddy, S. R. 
[Accelerated Profile HMM Searches.](https://doi.org/10.1371/journal.pcbi.1002195) 
PLoS Comput. Biol. 7, e1002195 (2011)_  

* __[pplacer](https://github.com/matsen/pplacer)__  
_Matsen, F. A., Kodner, R. B. & Armbrust, E. V. 
[pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree.](https://doi.org/10.1186/1471-2105-11-538)
 BMC Bioinformatics 11, 538 (2010)_  

* __[Prodigal](https://github.com/hyattpd/Prodigal)__  
_Hyatt, D. et al. 
[Prodigal: prokaryotic gene recognition and translation initiation site identification.](https://doi.org/10.1186/1471-2105-11-119)
 BMC Bioinformatics 11, 119 (2010)_  
