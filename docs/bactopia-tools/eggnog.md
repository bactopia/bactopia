# Bactopia Tools - *eggnog*
The `eggnog` tool uses [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) to 
assign functional annotation to protein sequences. eggNOG-mapper uses orthologous groups
and phylogenies from the eggNOG database to more precisely functionally annotate than 
traditional homology methods.

## Example
The following command will run `eggnog` on each available sample.
```
bactopia tools eggnog \
    --bactopia ~/bactopia-tutorial/bactopia \
    --eggnog ~/bactopia-tutorial/eggnogdb
```

## Output Overview
Below is the default output structure for the `eggnog` tool. Where possible the 
file descriptions below were modified from a tools description.

```
bactopia-tools/
└── eggnog/
    └── ${PREFIX}
        └── bactopia-info
        │   ├── eggnog-report.html
        │   ├── eggnog-timeline.html
        │   └── eggnog-trace.txt
        ├── ${SAMPLE_NAME}.emapper.annotations
        └── ${SAMPLE_NAME}.emapper.seed_orthologs
```

| Filename | Description |
|-----------|-------------|
| ${SAMPLE_NAME}.emapper.annotations | The final [eggNOG functional annotations](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2#v200) |
| ${SAMPLE_NAME}.emapper.seed_orthologs | A list of best match for each query against the whole eggNOG protein space |

### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| eggnog-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| eggnog-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| eggnog-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

    --eggnog STR            Directory containing the the eggNOG database files:
                                eggnog.db and eggnog_proteins.dmnd.  If the database is not
                                found, you must use '--download_eggnog'.
                                WARNING: eggNOG databases stored on NFS will see a significant
                                         increase in runtimes. If possible, SSD or Ramdisk 
                                         should be used.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: emapper

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

eggNOG-mapper Parameters:
*Note: Unless specified the eggNOG-mapper defaults are used.*

    --download_eggnog       Download the latest eggNOG database, even if it exists.

eggNOG Annotation Parameters:
    --tax_scope STR         Fix the taxonomic scope used for annotation, so only orthologs
                                from a particular clade are used for functional transfer.

    --target_orthologs STR  Defines what type of orthologs should be used for functional
                                transfer.
                                Choices are: one2one, many2one, one2many, many2many, all

    --go_evidence STR       Defines what type of GO terms should be used for annotation. Choices are:
                                'experimental': Use only terms inferred from experimental evidence
                                'non-electronic': Use only non- electronically curated terms

eggNOG HMM Search Parameters:
    --hmm_maxhits INT       Max number of hits to report.

    --hmm_evalue FLOAT      E-value threshold.

    --hmm_score INT         Bit score threshold.

    --hmm_maxseqlen INT     Ignore query sequences larger than `maxseqlen`.

    --hmm_qcov FLOAT        Min query coverage (from 0 to 1).

    --Z INT                 Fixed database size used in phmmer/hmmscan (allows comparing e-values
                                among databases).

eggNOG DIAMOND Search Parameters:
    --use_diamond           Use DIAMOND instead of HMMER.

    --matrix STR            Scoring matrix. Choices are: BLOSUM62, BLOSUM90, BLOSUM80, BLOSUM50,
                                                         BLOSUM45, PAM250, PAM70, PAM30

    --gapopen INT           Gap open penalty

    --gapextend INT         Gap extend penalty

    --query_cover FLOAT     Report only alignments above the given percentage of query cover.

    --subject_cover FLOAT   Report only alignments above the given percentage of subject cover.

eggNOG Seed Ortholog Search Parameters:
    --seed_ortholog_evalue FLOAT
                            Min E-value expected when searching for seed eggNOG ortholog.
                                Applies to phmmer/diamond searches. Queries not having a
                                --significant seed orthologs will not be annotated.

    --seed_ortholog_score INT
                            Min bit score expected when searching for seed eggNOG ortholog.
                                Applies to phmmer/diamond searches. Queries not having a
                                --significant seed orthologs will not be annotated.

eggNOG Output Parameters:
    --keep_mapping_files    Do not delete temporary mapping files used for annotation (i.e.
                                HMMER and DIAMOND search outputs)

    --no_annot              Skip functional annotation, reporting only hits

    --no_file_comments      No header lines nor stats are included in the output files

    --no_refine             Skip hit refinement, reporting only HMM hits.

    --no_search             Skip HMM search mapping. Use existing hits file

eggNOG Predict Orthologs Parameters:
    --target_taxa STR       Taxa that will be searched for orthologs

    --predict_output_format STR
                            Choose the output format among: per_query, per_species.

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
conda create -n bactopia-eggnog -c conda-forge -c bioconda \
    eggnog-mapper
```

## References
* __[DIAMOND]()__  
_B. Buchfink, C. Xie, D. H. Huson, [Fast and sensitive protein alignment using DIAMOND.](http://dx.doi.org/10.1038/nmeth.3176) Nat. Methods. 12, 59–60 (2015)._  

* __[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)__  
_J. Huerta-Cepas, K. Forslund, L. P. Coelho, D. Szklarczyk, L. J. Jensen, C. von Mering, P. Bork, [Fast Genome-Wide Functional Annotation through Orthology Assignment by eggNOG-Mapper.](http://dx.doi.org/10.1093/molbev/msx148) Mol. Biol. Evol. 34, 2115–2122 (2017)._  

* __[eggNOG 5.0 Database](http://eggnog.embl.de/)__  
_J. Huerta-Cepas, D. Szklarczyk, D. Heller, A. Hernández-Plaza, S. K. Forslund, H. Cook, D. R. Mende, I. Letunic, T. Rattei, L. J. Jensen, C. von Mering, P. Bork, [eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses.](https://doi.org/10.1093/nar/gky1085) Nucleic Acids Res. 47, D309–D314 (2019)._  
