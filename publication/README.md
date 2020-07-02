
Below is a description of the files in this directory and subdirectories.

# `data` Folder
```
└── data
    ├── bactopia-analysis.html
    ├── fastani
    │   └── crispatus-include.txt
    ├── gtdb
    │   ├── exclude.txt
    │   ├── gtdbtk.filtered.tsv
    │   └── gtdbtk.summary.tsv
    ├── lactobacillus-accessions.txt
    ├── lactobacillus-results.txt
    ├── lactobacillus-summary.txt
    ├── phyloflash
    │   ├── phyloflash-alignment.fasta.gz
    │   ├── phyloflash-contree.txt
    │   ├── phyloflash-iqtree.txt
    │   ├── phyloflash-merged.fasta.gz
    │   └── phyloflash-summary.txt
    ├── roary
    │   ├── core-genome.aligned.fa.gz
    │   ├── core-genome.contree
    │   ├── core-genome.distance.txt
    │   └── core-genome.iqtree
    └── summary
        ├── amrfinder
        │   ├── amrfinder-gene-detailed-summary.txt
        │   ├── amrfinder-gene-summary.txt
        │   ├── amrfinder-protein-detailed-summary.txt
        │   └── amrfinder-protein-summary.txt
        ├── ariba
        │   ├── ariba-card-detailed-summary.txt
        │   ├── ariba-card-summary.txt
        │   ├── ariba-vfdb_core-detailed-summary.txt
        │   └── ariba-vfdb_core-summary.txt
        ├── lactobacillus-exclude.txt
        ├── lactobacillus-report.txt
        └── lactobacillus-summary.txt
```
This directory contains the files used to create the final results and phylogenies.

| Filename | Description |
|----------|-------------|
| bactopia-analysis.html | HTML output created from the R Markdown script bactopia-analysis.Rmd in the scripts directory |

### Directories
#### `fastani`

| Filename | Description |
|----------|-------------|
| crispatus-include.txt| List of genomes with > 95% ANI to *Lactobacillus crispatus* |

#### `gtdb`

| Filename | Description |
|----------|-------------|
| exclude.txt | Genomes not classified as Lactobacillus |
| gtdbtk.filtered.tsv| List of genomes with an insufficient number of amino acids in MSA |
| gtdbtk.summary.tsv| A summary of classifications provided by GTDB-Tk, see classification summary for more details |

#### `phyloflash`

| Filename | Description |
|----------|-------------|
| phyloflash-alignment.fasta.gz | The multiple sequence alignment of 16S genes |
| phyloflash-contree.txt| Consensus tree with assigned branch supports created from 16S alignments  |
| phyloflash-iqtree.txt| Full result of the run, this is the main report file |
| phyloflash-merged.fasta.gz| All 16S genes used in the multiple sequence alignment |
| phyloflash-summary.txt| The aggregated phyloFlash results of all samples |


#### `roary`
The results of the *Lactobacillus crispatus* pan-genome

| Filename | Description |
|----------|-------------|
| core-genome.aligned.fa.gz| The multiple sequence alignment of core genes |
| core-genome.contree | Consensus tree with assigned branch supports created from 16S alignments  |
| core-genome.distance.txt | Pairwise core genome SNP distance matrix |
| core-genome.iqtree | Full result of the IQTree run, this is the main report file |

#### `summary`

| Filename | Description |
|----------|-------------|
| {amrfinder\|ariba}-{gene\|protein\|card\|vfdb}-detailed-summary.txt | Detailed information about each hit against a specific antimicrobial resistance or Ariba dataset |
| {amrfinder\|ariba}-{gene\|protein\|card\|vfdb}-summary.txt | A presence/absence matrix for hits against a specific antimicrobial resistance or Ariba dataset  |
| lactobacillus-exclude.txt | A list of samples and the reason they failed quality cutoffs |
| lactobacillus-report.txt| A tab-delimited file containing sequence, assembly and annotation stats for all samples|
| lactobacillus-summary.txt| Brief breakdown of ranks and qc-failures |


# `figures\files\tables` Folders
```

└── figures
    ├── figure-1a-bactopia-overview.png
    ├── figure-1b-bactopia-workflow.pdf
    ├── figure-1b-bactopia-workflow.png
    ├── figure-1b-bactopia-workflow.svg
    ├── figure-2a-lactobacillus-16s.png
    ├── figure-2a-lactobacillus-16s.svg
    ├── figure-2b-lactobacillus-only-16s-annotated.png
    ├── figure-2b-lactobacillus-only-16s-annotated.svg
    ├── figure-2b-lactobacillus-only-16s.svg
    ├── figure-3-lcrispatus-core-genome-annotated.png
    ├── figure-3-lcrispatus-core-genome-annotated.svg
    ├── figure-3-lcrispatus-core-genome.svg
    ├── supplementary-figure-1-bactopia-workflow.pdf
    ├── supplementary-figure-1-bactopia-workflow.png
    ├── supplementary-figure-1-bactopia-workflow.svg
    ├── supplementary-figure-2-quality-by-year.pdf
    ├── supplementary-figure-2-quality-by-year.png
    ├── supplementary-figure-3-genome-size-assembly-vs-estimate.pdf
    ├── supplementary-figure-3-genome-size-assembly-vs-estimate.png
    ├── supplementary-figure-4-consistent-genome-size.pdf
    └── supplementary-figure-4-consistent-genome-size.png
└──  files
    ├── supplementary-data-1-lactobacillus-query-results.txt
    ├── supplementary-data-2-illumina-accessions.txt
    └── supplementary-data-3-nextflow-report.html
└── tables
    ├── supplementary-table-1-samples-excluded.txt
    ├── supplementary-table-2-non-lactobacillus-by-gtdb.txt
    ├── table-1-list-of-bioinformatic-tools.txt
    ├── table-2-comparison-of-workflows.txt
    ├── table-3-lactobacillus-sequence-summary.txt
    └── table-4-lactobacillus-crispatus-metadata.txt
```

| Filename / Directory | Description |
|----------|-------------|
| figures | Figures used in preprint |
| files | Supplementary data in the preprint |
| tables | Tab-delimited representations of tables in the preprint |

# `scripts` Folder
```
└── scripts
    ├── bactopia-analysis.Rmd
    ├── bactopia-workflow-key.R
    ├── bactopia-workflow.R
    └── lactobacillus-analysis.sh
```

This directory contains the scripts used in this analysis.

| Filename | Description |
|----------|-------------|
| bactopia-analysis.Rmd | This is the primary script for analysis results, used to create bactopia-analysis.html |
| bactopia-workflow-key.R | Used to create the key for the Bactopia Workflow diagram |
| bactopia-workflow.R | Used to create the Bactopia Workflow diagram |
| lactobacillus-analysis.sh | Commands used to run Bactopia |
