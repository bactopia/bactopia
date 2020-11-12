# Bactopia Tools - *pirate*
The `pirate` tool allows you to create a pan-genome of your samples using 
[PIRATE](https://github.com/SionBayliss/PIRATE). 

Often times, you may also want to include completed genomes in your pan-genome 
analysis. This is possible with the `pirate` tool. If you use the `--species` parameter,
all completed genomes available from RefSeq will be downloaded with 
[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) and 
reannotated with [Prokka](https://github.com/tseemann/prokka) (to make 
compatible gffs). 

You can also use the core genome alignment, produced by PIRATE, to create a 
core genome phylogeny using [ClonalFrameMl](https://github.com/xavierdidelot/ClonalFrameML), 
[maskrc-svg](https://github.com/kwongj/maskrc-svg), and [IQ-TREE](http://www.iqtree.org/).
The core genome pair-wise SNP distance for each sample is also calculated with
[snp-dists](https://github.com/tseemann/snp-dists).

## Example
The following command will run PIRATE, on a set of samples in the include file. 
Then it will create a phylogenetic tree based on the core-genome alignment.
```
bactopia tools pirate \
    --bactopia ~/bactopia-tutorial/bactopia \
    --include ~/bactopia-tutorial/GCF_900475245-include.txt \
    --cpus 4
```

## Output Overview
Below is the default output structure for the `pirate` tool. Where possible the 
file descriptions below were modified from a tools description.
```
bactopia-tools/
└── pirate/
    └── ${PREFIX}
        ├── bactopia-info
        │   ├── pirate-report.html
        │   ├── pirate-timeline.html
        │   └── pirate-trace.txt
        ├── clonalframe
        │   ├── clonalframe.emsim.txt
        │   ├── clonalframe.em.txt
        │   ├── clonalframe.importation_status.txt
        │   ├── clonalframe.labelled_tree.newick
        │   ├── clonalframe.ML_sequence.fasta
        │   ├── clonalframe.position_cross_reference.txt
        │   ├── core_gene_alignment-masked.aln.gz
        │   ├── start-tree.bionj
        │   ├── start-tree.ckp.gz
        │   ├── start-tree.iqtree
        │   ├── start-tree.log
        │   ├── start-tree.mldist
        │   ├── start-tree.model.gz
        │   └── start-tree.treefile
        ├── iqtree
        │   ├── core-genome.alninfo
        │   ├── core-genome.bionj
        │   ├── core-genome.ckp.gz
        │   ├── core-genome.contree
        │   ├── core-genome.iqtree
        │   ├── core-genome.log
        │   ├── core-genome.mldist
        │   ├── core-genome.model.gz
        │   ├── core-genome.splits.nex
        │   ├── core-genome.treefile
        │   └── core-genome.ufboot
        ├── pirate
        │   ├── binary_presence_absence.{fasta|nwk}
        │   ├── cluster_alleles.tab
        │   ├── co-ords
        │   │   └── ${SAMPLE_NAME}.co-ords.tab
        │   ├── core_alignment.fasta
        │   ├── core_alignment.gff
        │   ├── feature_sequences
        │   │   └── ${GENE_FAMILY}.{aa|nucleotide|.fasta
        │   ├── genome2loci.tab
        │   ├── genome_list.txt
        │   ├── loci_list.tab
        │   ├── loci_paralog_categories.tab
        │   ├── loci_paralog_categories.tab.idx
        │   ├── modified_gffs
        │   │   └── ${SAMPLE_NAME}.gff
        │   ├── pangenome_alignment.fasta.gz
        │   ├── pangenome_alignment.gff
        │   ├── pangenome.connected_blocks.tsv
        │   ├── pangenome.edges
        │   ├── pangenome.gfa
        │   ├── pangenome_iterations
        │   │   ├── pan_sequences.{50|60|70|80|90|95|98}.reclustered.reinflated
        │   │   ├── pan_sequences.blast.output
        │   │   ├── pan_sequences.cdhit_clusters
        │   │   ├── pan_sequences.core_clusters.tab
        │   │   ├── pan_sequences.mcl_log.txt
        │   │   └── pan_sequences.representative.fasta
        │   ├── pangenome.order.tsv
        │   ├── pangenome.reversed.tsv
        │   ├── pangenome.syntenic_blocks.tsv
        │   ├── pan_sequences.fasta
        │   ├── paralog_clusters.tab
        │   ├── PIRATE.gene_families.ordered.tsv
        │   ├── PIRATE.gene_families.tsv
        │   ├── PIRATE.genomes_per_allele.tsv
        │   ├── PIRATE.log
        │   ├── PIRATE.pangenome_summary.txt
        │   ├── PIRATE_plots.pdf
        │   ├── PIRATE.unique_alleles.tsv
        │   └── split_groups.log
        ├── refseq
        │    ├── fasta
        │    │   └── *.fna
        │    └── gff
        │        └── *.gff 
        ├── ${PREFIX}.aligned.fa.gz
        ├── ${PREFIX}.distance.txt
        └── ${PREFIX}.iqtree
```

| Filename | Description |
|-----------|-------------|
| ${PREFIX}.aligned.fa.gz | A multiple sequence alignment FASTA of the core genome |
| ${PREFIX}.distance.txt | Core genome Pair-wise SNP distance for each sample |
| ${PREFIX}.iqtree | Full result of the IQ-TREE core genome phylogeny |

### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| pirate-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| pirate-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| pirate-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### clonalframe
Where possible descriptions were taken from the [ClonalFrameML Wiki](https://github.com/xavierdidelot/clonalframeml/wiki), 
IQ-TREE's [Command Reference](https://github.com/Cibiv/IQ-TREE/wiki/Command-Reference)
page, [Web Server Tutorial](https://github.com/Cibiv/IQ-TREE/wiki/Web-Server-Tutorial) page, 
and the [Tutorial](http://www.iqtree.org/doc/Tutorial) page.

| Filename | Description |
|----------|-------------|
| clonalframe.emsim.txt | The bootstrapped values for the three parameters R/theta, nu and delta |
| clonalframe.em.txt | The point estimates for R/theta, nu, delta and the branch lengths |
| clonalframe.importation_status.txt | The list of reconstructed recombination events |
| clonalframe.labelled_tree.newick | The output tree with all nodes labelled so that they can be referred to in other files |
| clonalframe.ML_sequence.fasta | The sequence reconstructed by maximum likelihood for all internal nodes of the phylogeny, as well as for all missing data in the input sequences |
| clonalframe.position_cross_reference.txt | A vector of comma-separated values indicating for each location in the input sequence file the corresponding position in the sequences in the output *ML_sequence.fasta* file |
| core_gene_alignment-masked.aln.gz | A core-genome alignment with the recomination masked |
| start-tree.bionj | A neighbor joining tree produced by BIONJ |
| start-tree.ckp.gz | IQ-TREE writes a checkpoint file |
| start-tree.iqtree | Full result of the run, this is the main report file |
| start-tree.log | Run log |
| start-tree.mldist | Contains the likelihood distances |
| start-tree.model.gz | Information about all models tested |
| start-tree.treefile | Maximum likelihood tree in NEWICK format, can be visualized with treeviewer programs |

#### iqtree
Where possible descriptions were taken from 
IQ-TREE's [Command Reference](https://github.com/Cibiv/IQ-TREE/wiki/Command-Reference)
page, [Web Server Tutorial](https://github.com/Cibiv/IQ-TREE/wiki/Web-Server-Tutorial) page, 
and the [Tutorial](http://www.iqtree.org/doc/Tutorial) page.

| Filename | Description |
|----------|-------------|
| core-genome.alninfo | Alignment site statistics |
| core-genome.bionj | A neighbor joining tree produced by BIONJ |
| core-genome.ckp.gz | IQ-TREE writes a checkpoint file |
| core-genome.contree | Consensus tree with assigned branch supports where branch lengths are optimized on the original alignment; printed if Ultrafast Bootstrap is selected |
| core-genome.iqtree | Full result of the run, this is the main report file |
| core-genome.log | Run log |
| core-genome.mldist | Contains the likelihood distances |
| core-genome.model.gz | Information about all models tested |
| core-genome.splits.nex | Support values in percentage for all splits (bipartitions), computed as the occurence frequencies in the bootstrap trees |
| core-genome.treefile | Maximum likelihood tree in NEWICK format, can be visualized with treeviewer programs |
| core-genome.ufboot | Trees created during the bootstrap steps |

#### pirate
Where possible descriptions were taken from 
PIRATE's [Output files](https://github.com/SionBayliss/PIRATE#output-files)
page.

| Filename | Description |
|----------|-------------|
| binary_presence_absence.{fasta\|nwk} | A tree (.nwk) generated by fasttree from binary gene_family presence-absence data and the fasta file used to create it |
| cluster_alleles.tab | List of alleles in paralogous clusters |
| co-ords/${SAMPLE_NAME}.co-ords.tab | Gene feature co-ordinates for each sample |
| core_alignment.fasta | Gene-by-gene nucleotide alignments of the core genome created using MAFFT |
| core_alignment.gff | Annotation containing the position of the gene family within the core genome alignment |
| feature_sequences/${GENE_FAMILY}.{aa\|nucleotide}.fasta | Amino acid and nucleotide sequences for each gene family |
| genome2loci.tab | List of loci for each genome |
| genome_list.txt | List of genomes in the analysis |
| loci_list.tab | List of loci and their associated genomes |
| loci_paralog_categories.tab | Concatenation of classified paralogs |
| loci_paralog_categories.tab.idx | Index of the classified paralogs |
| modified_gffs/${SAMPLE_NAME}.gff | GFF3 files which have been standardised for PIRATE |
| pangenome_alignment.fasta.gz | Gene-by-gene nucleotide alignments of the full pangenome created using MAFFT |
| pangenome_alignment.gff | Annotation containing the position of the gene family within the pangenome alignment |
| pangenome.connected_blocks.tsv | List of connected blocks in the pangenome graph |
| pangenome.edges | List of classified edges in the pangenome graph |
| pangenome.gfa | GFA network file representing all unique connections between gene families |
| pangenome_iterations/pan_sequences.{50\|60\|70\|80\|90\|95\|98}.reclustered.reinflated | List of clusters for each reinflation threshold  |
| pangenome_iterations/pan_sequences.blast.output | BLAST output of sequences against representatives and self hits. |
| pangenome_iterations/pan_sequences.cdhit_clusters | A list of CDHIT representative clusters |
| pangenome_iterations/pan_sequences.core_clusters.tab | A list of core clusters. |
| pangenome_iterations/pan_sequences.mcl_log.txt | A log file from `mcxdeblast` and `mcl` |
| pangenome_iterations/pan_sequences.representative.fasta | FASTA file with sequences for each representative cluster |
| pangenome.order.tsv | Sorted list gene_families file on pangenome graph |
| pangenome.reversed.tsv | List of reversed blocks in the pangenome graph |
| pangenome.syntenic_blocks.tsv | List of syntenic blocks in the pangenome graph |
| pan_sequences.fasta | All representative sequences in the pangenome  |
| paralog_clusters.tab | List of paralogous clusters |
| PIRATE.gene_families.ordered.tsv | Tabular summary of all gene families ordered on syntenic regions in the pangenome graph |
| PIRATE.gene_families.tsv | Tabular summary of all gene families |
| PIRATE.genomes_per_allele.tsv | A list of genomes associated with each allele |
| PIRATE.log | PIRATE log file |
| PIRATE.pangenome_summary.txt | Short summary of the number and frequency of genes in the pangenome |
| PIRATE_plots.pdf | Summary plots of the PIRATE pangenome |
| PIRATE.unique_alleles.tsv | Tabular summary of all unique alleles of each gene family |
| split_groups.log | Concatenation of log files from splitting paralogs |


#### refseq
| Extension | Description |
|----------|-------------|
| .fna | FASTA formated genome downloaded from NCBI Assembly database. |
| .gff | GFF output from the Prokka re-annotation of the reference FASTA |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.
                                This parameter is not required if '--only_completed' is used.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: core-genome

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 2880 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 64 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 10

RefSeq Assemblies Related Parameters:
    --assembly              A single assembly, or directory of assemblies to be included in the
                                pan-genome analysis. If compressed, gzip and the ".gz" extension
                                must be used.

    --assembly_pattern      If a directory is given, use the given pattern to match assemblies.
                                Default: *.fna

    --species STR           The name of the species to download RefSeq assemblies for. This
                                is a completely optional step and is meant to supplement
                                your dataset with high-quality completed genomes.

    --accession STR         A NCBI Assembly database RefSeq accession to be downloaded and included
                                in the pan-genome analysis.

    --accessions STR        A file with Assembly accessions (e.g. GCF*.*) to download from RefSeq.

    --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                number of available genomes exceeds the given limit, a
                                random subset will be selected.
                                Default: Download all available genomes

    --only_completed        Pan-genome will be created using only the completed RefSeq genomes. Requires
                                either '--accessions' and/or '--species'

    --prokka_evalue STR     Similarity e-value cut-off
                                Default: 1e-09

    --prokka_coverage INT   Minimum coverage on query protein
                                Default: 80

PIRATE Related Parameters:
    --steps STR                 Percent identity thresholds to use for pangenome construction
                                    Default: 50,60,70,80,90,95,98

    --features STR              Choose features to use for pangenome construction. Multiple may be
                                    entered, separated by a comma.
                                    Default: CDS

    --nucl                      Input CDS are nucleotides (e.g. not translated to AA sequence)

    --para_off                  Switch off paralog identification

    --keep_all_files            Retain all intermediate files

PIRATE Advanced Parameters:
    --perc INT                  Single percent identity threshold to use for pangenome
                                    Default: 98 %

    --cd_low INT                CDHIT lowest percentage identity threshold
                                    Default: 98 %

    --cd_step FLOAT             CDHIT step size
                                    Default: 0.5

    --evalue STR                E-value used for BLAST hit filtering
                                    Default: 1E-6

    --hsp_len FLOAT             Remove BLAST hsps that are < hsp_len proportion of query length
                                    Default: 0

    --mcl_inflation FLOAT       MCL inflation value
                                    Default: 1.5

    --use_diamond               Use Diamond instead of BLAST - incompatible with --nucl

    --split_diamond             Split diamond files into batches for processing


IQ-TREE Related Parameters:
    --skip_phylogeny        Skip the creation a core-genome based phylogeny

    --m STR                 Substitution model name
                                Default: MFP

    --bb INT                Ultrafast bootstrap replicates
                                Default: 1000

    --alrt INT              SH-like approximate likelihood ratio test replicates
                                Default: 1000

    --asr                   Ancestral state reconstruction by empirical Bayes
                                Default: false

    --iqtree_opts STR       Extra IQ-TREE options in quotes.
                                Default: ''

ClonalFrameML Related Parameters:
    --skip_clonalframe      Skip the ClonalFrameML and use the original core-genome
                                alignment for the final tree.

    --emsim INT             Number of simulations to estimate uncertainty in the EM results.
                                Default: 100

    --clonal_opts STR       Extra ClonalFrameML options in quotes.
                                Default: ''

SNP-Dists Related Parameters:
    --a                     Count all differences not just [AGTC]
                                Default: false

    --b                     Blank top left corner cell
                                Default: false

    --c                     Output CSV instead of TSV
                                Default: false

    --k                     Keep case, don't uppercase all letters
                                Default: false

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
                                Default: true

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
conda create -y -n bactopia-pirate -c r -c conda-forge -c bioconda \
    bioconductor-ggtree \
    clonalframeml \
    iqtree \
    maskrc-svg \
    ncbi-genome-download \
    pigz \
    pirate \
    prokka \
    r-dplyr \
    r-ggplot2 \
    r-gridextra \
    r-phangorn \
    rename \
    snp-dists \
    tbl2asn-forever
```

## References
* __[ClonalFramML](https://github.com/xavierdidelot/ClonalFrameML)__  
_Didelot, X. & Wilson, D. J. 
[ClonalFrameML: efficient inference of recombination in whole bacterial genomes.](http://dx.doi.org/10.1371/journal.pcbi.1004041)
 PLoS Comput. Biol. 11, e1004041 (2015)_  

* __[IQ-TREE](https://github.com/Cibiv/IQ-TREE)__  
_L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) 
[IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.](https://doi.org/10.1093/molbev/msu300)
 Mol. Biol. Evol., 32:268-274._  
_S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, L.S. Jermiin (2017) 
[ModelFinder: Fast model selection for accurate phylogenetic estimates.](https://doi.org/10.1038/nmeth.4285) 
Nat. Methods, 14:587-589._  
_D.T. Hoang, O. Chernomor, A. von Haeseler, B.Q. Minh, L.S. Vinh (2018) [
UFBoot2: Improving the ultrafast bootstrap approximation.](https://doi.org/10.1093/molbev/msx281)
 Mol. Biol. Evol., 35:518–522._  

* __[maskrc-svg](https://github.com/kwongj/maskrc-svg)__  
_Kwong, J. [maskrc-svg - Masks recombination as detected by ClonalFrameML or Gubbins and draws an SVG.](https://github.com/kwongj/maskrc-svg)_  

* __[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)__  
_Blin, K. [ncbi-genome-download: Scripts to download genomes from the NCBI FTP 
servers](https://github.com/kblin/ncbi-genome-download)_  

* __[PIRATE](http://github.com/SionBayliss/PIRATE)__  
_S. C. Bayliss, H. A. Thorpe, N. M. Coyle, S. K. Sheppard, E. J. Feil (2019)
[PIRATE: A fast and scalable pangenomics toolbox for clustering diverged orthologues in bacteria.](https://doi.org/10.1093/gigascience/giz119) 
Gigascience. 8_  

* __[Prokka](https://github.com/tseemann/prokka)__  
_Seemann, T. [Prokka: rapid prokaryotic genome annotation](http://dx.doi.org/10.1093/bioinformatics/btu153). 
Bioinformatics 30, 2068–2069 (2014)._  

* __[snp-dists](https://github.com/tseemann/snp-dists)__  
_Seemann, T. [snp-dists - Pairwise SNP distance matrix from a FASTA sequence alignment.](https://github.com/tseemann/snp-dists)_  
