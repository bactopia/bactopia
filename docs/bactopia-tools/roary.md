# Bactopia Tools - *roary*
The `roary` tool allows you to create a pan-genome of your samples using 
[Roary](https://github.com/sanger-pathogens/Roary/). 

Often times, you may
also want to include completed genomes in your pan-genome analysis. This 
is possible with the `roary` tool. If you use the `--species` parameter,
all completed genomes available from RefSeq will be downloaded with 
[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) and 
reannotated with [Prokka](https://github.com/tseemann/prokka) (to make 
compatible gffs). 

You can also use the core genome alignment, produced by Roary, to create a 
core genome phylogeny using [ClonalFrameMl](https://github.com/xavierdidelot/ClonalFrameML), 
[maskrc-svg](https://github.com/kwongj/maskrc-svg), and [IQ-TREE](http://www.iqtree.org/).
The core genome pair-wise SNP distance for each sample is also calculated with
[snp-dists](https://github.com/tseemann/snp-dists).

## Example
The following command will run Roary, on a set of samples in the *include* file. Then it 
will create a phylogenetic tree based on the core-genome alignment.
```
bactopia tools roary \
    --bactopia ~/bactopia-tutorial/bactopia \
    --include ~/bactopia-tutorial/GCF_900475245-include.txt \
    --cpus 4 \
    --n
```

## Output Overview
Below is the default output structure for the `roary` tool. Where possible the 
file descriptions below were modified from a tools description.
```
bactopia-tools
└──roary/
   └── ${PREFIX}
       ├── bactopia-info
       │   ├── roary-report.html
       │   ├── roary-timeline.html
       │   └── roary-trace.txt
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
       ├── core-genome.aligned.fa.gz
       ├── core-genome.distance.txt
       ├── core-genome.iqtree
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
       ├── refseq
       │   ├── fasta
       │   │   └── GCF_900475245.fna
       │   └── gff
       │       └── GCF_900475245.gff
       └── roary
           ├── accessory_binary_genes.fa
           ├── accessory_binary_genes.fa.newick
           ├── accessory_graph.dot
           ├── accessory.header.embl
           ├── accessory.tab
           ├── blast_identity_frequency.Rtab
           ├── clustered_proteins
           ├── conserved_vs_total_genes.png
           ├── core_accessory_graph.dot
           ├── core_accessory.header.embl
           ├── core_accessory.tab
           ├── core_alignment_header.embl
           ├── core_gene_alignment.aln.gz
           ├── gene_presence_absence.csv
           ├── gene_presence_absence.Rtab
           ├── number_of_conserved_genes.Rtab
           ├── number_of_genes_in_pan_genome.Rtab
           ├── number_of_new_genes.Rtab
           ├── number_of_unique_genes.Rtab
           ├── pan_genome_reference.fa
           ├── Rplots.pdf
           ├── summary_statistics.txt
           └── unique_vs_new_genes.png
```

| Filename | Description |
|-----------|-------------|
| core-genome.aligned.fa.gz | A multiple sequence alignment FASTA of the core genome |
| core-genome.distance.txt | Core genome Pair-wise SNP distance for each sample |
| core-genome.iqtree | Full result of the IQ-TREE core genome phylogeny |

### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| roary-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| roary-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| roary-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

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

#### refseq
| Extension | Description |
|----------|-------------|
| .fna | FASTA formated genome downloaded from NCBI Assembly database. |
| .gff | GFF output from the Prokka re-annotation of the reference FASTA |

#### roary
Where possible descriptions were taken from the [Roary Documentation](http://sanger-pathogens.github.io/Roary/).

| Filename | Description |
|----------|-------------|
| accessory_binary_genes.fa | A FASTA file with binary presence and absence of accessory genes |
| accessory_binary_genes.fa.newick | A tree created using the binary presence and absence of accessory genes |
| accessory_graph.dot | A graph in DOT format of how genes are linked together at the contig level in the accessory genome |
| accessory.header.embl | Tab/EMBL formatted file of accessory genes |
| accessory.tab | Tab/EMBL formatted file of accessory genes |
| blast_identity_frequency.Rtab | Blast results for percentage idenity graph |
| clustered_proteins | Groups file where each line lists the sequences in a cluster |
| conserved_vs_total_genes.png | Plot compairing conserved genes and total genes |
| core_accessory_graph.dot | A graph in DOT format of how genes are linked together at the contig level in the pan genome |
| core_accessory.header.embl | Tab/EMBL formatted file of core genes |
| core_accessory.tab | Tab/EMBL formatted file of core genes |
| core_alignment_header.embl | Tab/EMBL formatted file of core genome alignment |
| core_gene_alignment.aln.gz | A multi-FASTA alignment of all of the core genes |
| gene_presence_absence.csv | Lists each gene and which samples it is present in |
| gene_presence_absence.Rtab | Tab delimited binary matrix with the presence and absence of each gene in each sample |
| number_of_conserved_genes.Rtab | Graphs on how the pan genome varies as genomes are added (in random orders) |
| number_of_genes_in_pan_genome.Rtab | Graphs on how the pan genome varies as genomes are added (in random orders) |
| number_of_new_genes.Rtab | Graphs on how the pan genome varies as genomes are added (in random orders) |
| number_of_unique_genes.Rtab | Graphs on how the pan genome varies as genomes are added (in random orders) |
| pan_genome_reference.fa | FASTA file which contains a single representative nucleotide sequence from each of the clusters in the pan genome (core and accessory) |
| Rplots.pdf | PDF containing each plot |
| summary_statistics.txt | Number of genes in the core and accessory |
| unique_vs_new_genes.png | Plot compairing new vs old genes |


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
                                Default: core-genome

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 2880 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 64 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

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

    --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                number of available genomes exceeds the given limit, a 
                                random subset will be selected.
                                Default: Download all available genomes
    
    --only_completed        Pan-genome will be created using only the completed RefSeq genomes.    

    --prokka_evalue STR     Similarity e-value cut-off
                                Default: 1e-09

    --prokka_coverage INT   Minimum coverage on query protein
                                 Default: 80

Roary Related Parameters:
    --o STR                 Clusters output filename
                                Default: clustered_proteins

    --n                     Execute a fast core gene alignment with MAFFT
                                Default: Use PRANK

    --i INT                 Minimum percentage identity for blastp
                                Default: 95

    --cd INT                Percentage of isolates a gene must be in to be core
                                Default: 99%

    --g INT                 Maximum number of clusters
                                Default: 50000

    --s                     Do not split paralogs
                                Default: false

    --ap                    Allow paralogs in core alignment
                                Default: false

    --iv STR                Change the MCL inflation value
                                Default: 1.5

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
                                Default: null

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
Below is the command that was used to create the Conda environment.
```
conda create -y -n bactopia-roary -c r -c conda-forge -c bioconda \
    clonalframeml \
    iqtree \
    maskrc-svg \
    ncbi-genome-download \
    pigz \
    prokka \
    r-ggplot2 \
    rename \
    roary \
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

* __[Prokka](https://github.com/tseemann/prokka)__  
_Seemann, T. [Prokka: rapid prokaryotic genome annotation](http://dx.doi.org/10.1093/bioinformatics/btu153). 
Bioinformatics 30, 2068–2069 (2014)._  

* __[Roary](https://github.com/sanger-pathogens/Roary)__  
_Page, A. J. et al. 
[Roary: rapid large-scale prokaryote pan genome analysis.](https://doi.org/10.1093/bioinformatics/btv421)
 Bioinformatics 31, 3691–3693 (2015)_  

* __[snp-dists](https://github.com/tseemann/snp-dists)__  
_Seemann, T. [snp-dists - Pairwise SNP distance matrix from a FASTA sequence alignment.](https://github.com/tseemann/snp-dists)_  
