# Bactopia Tools - *phyloflash*
The `phyloflash` tool uses [phyloFlash](https://github.com/HRGV/phyloFlash) to 
resconstruct 16S rRNA genes from your input samples. Optionally these reconstructed
genes can then be aligned to one another with [MAFFT](https://mafft.cbrc.jp/alignment/software/) 
and a phylogenetic representation created using [IQ-TREE](http://www.iqtree.org/)

## Example
The following command will reconstruct the 16S rRNA gene for each sample except those listed in the *exclude* file.
```
bactopia tools phyloflash \
    --bactopia ~/bactopia-tutorial/bactopia \
    --phyloflash ~/bactopia-tutorial/bactopia-datasets/16s/138 \
    --exclude ~/bactopia-tutorial/bactopia-tools/summary/bactopia-exclude.txt
```

## Output Overview
Below is the default output structure for the  `phyloflash` tool. Where possible the 
file descriptions below were modified from a tools description.

```
bactopia-tools/
└── phyloflash/
    └── ${PREFIX}
        ├── alignment
        │   ├── phyloflash-alignment.fasta
        │   └── phyloflash-matches.txt
        ├── bactopia-info
        │   ├── phyloflash-report.html
        │   ├── phyloflash-timeline.html
        │   └── phyloflash-trace.txt
        ├── iqtree
        │   ├── 16s.alninfo
        │   ├── 16s.bionj
        │   ├── 16s.ckp.gz
        │   ├── 16s.iqtree
        │   ├── 16s.log
        │   ├── 16s.mldist
        │   ├── 16s.model.gz
        │   ├── 16s.treefile
        │   └── 16s.uniqueseq.phy
        ├── phyloflash.iqtree
        ├── phyloflash-summary.txt
        └── samples
            └── ${SAMPLE_NAME}
                ├── ${SAMPLE_NAME}.all.dbhits.NR97.fa
                ├── ${SAMPLE_NAME}.all.final.fasta
                ├── ${SAMPLE_NAME}.all.final.phyloFlash.dbhits.fa
                ├── ${SAMPLE_NAME}.all.final.phyloFlash.notmatched.fa
                ├── ${SAMPLE_NAME}.all.vsearch.csv
                ├── ${SAMPLE_NAME}.assemratio.csv
                ├── ${SAMPLE_NAME}.assemratio.csv.svg
                ├── ${SAMPLE_NAME}.bbmap.out
                ├── ${SAMPLE_NAME}.bbmap.sam
                ├── ${SAMPLE_NAME}.hitstats
                ├── ${SAMPLE_NAME}.idhistogram
                ├── ${SAMPLE_NAME}.idhistogram.svg
                ├── ${SAMPLE_NAME}.inserthistogram
                ├── ${SAMPLE_NAME}.inserthistogram.svg
                ├── ${SAMPLE_NAME}.mapratio.csv
                ├── ${SAMPLE_NAME}.mapratio.csv.svg
                ├── ${SAMPLE_NAME}.phyloFlash
                ├── ${SAMPLE_NAME}.phyloFlash.extractedSSUclassifications.csv
                ├── ${SAMPLE_NAME}.phyloFlash.html
                ├── ${SAMPLE_NAME}.phyloFlash.json
                ├── ${SAMPLE_NAME}.phyloFlash.NTUabundance.csv
                ├── ${SAMPLE_NAME}.phyloFlash.NTUabundance.csv.svg
                ├── ${SAMPLE_NAME}.phyloFlash.NTUfull_abundance.csv
                ├── ${SAMPLE_NAME}.phyloFlash.report.csv
                ├── ${SAMPLE_NAME}.phyloFlash.unassembled.NTUabundance.csv
                ├── ${SAMPLE_NAME}.remap_spades.bbmap.out
                ├── ${SAMPLE_NAME}.spades.out
                ├── ${SAMPLE_NAME}.spades_rRNAs.final.fasta
                ├── ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.1.fq
                ├── ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.2.fq
                ├── ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.sam
                ├── ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU_spades.sam
                ├── ${SAMPLE_NAME}.SSU.collection.alignment.fasta
                ├── ${SAMPLE_NAME}.SSU.collection.fasta
                ├── ${SAMPLE_NAME}.SSU.collection.fasta.tree
                ├── ${SAMPLE_NAME}.SSU.collection.fasta.tree.svg
                ├── ${SAMPLE_NAME}.toalign.fasta
                └── ${SAMPLE_NAME}-unprocessed.txt
```

| Filename | Description |
|----------|-------------|
| phyloflash.iqtree | Full result of the run, this is the main report file (a copy of *iqtree/16s.iqtree*) |
| phyloflash-summary.txt | The aggregated phyloFlash results of all samples |

### Directory Description
#### alignment
| Filename | Description |
|----------|-------------|
| phyloflash-alignment.fasta | The multiple sequence alignment produced by MAFFT. |
| phyloflash-matches.txt     | A list of reconstructed 16S genes and their match |

#### bactopia-info
| Filename | Description |
|----------|-------------|
| phyloflash-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| phyloflash-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| phyloflash-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### iqtree
Where possible descriptions were taken from 
IQ-TREE's [Command Reference](https://github.com/Cibiv/IQ-TREE/wiki/Command-Reference)
page, [Web Server Tutorial](https://github.com/Cibiv/IQ-TREE/wiki/Web-Server-Tutorial) page, 
and the [Tutorial](http://www.iqtree.org/doc/Tutorial) page.

| Filename | Description |
|----------|-------------|
| 16s.alninfo | Alignment site statistics |
| 16s.bionj | A neighbor joining tree produced by BIONJ |
| 16s.ckp.gz | IQ-TREE writes a checkpoint file |
| 16s.contree | Consensus tree with assigned branch supports where branch lengths are optimized on the original alignment; printed if Ultrafast Bootstrap is selected |
| 16s.iqtree | Full result of the run, this is the main report file |
| 16s.log | Run log |
| 16s.mldist | Contains the likelihood distances |
| 16s.model.gz | Information about all models tested |
| 16s.splits.nex | Support values in percentage for all splits (bipartitions), computed as the occurence frequencies in the bootstrap trees |
| 16s.treefile | Maximum likelihood tree in NEWICK format, can be visualized with treeviewer programs |
| 16s.ufboot | Trees created during the bootstrap steps |
| 16s.uniqueseq.phy | Unique sequences indentified by IQ-TREE |


#### samples
Where possible descriptions were taken from phyloFlash's [Output Summary](https://hrgv.github.io/phyloFlash/output.html) 
and the phyloFlash source [PhyloFlash.pm](https://github.com/HRGV/phyloFlash/blob/master/PhyloFlash.pm)

| Filename | Description |
|----------|-------------|
| ${SAMPLE_NAME}.all.dbhits.NR97.fa | Reference sequences from database with hits from the supplied reads, clustered at 97% identity |
| ${SAMPLE_NAME}.all.final.fasta | All assembled and reconstructed sequences from SPAdes in a single file |
| ${SAMPLE_NAME}.all.final.phyloFlash.dbhits.fa | FASTA file of all sequences in database with hits to reconstructed sequences |
| ${SAMPLE_NAME}.all.final.phyloFlash.notmatched.fa | FASTA file of full-length sequences without any database hits |
| ${SAMPLE_NAME}.all.vsearch.csv | CSV file of Vsearch output |
| ${SAMPLE_NAME}.assemratio.csv | CSV file of ratio assembled to unassembled |
| ${SAMPLE_NAME}.assemratio.csv.svg | A SVG image of the above ratios |
| ${SAMPLE_NAME}.bbmap.out | The bbmap log |
| ${SAMPLE_NAME}.bbmap.sam | The alignment of reads against 16S genes |
| ${SAMPLE_NAME}.hitstats | A SVG image of the above ratios |
| ${SAMPLE_NAME}.idhistogram | Histogram of the % identity of reads vs. reference database sequences, in tab-separated format |
| ${SAMPLE_NAME}.idhistogram.svg | A SVG image of the histogram above |
| ${SAMPLE_NAME}.inserthistogram | Histogram of detected insert sizes in tab-separated format, if paired-end reads were input |
| ${SAMPLE_NAME}.inserthistogram.svg | A SVG image of the above histogram |
| ${SAMPLE_NAME}.mapratio.csv | Ratios of mapped vs unmapped to report |
| ${SAMPLE_NAME}.mapratio.csv.svg | A SVG image of the above ratio |
| ${SAMPLE_NAME}.phyloFlash | Plain text file version of the HTML report |
| ${SAMPLE_NAME}.phyloFlash.extractedSSUclassifications.csv | Taxonomic classification of full-length sequences, in CSV format |
| ${SAMPLE_NAME}.phyloFlash.html | phyloFlash report file in HTML format, with a report on the taxonomic composition of SSU rRNA reads, quality metrics for the library, and affiliation of the reconstructed/assembled full-length sequences |
| ${SAMPLE_NAME}.phyloFlash.json | JSON version of *${SAMPLE_NAME}.phyloFlash* |
| ${SAMPLE_NAME}.phyloFlash.NTUabundance.csv | The list of uniqe higher level taxa (e.g. orders for bacteria) in the order of their appearance |
| ${SAMPLE_NAME}.phyloFlash.NTUabundance.csv.svg | A SVG image depicting the NTU abundances |
| ${SAMPLE_NAME}.phyloFlash.NTUfull_abundance.csv | NTU abundances (untruncated) from initial mapping, in CSV format |
| ${SAMPLE_NAME}.phyloFlash.report.csv | phyloFlash report in CSV format |
| ${SAMPLE_NAME}.phyloFlash.unassembled.NTUabundance.csv | Taxonomic composition of unassembled SSU reads in CSV format |
| ${SAMPLE_NAME}.remap_spades.bbmap.out | SAM file of re-mapping extracted reads to SPAdes full-length sequences |
| ${SAMPLE_NAME}.spades.out | The SPAdes log |
| ${SAMPLE_NAME}.spades_rRNAs.final.fasta | Assembled OTUs from SPAdes with phyloFlash simplified headers |
| ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.1.fq | The filtered SSU reads and their paired read, forward read file |
| ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.2.fq | The filtered SSU reads and their paired read, reverse read file |
| ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU.sam | SAM file of initial read mapping to SSU rRNA database |
| ${SAMPLE_NAME}.${SAMPLE_NAME}_R1.fastq.gz.SSU_spades.sam | SAM file of re-mapping extracted reads to SPAdes full-length sequences |
| ${SAMPLE_NAME}.SSU.collection.alignment.fasta | An aligned multifasta of all the predicted OTUs and the references |
| ${SAMPLE_NAME}.SSU.collection.fasta | A multifasta of all the predicted OTUs and the references |
| ${SAMPLE_NAME}.SSU.collection.fasta.tree | An NJ tree of the MAFFT alignment of all the predicted OTUs and the references |
| ${SAMPLE_NAME}.SSU.collection.fasta.tree.svg | An SVG image of the tree above |
| ${SAMPLE_NAME}.toalign.fasta | Sequences from the sample that were used in the MAFFT alignment |
| ${SAMPLE_NAME}-unprocessed.txt | Text file with reason for not processing sample |

## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

    --phyloflash STR        Directory containing a pre-built phyloFlash database.

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: phyloflash

    --outdir DIR            Directory to write results to
                                Default: ./

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 1440 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 4

phyloFlash Related Parameters:
    --download_phyloflash   Download the latest phyloFlash database, even it exists.
    
    --yes                   You acknowledge SILVAs license.

    --taxlevel INT          Level in the taxonomy string to summarize read counts per taxon.
                                Numeric and 1-based (i.e. "1" corresponds to "Domain").
                                Default: 6

    --phyloflash_opts STR   Extra phyloFlash options in quotes.
                                Default: ''

    --allow_multiple_16s    Include samples with multiple reconstructed 16S genes. Due to
                                high sequence similarity in true multi-copy 16S genes, it
                                is unlikely each copy will be reconstructed, instead only
                                one. In order to get more than one reconstructed 16S gene
                                there must be a significant difference in the sequence
                                identity. As a consequence, any samples that have multiple
                                16S genes reconstructed contain multiple different species
                                within their sequencing.
                                Default: Exclude samples with multiple 16S genes


MAFFT Related Parameters:
    --align_all             Include reconstructed 16S genes as well as the corresponding
                                reference 16S genes in the alignment.

    --mafft_opts STR        MAFFT options to include (in quotes).
                                Default: ''

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
Below is the command that was used to create the Conda environment.
```
conda create -y -n bactopia-phyloflash -c conda-forge -c bioconda \
    phyloflash \
    mafft \
    iqtree \
    pigz
```

## References
* __[Barrnap](https://github.com/tseemann/barrnap)__  
_Seemann, T. [Barrnap: Bacterial ribosomal RNA predictor](https://github.com/tseemann/barrnap)._  

* __[BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)__  
_Bushnell, B. [BBMap short read aligner, and other bioinformatic tools.](http://sourceforge.net/projects/bbmap/)_  

* __[Bedtools](https://github.com/arq5x/bedtools2)__  
_Quinlan, A. R. & Hall, I. M. [BEDTools: a flexible suite of utilities for 
comparing genomic features](http://dx.doi.org/10.1093/bioinformatics/btq033). 
Bioinformatics 26, 841–842 (2010)._  

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

* __[MAFFT](https://mafft.cbrc.jp/alignment/software/)__  
_Katoh, K. & Standley, D. M. 
[MAFFT multiple sequence alignment software version 7: improvements in performance and usability.](https://doi.org/10.1093/molbev/mst010) 
Mol. Biol. Evol. 30, 772–780 (2013)_  

* __[nhmmer](http://hmmer.org/)__  
_Wheeler, T. J. & Eddy, S. R. 
[nhmmer: DNA homology search with profile HMMs.](https://doi.org/10.1093/bioinformatics/btt403)
 Bioinformatics 29, 2487–2489 (2013)_  

* __[phyloFlash](https://github.com/HRGV/phyloFlash)__  
_H. R. Gruber-Vodicka, B.KB. Seah, E. Pruesse. 
[phyloFlash — Rapid SSU rRNA profiling and targeted assembly from metagenomes.](https://doi.org/10.1101/521922) 
bioRxiv 521922_  

* __[SILVA rRNA Database](https://www.arb-silva.de/)__  
_Quast, C. et al. 
[The SILVA ribosomal RNA gene database project: improved data processing and web-based tools.](https://doi.org/10.1093/nar/gks1219) 
Nucleic Acids Res. 41, D590–6 (2013)_

* __[SPAdes](https://github.com/ablab/spades)__  
_Bankevich, A., et al. 
[SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing.](https://doi.org/10.1089/cmb.2012.0021) 
Journal of computational biology 19.5 (2012): 455-477._  

* __[VSEARCH](https://github.com/torognes/vsearch)__  
_Rognes, T., Flouri, T., Nichols, B., Quince, C. & Mahé, F. 
[VSEARCH: a versatile open source tool for metagenomics.](https://doi.org/10.7717/peerj.2584)
 PeerJ 4, e2584 (2016)_  
