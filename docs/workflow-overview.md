# Workflow Overview
Bactopia is an extensive workflow integrating numerous steps in bacterial genome analysis. Through out the workflow there are steps that are *always enabled* and *dataset enabled*. Each of the steps depicted in the image below are described in this section. 

A list of software **directly** used in each step is also listed. Please check out the [Acknowledgements](/acknowledgements/) section to get the full list of software as well how to download and cite said software.

![Bactopia Workflow](data/bactopia-workflow.png)

## Always Enabled Steps
The *Always Enabled Steps* are always executed by Bactopia. These steps do not depend of external datasets and thus are always enabled.

#### Gather FASTQs
Specifies exactly where the input FASTQ/FASTAs are coming from. If you are using local inputs (e.g. `--R1/--R2`, `--fastqs`) it will verify they can be accessed. 

If an accession(s) (`--accession` or `--accessions`) was given, the corresponding FASTQs (SRA/ENA) or assemblies (NCBI Assembly) are downloaded in this step. All assemblies will have 2x250bp Illumina reads simulated withour insertions or deletions and a minimum PHRED score of Q33.


| Software | Usage   |
|----------|:--------------|
| [ART](https://github.com/rpetit3/ena-dl)   | Generate simulated Illumina reads from an assembly |
| [ena-dl](https://github.com/rpetit3/ena-dl)   | Download FASTQ files from ENA |
| [ncbi-genome-download](https://github.com/rpetit3/ena-dl)   | Download GenBank/RefSeq assemblies from NCBI Assembly database |

#### Validate FASTQs
Determines if the FASTQ file contains enough sequencing to continue processing. The `--min_reads` and `--min_basepairs` parameters adjust the minimum amount of sequencing required to continue processing. This step does not *directly* test the validity of the FASTQ format (although, it would fail if the format is invalid!).

| Software | Usage |
|----------|:--------------|
| [fastq-scan](https://github.com/rpetit3/fastq-scan) | Determine total read and basepairs of FASTQ | 

#### Original Summary
Produces summary statistics (read lengths, quality scores, etc...) based on the original input FASTQs.

| Software | Usage |
|----------|:--------------|
| [FastQC](https://github.com/s-andrews/FastQC)   | Generates a HTML report of original FASTQ summary statistics | 
| [fastq-scan](https://github.com/rpetit3/fastq-scan) | Generates original FASTQ summary statistics in JSON format |

#### Genome Size
The genome size is by various programs in the Bactopia workflow. By default, if no genome size is given one is estimated using Mash. Otherwise, a specific genome size can be specified or completely disabled using the `--genome_size` parameter. See [Genome Size Parameter](/usage-basic/#-genome_size) to learn more about specifying the genome size. 

| Software | Usage |
|----------|:--------------|
| [Mash](https://github.com/marbl/Mash)     | If not given, estimates genome size of sample | 

#### Quality Control
The input FASTQs go through a few clean up steps. First, Illumina related adapters and phiX contaminants are removed. Then reads that fail to pass length and/or quality requirements are filtered out. If the genome size is available, sequence error-corrections are made and the total sequencing is reduced to a specified coverage. After this step, all downstream analyses are based on the QC'd FASTQ and the original is no longer used.

| Software | Description |
|----------|:--------------|
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)  | Removes Illumina adapters, phiX contaminants, filters reads based on length and quality score, and reduces inputs to a specified coverage. | 
| [Lighter](https://github.com/mourisl/Lighter)  | Corrects sequencing errors | 

#### QC Summary
Produces summary statistics (read lengths, quality scores, etc...) based on the final set of QC'd FASTQs.

| Software | Usage |
|----------|:--------------|
| [FastQC](https://github.com/s-andrews/FastQC)   | Generates a HTML report of QC'd FASTQ summary statistics | 
| [fastq-scan](https://github.com/rpetit3/fastq-scan) | Generates QC'd FASTQ summary statistics in JSON format |

#### Count 31-mers
All 31 basepair (31-mers) sequences are counted and the singletons (those 31-mers only counted once) are filtered out.

| Software | Description |
|----------|:--------------|
| [McCortex](https://github.com/mcveanlab/mccortex) | Counts 31-mers in the input FASTQ |

#### Minmer Sketch
A minmer sketch and signature is created based on the QC'd FASTQs for multiple values of *k*. If datasets are available, the sketches/signatures are used for further downstream analysis.

| Software | Usage   |
|----------|:--------------|
| [Mash](https://github.com/marbl/Mash) | Produces a sketch (*k*=21,31) of tje QC'd FASTQ |
| [Sourmash](https://github.com/dib-lab/sourmash) | Produces a signature (*k*=21,31,51) of the QC'd FASTQ |

#### De novo Assembly
The QC'd FASTQs are assembled using the [Shovill](https://github.com/tseemann/shovill) pipeline. This allows for a seamless assembly process using [MEGAHIT](https://github.com/voutcn/megahit), [SKESA](https://github.com/ncbi/SKESA), [SPAdes](https://github.com/ablab/spades) or [Velvet](https://github.com/dzerbino/velvet). Alternatively, if long reads are available to complement Illumina paired-end reads, hybrid assembly is available through [Unicycler](https://github.com/rrwick/Unicycler).

| Software | Usage   |
|----------|:--------------|
| [assembly-scan](https://github.com/rpetit3/assembly-scan) | Generates summary statistics of the final assembly |
| [Shovill](https://github.com/tseemann/shovill)  | Manages multiple steps in the Illumina assembly process |
| [Unicycler](https://github.com/rrwick/Unicycler) | Manages multiple steps in the hybrid assembly process |

#### Assembly Quality Assessment
After assembly, the *de novo* assembly is assessed for its biological (e.g. containment & contamination) as well as its technical (e.g. misassemblies and errors) quality using [CheckM](https://github.com/Ecogenomics/CheckM) and [QUAST](http://quast.sourceforge.net/).

| Software | Usage   |
|----------|:--------------|
| [CheckM](https://github.com/Ecogenomics/CheckM) | Assess the biological quality of a *de novo* assembly based on presence of marker genes |
| [QUAST](http://quast.sourceforge.net/) | Gives a summary on the technical (e.g. misassemblies etc) quality of a *de novo* assembly |


#### Genome Annotation
Genes are predicted and annotated from the assembled genome using [Prokka](https://github.com/tseemann/prokka). If available, a [clustered RefSeq protein set](/datasets/#species-specific) is used for the first pass of annotation.

| Software | Usage   |
|----------|:--------------|
| [Prokka](https://github.com/tseemann/prokka)   | Predicts and annotates assembled genomes |

#### Antimicrobial Resistance
Searches for antimicrobial resistance genes and assosiated point mutations in the annotated gene and protein sequences. If datasets are available, [local assemblies](/workflow-overview/#local-assembly) can also be used to predict antibiotic resistance.

| Software | Usage   |
|----------|:--------------|
| [AMRFinderPlus](https://github.com/ncbi/amr) | Predicts antimicrobial resistance based on genes and point mutations |

## Dataset Enabled Steps
The remaining *Dataset Enabled Steps* require supplemental datasets to be available to be executed. There are many datasets available that Bactopia can take advantage of. To learn more about setting up these datasets, check out [Build Datasets](/datasets/#build-datasets). These datasets can be broken into two groups, *Public Datasets* and *User Datasets*.

### Public Datasets
[Publicly available datasets](/datasets/#general) can be used for further analysis.

#### Call Variants (Auto)
Variants are predicted using [Snippy](https://github.com/tseemann/snippy). The QC'd FASTQs are aligned to the nearest (based on Mash distance) RefSeq completed genome. By default, only the nearest genome is selected, but multiple genomes can be selected (`--max_references`) or this feature can be completely disabled (`disable_auto_variants`).

| Software | Usage   |
|----------|:--------------|
| [Bedtools](https://github.com/arq5x/bedtools2) | Generates the per-base coverage of the reference alignment |
| [NCBI Genome Download](https://github.com/kblin/ncbi-genome-download) | Downloads the RefSeq completed genome |
| [Snippy](https://github.com/tseemann/snippy)   | Manages multiple steps in the haploid variant calling process |
| [vcf-annotator](https://github.com/rpetit3/vcf-annotator) | Adds annotations from reference GenBank to the final VCF |

#### Local Assembly
Using available [Ariba reference datasets](/acknowledgements/#ariba-reference-datasets), determines which reference sequences were found with an additional detailed report summarizing the results.

| Software | Usage   |
|----------|:--------------|
| [Ariba](https://github.com/sanger-pathogens/ariba)    | Creates local assemblies of reference sequences |

#### Minmer Query
Screens QC'd FASTQs and signatures against available [Minmer Datasets](/acknowledgements/#minmer-datasets).

| Software | Usage   |
|----------|:--------------|
| [Mash](https://github.com/marbl/Mash) | Screens against RefSeq and/or PLSDB sketches |
| [Sourmash](https://github.com/dib-lab/sourmash) | Screens signature against GenBank |

#### Sequence Type
Uses a [PubMLST.org](https://pubmlst.org/) MLST schema to determine the sequence type of the sample.

| Software | Usage   |
|----------|:--------------|
| [Ariba](https://github.com/sanger-pathogens/ariba)    | Runs QC'd FASTQ against a MLST database |
| [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)   | Aligns MLST loci against the assembled genome |

### User Datasets
Another option is for users to [provide their own data](/datasets/#user-populated-folders) to include in the analysis.

#### BLAST Alignment
Each gene, protein, or primer sequence provided by the user is aligned against the assembled genome.

| Software | Usage   |
|----------|:--------------|
| [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)   | Aligns reference sequences against the assembled genome |

#### Call Variants (User)
Uses the same procedure as [Call Variants (Auto)](/workflow-overview/#call-variants-auto), except variants are called against each reference provided by the user.

| Software | Usage   |
|----------|:--------------|
| [Bedtools](https://github.com/arq5x/bedtools2) | Generates the per-base coverage of the reference alignment |
| [Snippy](https://github.com/tseemann/snippy)   | Manages multiple steps in the haploid variant calling process |
| [vcf-annotator](https://github.com/rpetit3/vcf-annotator) | Adds annotations from reference GenBank to the final VCF |

#### Reference Mapping
Aligns the QC'd FASTQs to each sequence provided by the user.

| Software | Usage   |
|----------|:--------------|
| [Bedtools](https://github.com/arq5x/bedtools2) | Generates the per-base coverage of the reference alignment |
| [BWA](https://github.com/lh3/bwa/)      | Aligns QC'd FASTQ to a reference sequence |
| [Samtools](https://github.com/samtools/samtools) | Converts alignment from SAM to BAM |
