[![Gitlab pipeline status (branch)](https://img.shields.io/gitlab/pipeline/bactopia/bactopia/master?style=flat-square&logo=appveyor)](https://gitlab.com/bactopia/bactopia/pipelines/latest)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bactopia/badges/installer/conda.svg)](https://bioconda.github.io/recipes/bactopia/README.html) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bactopia/badges/downloads.svg)](https://anaconda.org/bioconda/bactopia)
[![Gitter](https://badges.gitter.im/bactopia/bactopia.svg)](https://gitter.im/bactopia/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![GitHub](https://img.shields.io/github/license/bactopia/bactopia)](https://raw.githubusercontent.com/bactopia/bactopia/master/LICENSE)

# Bactopia
Bactopia is an extensive workflow for processing Illumina sequencing of bacterial genomes. The goal of Bactopia is to process your data with a broad set of tools, so that you can get to the fun part of analyses quicker!

Bactopia can be split into three main parts: [Bactopia Datasets](https://bactopia.github.io/datasets/), [Bactopia Analysis Pipeline](https://bactopia.github.io/#bactopia-workflow), and [Bactopia Tools](https://bactopia.github.io/tools/bactopia-tools/).

![Bactopia Overview](docs/data/bactopia-overview.png)

Bactopia Datasets provide a framework for including many existing public datasets, as well as private datasets, into your analysis The process of downloading, building, and (or) configuring these datasets for Bactopia has been automated.

Bactopia Analysis Pipeline is the main *per-isolate* workflow in Bactopia. Built with  [Nextflow](https://www.nextflow.io/), input FASTQs (local or available from SRA/ENA) are pu through numerous analyses including: quality control, assembly, annotation, reference mapping, variant calling, minmer sketch queries, blast alignments, insertion site prediction, sequence typing, and more. The Bactopia Analysis Pipeline automatically selects which analyses to include based on the available Bactopia Datasets.

Bactopia Tools are a set a independent workflows for comparative analyses. The comparative analyses may include summary reports, pan-genome, or phylogentic tree construction. Using the [predictable output structure](https://bactopia.github.io/output-overview/) of Bactopia you can pick and choose which samples to include for processing with a Bactopia Tool.

Bactopia was inspired by [Staphopia](https://staphopia.emory.edu/), a workflow we (Tim Read and myself) released that targets *Staphylococcus aureus* genomes. Using what we learned from Staphopia and user feedback, Bactopia was developed from scratch with usability, portability, and speed in mind from the start.

# Documentation
Documentation for Bactopia is available at https://bactopia.github.io/. The documentation includes a tutorial replicating [Staphopia](https://staphopia.emory.edu) and a complete overview of Bactopia. I highly encourage you check it out!

# Quick Start
```
conda create -y -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
bactopia datasets datasets

# Paired-end
bactopia --R1 R1.fastq.gz --R2 R2.fastq.gz --sample SAMPLE_NAME \
         --dataset datasets/ --outdir OUTDIR

# Single-End
bactopia --SE SAMPLE.fastq.gz --sample SAMPLE --dataset datasets/ --outdir OUTDIR

# Multiple Samples
bactopia prepare MY-FASTQS/ > fastqs.txt
bactopia --fastqs fastqs.txt --dataset datasets --outdir OUTDIR

# Single ENA/SRA Experiment
bactopia --accession SRX000000 --dataset datasets --outdir OUTDIR

# Multiple ENA/SRA Experiments
bactopia search "staphylococcus aureus" > accessions.txt
bactopia --accessions accessions.txt --dataset datasets --outdir ${OUTDIR}
```

# Installation
Bactopia has **a lot** of tools built into its workflow. As you can imagine, all these tools lead to numerous dependencies, and navigating dependencies can often turn into a very frustrating process. With this in mind, from the onset Bactopia was developed to only include programs that are installable using [Conda](https://conda.io/en/latest/).

Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. In other words, it makes it super easy to get the tools you need installed! The [official Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is a good starting point for getting started with Conda. Bactopia has been tested using the [Miniconda installer](https://conda.io/en/latest/miniconda.html), but the [Anaconda installer](https://www.anaconda.com/distribution/) should work the same.

Once you have Conda all set up, you are ready to create an environment for Bactopia. To do so, you can use the following command:

```
conda create -n bactopia -c conda-forge -c bioconda bactopia
```

After a few minutes you will have a new conda environment suitably named *bactopia*. To activate this environment, you will can use the following command:

```
conda activate bactopia
```

And voilà, you are all set to get started processing your data!

# Please Cite Datasets and Tools
If you have used Bactopia in your work, please be sure to cite any datasets or tools you may have used. [A list of each dataset/tool used by Bactopia has been made available](https://bactopia.github.io/acknowledgements/). 

*If a citation needs to updated please let me know!*

A BibTeX file of each citation is also available at [Bactopia Datasets and Software BibTeX](docs/data/bactopia-datasets-software.bib)

# Acknowledgements
Bactopia is truly a case of *"standing upon the shoulders of giants"*. Nearly every component of Bactopia was created by others and made freely available to the public.

I would like to personally extend my many thanks and gratitude to the authors of these software packages and public datasets. If you've made it this far, I owe you a beer 🍻 (or coffee ☕!) if we ever encounter one another in person. Really, thank you very much!


# Feedback
Your feedback is very valuable! If you run into any issues using Bactopia, have questions, or have some ideas to improve Bactopia, I highly encourage you to submit it to the [Issue Tracker](https://github.com/bactopia/bactopia/issues).


# License
[MIT License](https://raw.githubusercontent.com/bactopia/bactopia/master/LICENSE)

# Citation
Not published yet

# Author 

* Robert A. Petit III
* Twitter: [@rpetit3](https://twitter.com/rpetit3)
