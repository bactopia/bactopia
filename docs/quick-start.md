# Quick Start
Here we go! No time to waste, let's get the ball rolling! Why are you still reading this?!? Go! Go! Go!

## Installation
```
conda create -y -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
```

## Build Dataset
```
bactopia datasets
```

This create a folder `./datasets` and will build the following datasets:

- [CARD](https://card.mcmaster.ca/)
- [VFDB](http://www.mgc.ac.cn/VFs/)
- [RefSeq Mash Sketch](https://mash.readthedocs.io/en/latest/data.html)
- [GenBank Sourmash Signatures](https://sourmash.readthedocs.io/en/latest/datasets.html?highlight=--track-abundance#genbank-lca-dataset)
- [PLSDB Mash Sketch & BLAST](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)

More information about these datasets is available at [Build Datasets](/datasets/).

## Run Bactopia!
On the first launch of Bactopia it will install the Conda environments, so expect some delays in doing so!

### Single Sample
#### Paired-End
```
bactopia --R1 SEQS_R1.fastq.gz \
         --R2 SEQS_R2.fastq.gz \
         --sample SAMPLE_NAME \
         --datasets datasets/ \
         --outdir OUTDIR
```

In the command above, be sure to replace *SEQS_R1.fastq.gz* and *SEQS_R2.fastq.gz* with the name of your FASTQ files. You will also want to replace *SAMPLE_NAME* with your sample's name and *OUTDIR* with a directory name you would like to use for results.

#### Single-End
```
bactopia --SE SEQS.fastq.gz --sample SAMPLE_NAME --datasets datasets/ --outdir OUTDIR
```

In the command above, be sure to replace *SEQS.fastq.gz* with the name of your FASTQ file. You will also want to replace *SAMPLE_NAME* with your sample's name and *OUTDIR* with a directory name you would like to use for results.

### Multiple Samples
```
bactopia prepare directory-of-fastqs/ > fastqs.txt
bactopia --fastqs fastqs.txt --datasets datasets --outdir OUTDIR
```

In the command above, be sure to replace *OUTDIR* with a directory name you would like to use for results.
