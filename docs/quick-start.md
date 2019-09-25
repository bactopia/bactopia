# Quick Start
Here we go! No time to waste, let's get the ball rolling! Why are you still reading this?!? Go! Go! Go!

## Installation
```
conda create -y -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
```

## Build Dataset
```
bactopia datasets datasets/
```

This will build the following datasets:

- [CARD](https://card.mcmaster.ca/)
- [VFDB](http://www.mgc.ac.cn/VFs/)
- [RefSeq Mash Sketch](https://mash.readthedocs.io/en/latest/data.html)
- [GenBank Sourmash Signatures](https://sourmash.readthedocs.io/en/latest/datasets.html?highlight=--track-abundance#genbank-lca-dataset)
- [PLSDB Mash Sketch & BLAST](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)

More information about these datasets is available at [Build Datasets](/datasets/).

## Run Bactopia!
### Single Sample
#### Paired-End
```
bactopia --R1 ${SAMPLE}_R1.fastq.gz --R2 ${SAMPLE}_R2.fastq.gz --sample ${SAMPLE} \
         --dataset datasets/ --outdir ${OUTDIR}
```

#### Single-End
```
bactopia --SE ${SAMPLE}.fastq.gz --sample ${SAMPLE} --dataset datasets/ --outdir ${OUTDIR}
```

### Multiple Samples
```
bactopia prepare directory-of-fastqs/ > fastqs.txt
bactopia --fastqs fastqs.txt --dataset datasets --outdir ${OUTDIR}
```
