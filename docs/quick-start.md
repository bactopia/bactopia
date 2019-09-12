# Quick Start
Here we go! No time to waste, let's get the ball rolling! Why are you still reading this?!? Go! Go! Go!

## Installation
```
conda create -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
```

## Build Dataset
```
setup-datasets datasets
```


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
fastqs-fofn directory-of-fastqs/ > fastqs.txt
bactopia --fastqs fastqs.txt --dataset datasets --outdir ${OUTDIR}
```
