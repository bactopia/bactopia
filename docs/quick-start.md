# Quick Start
Here we go! No time to waste, let's get the ball rolling! Why are you still reading this?!?

## Installation
```
conda create -n bactopia -c rpetit3 bactopia
conda activate bactopia
```

## Build Dataset
```
setup-datasets.py datasets
```

## Make FOFN
```
fastqs-fofn.py fastqs/ > fastqs.txt
```

## Run Bactopia!
```
bactopia.nf --fastqs fastqs.txt --outdir processed
```
