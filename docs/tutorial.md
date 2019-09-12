# TODO Tutorial
For this tutorial, we will attempt to replicate the [Staphopia](https://staphopia.emory.edu) analysis pipeline as best we can in Bactopia. This will involve:

- Building datasets
- Acquiring Staphopia datasets
- Downloading *Staphylococcus aureus* genomes from ENA
- Single sample processing
- Multiple sample processing using FOFN

!!! error "Bactopia Should Be Installed"
    This tutorial assumes you have already installed Bactopia. If you have not, please check out how to at [Installation](installation.md).

## Getting Started
For starters, we'll create a empty directory to conduct this tutorial in.

```
mkdir ~/bactopia-tutorial
cd ~/bactopia-tutorial
```

## Building Datasets
The first thing we'll want to do is build our datasets!
```
setup-datasets.py datasets --ariba "vfdb_core,card" --species "Staphylococcus aureus" --include_genus
```

!!! info "Use CARD over MEGARes"
    Staphopia v1 made use of MEGAres, for the purposes of this tutorial we are going to use the CARD database instead.

You should now have a directory named `datasets` that has all the available datasets to be used by Bactopia.

## Staphopia Datasets
Staphopia includes a few *optional* datasets that we'll want to also include. These datasets include those related to variant calling and SCCmec typing.

We can acquire these files using the [Bactopia Datasets](https://github.com/bactopia/bactopia-datasets) GitHub repository. For this tutorial a [Staphopia v1](https://github.com/bactopia/bactopia-datasets/tree/staphopia-v1) branch has been created, which includes this optional datasets.

First let's clone the repository, then we'll move the files into our recently built datasets folder.

```
git clone https://github.com/bactopia/bactopia-datasets.git
cd bactopia-datasets
mv species-specific/ ~/bactopia-tutorial/datasets
```

~Voila! That should be it. You should not have the Staphopia v1 datasets included with your recentely built datasets (e.g. *S. aureus* protein clusters, RefSeq sketch, etc...)
## Example Dataset
For this tutorial, we will be using a few sequenced samples availble from the European Nucleotide Archive.

## Download Genomes
We will now use the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena) to download a few *S. aureus* genomes.


## Running Bactopia

### Single Sample

### Multiple Samples (FOFN)






