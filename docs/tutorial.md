You should now have a directory named `datasets` that has all the available datasets to be used by Bactopia.# Tutorial
For this tutorial, we will attempt to replicate the [Staphopia](https://staphopia.emory.edu) analysis pipeline with Bactopia. 

We will use *S. aureus* samples associated with cystic fibrosis lung infections that were recently published (details below, shameless self plug!) and are available from BioProject accession [PRJNA480016](https://www.ebi.ac.uk/ena/data/view/PRJNA480016).

* *Bernardy, Eryn E., et al. ["Whole-Genome Sequences of Staphylococcus aureus Isolates from Cystic Fibrosis Lung Infections."](https://doi.org/10.1128/MRA.01564-18) Microbiol Resour Announc 8.3 (2019): e01564-18.*

Overall the goal of the tutorial is to:

- Build datasets
- Acquire Staphopia datasets
- Use Bactopia to process:
    - A sample from SRA/ENA
    - Multiple samples from SRA/ENA
    - Single local sample
    - Multiple local samples using FOFN

!!! error "Bactopia Should Be Installed"
    This tutorial assumes you have already installed Bactopia. If you have not, please check out how to at [Installation](installation.md).

## Build Datasets
First let's create a directory to work in and activate our Bactopia environment.
```
mkdir bactopia-tutorial
cd bactopia-tutorial
conda activate bactopia
```

Now we are ready to build our datasets!

```
bactopia datasets \
    --ariba "vfdb_core,card" \
    --species "Staphylococcus aureus" \
    --include_genus \
    --cpus 4
```

Let's review what is happening here.

`--ariba "vfdb_core,card"` says to download and setup the [VFDB Core](http://www.mgc.ac.cn/VFs/) and the [CARD](https://card.mcmaster.ca/) databases to be used by Ariba.

`--species "Staphylococcus aureus"` will download MLST schemas associated with *S. aureus* it will also download completed *S. aureus* genomes (RefSeq only) that are used to create a set of protein set for annotation, a Mash sketch for automatic variant calling to the nearest neighbor, and calculate genome size statistics.

`--include_genus` will also download completed genomes from the *Staphylococcus* genus that will be used for the protein set. These completed genomes **are not** used for the sketch creation or genome size calculation.

`--cpus 4` will use 4 cpus for downloading and the clustering step. Adjust this number according to your setup!

These datasets will be built into the default `datasets/` folder which can be changed using `--outdir`

!!! info "Use CARD over MEGARes"
    Staphopia v1 made use of MEGARes, for the purposes of this tutorial we are going to use the CARD database instead.

If all goes well, the newly created datasets are now available in the folder `datasets/`.

We have now completed the dataset creation step! Pat yourself on the back! Next we'll supplement these datasets with some optional *S. aureus* specific datasets.

## Acquire Staphopia Datasets
Staphopia includes a few *optional* datasets such as *S. aureus* N315 reference genome and SCCmec sequences (primers, proteins, full cassettes).

We will acquire these files using the [Bactopia Datasets](https://github.com/bactopia/bactopia-datasets) GitHub repository. For this tutorial a [Staphopia v1](https://github.com/bactopia/bactopia-datasets/tree/staphopia-v1) branch has been created, which includes this optional dataset. Now let's clone the repository.

```
git clone -b staphopia-v1 https://github.com/bactopia/bactopia-datasets.git
```

Next we'll copy the files into our recently built datasets folder and delete the `bactopia-datasets` repository since we no longer need it.

```
cp -r bactopia-datasets/species-specific/ datasets/
rm -rf bactopia-datasets/
```

~Voilà! 

That should be it. You should now have the Staphopia v1 datasets included with your recentely built datasets (e.g. *S. aureus* protein clusters, RefSeq sketch, etc...)

## Running Bactopia
OK! Get your servers started up! It is time to get processing!

### Samples on SRA
#### Single Sample
Let's start this by downloading a single sample from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA), and processing it through Bactopia.

```
bactopia --accession SRX4563634 \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 100 \
         --genome_size median \
         --cpus 2 \
         --outdir ena-single-sample
```

So, what's happening here?

`--accession SRX4563634` is telling Bactopia to download FASTQs associated with Exeriment accession SRX4563634.

`--datasets datasets/` tells Bactopia your pre-built datasets are in the folder `datasets`.

`--species "Staphylococcus aureus"` tells Bactopia, within the datasets folder, use the species specific dataset for *S. aureus*.

`--coverage 100` will limit the cleaned up FASTQ file to an estimated 100x coverage based on the genome size.

`--genome_size median` tells Bactopia to use the median genome size of completed *S. aureus* genomes. The minimum, maximum, median, and mean genome sizes were calculated during the dataset building step. If a genome size was not given, it would have been estimated by Mash.

`--cpus 2` tells Bactopia to use a maximum of 2 cpus per process. Adjust this parameter to fit your setup!

`--outdir ena-single-sample` tells Bactopia to dump the results into the `ena-single-sample` folder. Please keep in mind, this will not stop Nextflow from creating files (.nextflow, trace.txt, etc...) and directories (work and .nextflow/) within your current directory.

!!! info "Use --use_ena to download from ENA"
    If you append `--use_ena` to the command above the FASTQ files for SRX4563634 will be downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) (ENA) instead of SRA.

Once you launch this command, sit back, relax and watch the Nextflow give realtime updates for SRX4563634's analysis! 

The **approximate completion time is ~15-30 minutes** depending on the number of cpus given and download times from ENA.

Once complete, the results from numerous tools available to you in `ena-single-sample/SRX4563634/`. 

#### Multiple Samples
Now we are going to have Bactopia download and process 5 samples from ENA. To do this we can use the `bactopia search` function.

```
bactopia search PRJNA480016 --limit 5
```

This will produce three files: `ena-accessions.txt`, `ena-results.txt` and `ena-summary.txt`. To learn more about these files see [Generating Accession List](/usage-basic/#generating-accession-list).

For this tutorial, `ena-accessions.txt` is the file we need. It contains five Experiment accessions, a single one per line. Just like this:
```
SRX4563688
SRX4563687
SRX4563686
SRX4563689
SRX4563690
```

*Note: you may have 5 different accessions from the PRJNA480016 project.*

To process these samples, we will adjust our command used previously.

```
bactopia --accessions ena-accessions.txt \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 100 \
         --genome_size median \
         --cpus 2 \
         --outdir ena-multiple-samples
```

Instead of `--accession` we are now using `--accessions ena-accessions.txt` which tells Bactopia to read `ena-accessions.txt`, and download each Experiment accession from SRA (for ENA add `--use_ena`) and then process it.

At this point, you might want to go for a walk or make yourself a coffee! This step has an **approximate completion time of ~45-120 minutes**, which again is fully dependent on the cpus used and the download times from SRA (or ENA).

Once this is complete, the results for all five samples will be found in the `ena-multiple-samples` directory. Each sample will have there own folder of results.

### Local Samples
So for the local samples, we're going to recycle some of the samples we downloaded from SRA/ENA.

First let's make a directory to put the FASTQs into:
```
mkdir fastqs
```

Now let's move some the FASTQs from our SRX4563634 sample into this folder.
```
cp ena-single-sample/SRX4563634/quality-control/SRX4563634*.fastq.gz fastqs/
```

Finally let's also make a single-end version of SRX4563634.
```
cat fastqs/SRX4563634*.fastq.gz > fastqs/SRX4563634-SE.fastq.gz
```

OK! Now we are ready to continue the tutorial!

#### Single Sample
Again we'll mostly be using the same parameters as previous, but with a few new ones. To process a single sample you can use the `--R1`/`--R2` (paired-end), `--SE` (single-end), and `--sample` parameters.

##### Paired-End
For paired-end reads you will want to use `--R1`, `--R2`, and `--sample`. For this paired-end example we'll use SRX4563634 again which we've copied to the `fastqs` folder.

```
bactopia --R1 fastqs/SRX4563634_R1.fastq.gz \
         --R2 fastqs/SRX4563634_R2.fastq.gz \
         --sample SRX4563634 \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 100 \
         --genome_size median \
         --cpus 2 \
         --outdir local-single-sample
```

Now Bactopia will recognize the `--R1` and `--R2` parameters as paired-end reads and process. The `--sample` is required and will be used for naming the output.

Similar to the single SRA/ENA sample, the **approximate completion time is ~15-30 minutes** depending on the number of cpus given.

Once complete, results can be found in `local-single-sample/SRX4563634/`. 

##### Single-End
In the case of Illumina reads, you're very unlikely to produce single-end reads, but they do exist in the wild (early days of Illumina!). Nevertheless, because single-end reads do exist, single-end support was built into Bactopia.

To analyze single-end reads, the `--SE` parameter will replace `--R1` and `--R2`. 
```
bactopia --SE fastqs/SRX4563634-SE.fastq.gz \
         --sample SRX4563634-SE \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 100 \
         --genome_size median \
         --cpus 2 \
         --outdir local-single-sample
```

Now SRX4563634-SE will be processed as a single-end sample. For single-end processing there are some paired-end only analyses (e.g. error correction, insertion sequences) that will be skipped. This leads to single-end samples being processed a little bit faster than pair-end samples. But, the **approximate completion time is still ~15-30 minutes**.

Once complete, you'll the results from numerous tools available to you in `local-single-sample/SRX4563634-SE/`. 

If you made it this far, you're almost done!

#### Multiple Samples (FOFN)
Here we go! The final way you can process samples in Bactopia!

Bactopia allows you to give a text file describing the input samples. This file of file names (FOFN), contains sample names and location to associated FASTQs. The Bactopia FOFN format is described in detail at [Basic Usage -> Multiple Samples](/usage-basic/#multiple-samples).

First we'll need to prepare a FOFN describing the FASTQ files in our `fastqs` folder. We can use `bactopia prepare` to do so:

```
bactopia prepare fastqs/ > fastqs.txt
```

This command will try to create a FOFN for you. For this turorial, the FASTQ names are pretty straight forward and should produce a correct FOFN (or at least it should! ... hopefully!). If that wasn't the case for you, there are ways to [tweak `bactopia prepare`](/usage-basic/#generating-a-fofn).

Now we can use the `--fastqs` parameters to process samples in the FOFN.

```
bactopia --fastqs fastqs.txt \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 100 \
         --genome_size median \
         --cpus 2 \
         --outdir local-multiple-samples
```

We no longer need `--R1`, `--R2`, `--SE`, or `--sample` as the values for these parameters can be determined from the FOFN. 

Here it is, the final wait! This step has an **approximate completion time of ~45-120 minutes**. So, you will definitely want to go for a walk or make yourself a coffee! You've earned it! 

Once this is complete, the results for each sample (within their own folder) will be found in the `local-multiple-samples` directory.

!!! info "FOFN is more cpu efficient, making it faster"
    The real benefit of using the FOFN method to process multiple samples is Nextflow's queue system will make better use of cpus. Processing multiple samples one at a time (via `--R1`/`--R2` or `--SE`) will lead more instances of jobs waiting on other jobs to finish, during which cpus aren't being used.

## What's next?
That should do it! Hopefully you have succeeded (yay! 🎉) and would like to use Bactopia on your own data! 

In this tutorial we covered how to build datasets (`bactopia datasets`) and how process samples. We also covered the `bactopia search` and `bactopia prepare` to prepare file for multiple sample processing.

If your ran into any issues, please let me know by submitting a [GitHub Issue](https://github.com/bactopia/bactopia/issues).
