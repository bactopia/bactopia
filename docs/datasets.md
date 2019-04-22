# Build Datasets
Bactopia can make use of many existing public datasets, as well as private datasets. The process of downloading, building, and (or) configuring these datasets for Bactopia has been automated for the user.

!!! info "Highly recommended to complete this step!"

    This step is completely optional, but it is highly recommended that you do not. By skipping this step of setting up public datasets, Bactopia will be limited to analyses like quality control, assembly, and 31-mer counting. 

## Included Datasets
Some datasets included are applicable to all bacterial species and some are specific to a bacterial species. If specified at runtime, Bactopia will recognize the datasets and execute the appropriate analyses.

### General
**[Ariba's `getref` Reference Datasets](https://github.com/sanger-pathogens/ariba/wiki/Task:-getref)**  
Allows reference datasets (resistance, virulence, and plamids) to be automatically downloaded and configured for usage by Ariba  

**[RefSeq Mash Sketch](https://mash.readthedocs.io/en/latest/data.html)**  
~100,000 genomes and plasmids from NCBI RefSeq, used to give an idea of what is your sequencing data (e.g. Are the sequences what you expected?)  

**[GenBank Sourmash Signatures](https://sourmash.readthedocs.io/en/latest/datasets.html?highlight=--track-abundance#genbank-lca-dataset)**  
~87,000 microbial genomes (includes viral and fungal) from NCBI GenBank, also gives an idea of what is your sequencing data.  

**[PLSDB Mash Sketch & BLAST](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/)**  
Includes meta data/annotations, Mash sketches, and BLAST database files of all plasmids stored in PLSDB.  

### Species Specific

**[PubMLST.org MLST Schemas](https://pubmlst.org/datasets/)**  
Multi-locus sequence typing (MLST) allelic profiles and seqeunces for a many different bacterial species (and even a few eukaryotes!).  

**Clustered RefSeq Proteins**  
For the given bacterial species, completed RefSeq genomes are downloaded and then the proteins are clustered and formatted for usage with Prokka.  

**Minmer Sketch of RefSeq Genomes**  
Using the completed genomes downloaded for clustering proteins a Mash sketch and Sourmash signature is created for these genomes. These sketches can then be used for automatic selection of reference genomes for variant calling.

**Optional User Populated Folders**  
A few folders for things such as calling variants, insertion sequences and primers are created that the user can manually populate. 

## Setting Up
Included in Bactopia is the `setup-datasets.py` script (located in the `bin` folder) to automate the process of downloading and/or building these datasets.

### Quick Start
``` bash
setup-datasets.py datasets
```

This will set up Ariba datasets (`card` and `vfdb_core`), RefSeq Mash sketch, GenBank Sourmash Signatures, and PLSDB in the newly created `datasets` folder.


### A Single Bacterial Species
``` bash
setup-datasets.py datasets --species "Haemophilus influenzae" --include_genus
```



### Multiple Bacterial Species
You can also set up datasets for multiple bacterial species at a time. There are two options to do so.

#### Comma-Separated 
At runtime, you can separate the the different species
``` bash
setup-datasets.py datasets --species "Haemophilus influenzae,Staphylococcus aureus" --include_genus
```
#### Text File

In order to do so, you will need to create a text file where each line is the name of a species to set up.

For example, you could create a `species.txt` file and include the following species in it.
``` bash
Haemophilus influenzae
Staphylococcus aureus
Mycobacterium tuberculosis
```

The new command becomes:

``` bash
setup-datasets.py datasets --species species.txt --include_genus
```

This will setup the MLST schema (if available) and a protein cluster FASTA file for each species in `species.txt`. 

## Usage
``` 
usage: setup-datasets.py [-h] [--ariba STR] [--species STR]
                              [--skip_prokka] [--include_genus]
                              [--identity FLOAT] [--overlap FLOAT]
                              [--max_memory INT] [--fast_cluster]
                              [--skip_minmer] [--skip_plsdb] [--cpus INT]
                              [--clear_cache] [--force] [--force_ariba]
                              [--force_mlst] [--force_prokka]
                              [--force_minmer] [--force_plsdb]
                              [--keep_files] [--list_datasets] [--depends]
                              [--version] [--verbose] [--silent]
                              OUTPUT_DIRECTORY

setup-datasets.py (v1.0.0) - Setup public datasets for Bactopia

positional arguments:
  OUTPUT_DIRECTORY  Directory to write output.

optional arguments:
  -h, --help        show this help message and exit

Ariba Reference Datasets:
  --ariba STR       Setup Ariba datasets for a given reference or a list of
                    references in a text file. (Default: card,vfdb_core)

Bacterial Species:
  --species STR     Download available (cg)MLST schemas and completed genomes
                    for a given species or a list of species in a text file.

Custom Prokka Protein FASTA:
  --skip_prokka     Skip creation of a Prokka formatted fasta for each species
  --include_genus   Include all genus members in the Prokka proteins FASTA
  --identity FLOAT  CD-HIT (-c) sequence identity threshold. (Default: 0.9)
  --overlap FLOAT   CD-HIT (-s) length difference cutoff. (Default: 0.8)
  --max_memory INT  CD-HIT (-M) memory limit (in MB). (Default: unlimited
  --fast_cluster    Use CD-HIT's (-g 0) fast clustering algorithm, instead of
                    the accurate but slow algorithm.

Minmer Datasets:
  --skip_minmer     Skip download of pre-computed minmer datasets (mash,
                    sourmash)

PLSDB (Plasmid) BLAST/Sketch:
  --skip_plsdb      Skip download of pre-computed PLSDB datbases (blast, mash)

Helpful Options:
  --cpus INT        Number of cpus to use. (Default: 1)
  --clear_cache     Remove any existing cache.
  --force           Forcibly overwrite existing datasets.
  --force_ariba     Forcibly overwrite existing Ariba datasets.
  --force_mlst      Forcibly overwrite existing MLST datasets.
  --force_prokka    Forcibly overwrite existing Prokka datasets.
  --force_minmer    Forcibly overwrite existing minmer datasets.
  --force_plsdb     Forcibly overwrite existing PLSDB datasets.
  --keep_files      Keep all downloaded and intermediate files.
  --list_datasets   List Ariba reference datasets and (cg)MLST schemas
                    available for setup.
  --depends         Verify dependencies are installed.

Adjust Verbosity:
  --version         show program's version number and exit
  --verbose         Print debug related text.
  --silent          Only critical errors will be printed.

example usage:
  setup-public-datasets.py outdir
  setup-public-datasets.py outdir --ariba 'card'
  setup-public-datasets.py outdir --species 'Staphylococcus aureus' --include_genus
```

### Useful Parameters
#### --clear_cache
To determine which MLST schemas are available, PubMLST.org is queryied. To prevent a query every run, a list of available schemas is cached to `$HOME/.bactopia/datasets.json`. The cache expires after 15 days, but in case a new species has been made available `--clear_cache` will force a requery of PubMLST.org.

#### --cpus
Increasing `--cpus` (it defaults to 1) is useful for speeding up the download and clustering steps.

#### --force*
If a dataset exists, it will only be overwritten if one of the `--force` parameters are used.
 
#### --include_genus
Completed RefSeq genomes are downloaded for a given species to be used for protein clustering. `--include_genus` will also download completed RefSeq genomes for each genus member.

#### --keep_files
Many intermediate files are downloaded/created (e.g. completed genomes) and deleted during the building process, use `--keep_files` to retain these files.

#### Tweaking CD-HIT
There are parameters (`--identity`, `--overlap`, `--max_memory`, and `--fast_cluster`) to tweak CD-HIT if you find it necessary. Please keep in mind, the only goal of the protein clustering step is to help speed up Prokka, by providing a decent set of proteins to annotate against first.
