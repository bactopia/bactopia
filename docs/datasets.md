# Build Datasets
Bactopia can make use of many existing public datasets, as well as private datasets. The process of downloading, building, and (or) configuring these datasets for Bactopia has been automated.

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
A few folders for things such as calling variants, insertion sequences and primers are created that the user can manually populate. More information is available below!

## Setting Up
Included in Bactopia is the `setup-datasets.py` script (located in the `bin` folder) to automate the process of downloading and/or building these datasets.

### Quick Start
``` bash
bactopia datasets
```
This will set up Ariba datasets (`card` and `vfdb_core`), RefSeq Mash sketch, GenBank Sourmash Signatures, and PLSDB in the newly created `datasets` folder. By default, `datasets` is used for the output directory, but this can be changed with `--outdir` .


### A Single Bacterial Species
``` bash
bactopia datasets --species "Haemophilus influenzae" --include_genus
```

### Multiple Bacterial Species
You can also set up datasets for multiple bacterial species at a time. There are two options to do so.

#### Comma-Separated 
At runtime, you can separate the the different species
``` bash
bactopia datasets --species "Haemophilus influenzae,Staphylococcus aureus" --include_genus
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
bactopia datasets --species species.txt --include_genus
```

This will setup the MLST schema (if available) and a protein cluster FASTA file for each species in `species.txt`. 

## Usage
``` 
bactopia datasets [-h] [--outdir STR] [--skip_ariba] [--ariba STR]
                  [--species STR] [--skip_mlst] [--skip_prokka]
                  [--include_genus]
                  [--assembly_level {all,complete,chromosome,scaffold,contig}]
                  [--limit INT] [--accessions STR] [--identity FLOAT]
                  [--overlap FLOAT] [--max_memory INT] [--fast_cluster]
                  [--skip_minmer] [--skip_plsdb] [--prodigal_tf STR]
                  [--reference STR] [--mapping STR] [--genes STR]
                  [--proteins STR] [--primers STR] [--force_optional]
                  [--cpus INT] [--clear_cache] [--force]
                  [--force_ariba] [--force_mlst] [--force_prokka]
                  [--force_minmer] [--force_plsdb] [--keep_files]
                  [--available_datasets] [--depends] [--version]
                  [--verbose] [--silent]

bactopia datasets - Setup public datasets for Bactopia

optional arguments:
  -h, --help            show this help message and exit
  --outdir STR          Directory to write output. (Default ./datasets)

Ariba Reference Datasets:
  --skip_ariba          Skip setup of Ariba datasets
  --ariba STR           Comma separated list of Ariba datasets to download and
                        setup. Available datasets include: argannot, card,
                        ncbi, megares, plasmidfinder, resfinder,
                        srst2_argannot, vfdb_core, vfdb_full, virulencefinder
                        (Default: "vfdb_core,card") Use --available_datasets
                        to see the full list.

Bacterial Species:
  --species STR         Download available MLST schemas and completed genomes
                        for a given species or a list of species in a text
                        file.
  --skip_mlst           Skip setup of MLST schemas for each species

Custom Prokka Protein FASTA:
  --skip_prokka         Skip creation of a Prokka formatted fasta for each
                        species
  --include_genus       Include all genus members in the Prokka proteins FASTA
  --assembly_level {all,complete,chromosome,scaffold,contig}
                        Assembly levels of genomes to download (Default:
                        complete).
  --limit INT           If available completed genomes exceeds a given limit,
                        a random subsample will be taken. (Default 1000)
  --accessions STR      A list of RefSeq accessions to download.
  --identity FLOAT      CD-HIT (-c) sequence identity threshold. (Default:
                        0.9)
  --overlap FLOAT       CD-HIT (-s) length difference cutoff. (Default: 0.8)
  --max_memory INT      CD-HIT (-M) memory limit (in MB). (Default: unlimited
  --fast_cluster        Use CD-HIT's (-g 0) fast clustering algorithm, instead
                        of the accurate but slow algorithm.

Minmer Datasets:
  --skip_minmer         Skip download of pre-computed minmer datasets (mash,
                        sourmash)

PLSDB (Plasmid) BLAST/Sketch:
  --skip_plsdb          Skip download of pre-computed PLSDB datbases (blast,
                        mash)

Optional User Provided Datasets:
  --prodigal_tf STR     A pre-built Prodigal training file to add to the
                        species annotation folder. Requires a single species
                        (--species) and will replace existing training files.
  --reference STR       A reference genome (FASTA/GenBank (preferred)) file or
                        directory to be added to the optional folder for
                        variant calling. Requires a single species
                        (--species).
  --mapping STR         A reference sequence (FASTA) file or directory to be
                        added to the optional folder for mapping. Requires a
                        single species (--species).
  --genes STR           A gene sequence (FASTA) file or directory to be added
                        to the optional folder for BLAST. Requires a single
                        species (--species).
  --proteins STR        A protein sequence (FASTA) file or directory to be
                        added to the optional folder for BLAST. Requires a
                        single species (--species).
  --primers STR         A primer sequence (FASTA) file or directory to be
                        added to the optional folder for BLAST. Requires a
                        single species (--species).
  --force_optional      Overwrite any existing files in the optional folders

Custom Options:
  --cpus INT            Number of cpus to use. (Default: 1)
  --clear_cache         Remove any existing cache.
  --force               Forcibly overwrite existing datasets.
  --force_ariba         Forcibly overwrite existing Ariba datasets.
  --force_mlst          Forcibly overwrite existing MLST datasets.
  --force_prokka        Forcibly overwrite existing Prokka datasets.
  --force_minmer        Forcibly overwrite existing minmer datasets.
  --force_plsdb         Forcibly overwrite existing PLSDB datasets.
  --keep_files          Keep all downloaded and intermediate files.
  --available_datasets  List Ariba reference datasets and MLST schemas
                        available for setup.
  --depends             Verify dependencies are installed.

Adjust Verbosity:
  --version             show program's version number and exit
  --verbose             Print debug related text.
  --silent              Only critical errors will be printed.

example usage:
  bactopia datasets
  bactopia datasets --ariba 'vfdb_core'
  bactopia datasets --species 'Staphylococcus aureus' --include_genus
```

### Useful Parameters
#### --clear_cache
To prevent a PubMLST.org query every run, a list of available schemas is cached to `$HOME/.bactopia/datasets.json`. The cache expires after 15 days, but in case a new species has been made available `--clear_cache` will force a query of PubMLST.org.

#### --cpus
Increasing `--cpus` (it defaults to 1) is useful for speeding up the download and clustering steps.

#### --force*
If a dataset exists, it will only be overwritten if one of the `--force` parameters are used.
 
#### --include_genus
Completed RefSeq genomes are downloaded for a given species to be used for protein clustering. `--include_genus` will also download completed RefSeq genomes for each genus member.

#### --assembly_level
By default, only completed genomes are downloaded. `--assembly_level` allows you to set the minimum assembly level (e.g. complete, scaffold, contigs, etc...) to download.

#### --limit
For some species of bacteria there might be thousands of completed genomes available. For dataset creation, downloading thousands of completed genomes will be time consuming and like take up a significant amount of storage. To help in such cases `--limit` can be used to limit the downloads to a random subset of genomes. The default value for `--limit` has been set to 1000 genomes. In cases where `--include_genus` is used, the random subsample will always include at least one genome from the given `--species` value.

#### --accessions
In cases where a random subset of completed genomes is not ideal, you can provide your own curated list of genomes to download with `--accessions`. The file should have a single NCBI RefSeq Assembly accession (E.g GCF_000008865) per line.

#### --keep_files
Many intermediate files are downloaded/created (e.g. completed genomes) and deleted during the building process, use `--keep_files` to retain these files.

#### Tweaking CD-HIT
There are parameters (`--identity`, `--overlap`, `--max_memory`, and `--fast_cluster`) to tweak CD-HIT if you find it necessary. Please keep in mind, the only goal of the protein clustering step is to help speed up Prokka, by providing a decent set of proteins to annotate against first.

## Datasets Folder Overview
After creating datasets you will have a directory structure that Bactopia recognizes. Based on the available datasets Bactopia will queue up the associated analyses.

Here is the directory structure for the Bactopia Datasets. Some of these include files from public datasets that can be used directly, but there are also other folders you can populate yourself to fit your needs.
```
${DATASET_FOLDER}
├── ariba
├── minmer
├── plasmid
└── species-specific
    └── ${SPECIES}
        ├── annotation
        │   ├── cdhit-stats.txt
        │   ├── genome_size.json
        │   ├── ncbi-metadata.txt
        │   ├── proteins.faa
        │   ├── proteins.faa.clstr
        │   └── proteins-updated.txt
        ├── minmer
        │   ├── minmer-updated.txt
        │   └── refseq-genomes.msh
        ├── mlst
        │   └── ${SCHEMA}
        │       ├── ariba.tar.gz
        │       ├── blastdb.tar.gz
        │       └── mlst-updated.txt
        └── optional
            ├── blast
            │   ├── genes
            │   │   └── ${NAME}.fasta
            │   ├── primers
            │   │   └── ${NAME}.fasta
            │   └── proteins
            │       └── ${NAME}.fasta
            ├── insertion-sequences
            │   └── ${NAME}.fasta
            ├── mapping-sequences
            │   └── ${NAME}.fasta
            └── reference-genomes
                └── ${NAME}.{gbk|fasta}
```
### General Datasets
General datasets can be used for all bacterial samples. There are three general dataset folders: `ariba`, `minmer` and `plasmid`.

The `ariba` folder contains pre-formatted datasets available from [Ariba's `getref` Reference Datasets](https://github.com/sanger-pathogens/ariba/wiki/Task:-getref). 

The `minmer` folder contains a [RefSeq Mash Sketch](https://mash.readthedocs.io/en/latest/data.html) and [GenBank Sourmash Signatures](https://sourmash.readthedocs.io/en/latest/datasets.html?highlight=--track-abundance#genbank-lca-dataset) of more than 100,000 genomes.

Finally, the `plasmid` folder contains a [PLSDB Mash Sketch & BLAST database](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/) from a curated set of plasmids.

!!! warning "Changing files in `ariba`, `minmer` and `plasmid` is not recommended"
    These directories are for general analysis and have been precomputed. Modifying these files
    may cause errors during analysis.

### Species Specific Datasets
Bactopia allows the datasets to be created for a specific species. The following sections outline the species specific datasets.

#### annotation
Completed RefSeq genomes are downloaded and then the proteins are clustered and formatted for usage with Prokka. The results from this clustering is stored in the `annotation` folder.  
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── annotation
            ├── cdhit-stats.txt
            ├── genome_size.json
            ├── ncbi-metadata.txt
            ├── prodigal.tf
            ├── proteins.faa
            ├── proteins.faa.clstr
            └── proteins-updated.txt
```

| Filename | Description |
|----------|-------------|
| cdhit-stats.txt | General statistics associated with CD-HIT clustering |
| genome_size.json | A list of genome size for each downloaded RefSeq genome |
| ncbi-metadata.txt | NCBI Assembly metadata associated with the downloaded RefSeq genomes |
| prodigal.tf | A pre-built species specific Prodigal training file provided with `--prodigal_tf` |
| proteins.faa | Set of [Prokka formatted proteins](https://github.com/tseemann/prokka#fasta-database-format) |
| proteins.faa.clstr | Description of the clusters created by CD-HIT |
| proteins-updated.txt | Information on the last time the protein set was updated |

!!! info "You can add your curated protein set here"
    If you have a set of proteins you would like to use for annotation, you can name it `proteins.faa` and place it 
    in the `annotation` folder. In order for your set of proteins to be used by Prokka, you must make sure you follow
    the [Prokka FASTA database format](https://github.com/tseemann/prokka#fasta-database-format).

    An alternative is to use the `--accessions` parameter and give `bactopia datasets` the list of RefSeq accessions when the
    dataset is created. In doing so the custom protein set will be automatically formatted using the genomes you specified.

#### minmer
By default, a [Mash sketch](https://github.com/marbl/Mash) is created for the completed genomes downloaded for clustering proteins. These sketches are then be used for automatic selection of reference genomes for variant calling.
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── minmer
            ├── minmer-updated.txt
            └── refseq-genomes.msh
```

| Filename | Description |
|----------|-------------|
| minmer-updated.txt | Information on the last time the mash sketch was updated |
| refseq-genomes.msh | A Mash sketch (k=31) of the RefSeq completed genomes |


!!! info "You can add your curated RefSeq sketch here"
    You can replace `refseq-genomes.msh` with a custom set of RefSeq genomes to be used for automatic reference selection. 
    The only requirements to do so are that only RefSeq genomes (start with GCF) are used and the `mash sketch` uses a k-mer
    length of 31 (`-k 31`). This will allow it to be compatible with Bactopia.

    An alternative is to use the `--accessions` parameter and give `bactopia datasets` the list of RefSeq accessions when the
    dataset is created. In doing so the mash sketch will be automatically created.
  

#### mlst
The `mlst` folder contains MLST schemas that have been formatted to be used by [Ariba](https://github.com/sanger-pathogens/ariba/) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── mlst
            └── ${SCHEMA}
                ├── ariba.tar.gz
                ├── blastdb.tar.gz
                └── mlst-updated.txt

```

| Filename | Description |
|----------|-------------|
| ariba.tar.gz | An Ariba formatted MLST dataset for a given schema |
| blastdb.tar.gz | A BLAST formatted MLST dataset for a given schema |
| mlst-updated.txt | Contains time stamp for the last time the MSLT dataset was updated |

!!! info "How does Bactopia handle organisms with multiple MLST schemas?"
    In a few cases, an organism might have multiple MLST schemas available (Example: [E. coli](https://pubmlst.org/escherichia/)). In such cases, each MLST schema is downloaded and set up. Bactopia will also call sequence types against each schema. 

!!! warning "Changing files in `mlst` is not recommended"
    The MLST schemas have been pre-formatted for your usage. There might be rare cases where you would like to provide your own schema. If this is the case it is recommended you take a look at: [What about MLST not hosted at pubmlst.org?](https://github.com/sanger-pathogens/ariba/issues/185) then follow the directory structure for `mlst`.

#### optional
Built into the Bactopia dataset structure is the `optional` folder that you, the user, can populate for species specific analysis. These could include specific genes you might want BLASTed against your samples or a specific reference you want all your samples mapped to and variants called.

##### blast
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── optional
            └── blast
                ├── genes
                │   └── ${NAME}.fasta
                ├── primers
                │   └── ${NAME}.fasta
                └── proteins
                    └── ${NAME}.fasta
```

In the `blast` directory there are three more directories! 

The `genes` folder is where you can place gene seqeunces (nucleotides) in FASTA format to query against assemblies using `blastn`.

The `primers` folder is where you can place primer sequences (nucleotides) in FASTA format to query against assemblies using `blastn`, but with primer-specific parameters and cut-offs.

Finally, the `proteins` (as you probably guessed!) is where you can place protein sequnces (amino acids) in FASTA format to query against assemblies using `blastp`.

##### mapping-sequences
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── optional
            └── mapping-sequences
                └── ${NAME}.fasta
```

In the `mapping-sequences` directory you can place FASTA files of any nucleotide sequence you would like FASTQ reads to be mapped against using [BWA](https://github.com/lh3/bwa). This can be useful if you are interested if whether a certain region or gene is covered or not.

##### reference-genomes
```
${DATASET_FOLDER}
└── species-specific
    └── ${SPECIES}
        └── optional
            └── reference-genomes
                └── ${NAME}.{gbk|fasta}
```

In the `reference-genomes` directory you can put a GenBank (preferred!) or FASTA file of a reference genome you would like variants to be called against using [Snippy](https://github.com/tseemann/snippy).
