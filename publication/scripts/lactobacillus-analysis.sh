#! /bin/bash
conda create -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia

# Build Lactobacillus dataset
bactopia datasets ${HOME}/bactopia-datasets \
    --species 'Lactobacillus' \
    --include_genus \
    --cpus 10

# Query ENA for all Lactobacillus sequence projects
bactopia search 1578 --prefix lactobacillus

# Process Lactobacillus samples
mkdir ${HOME}/bactopia
cd ${HOME}/bactopia
bactopia --accessions ${HOME}/lactobacillus-accessions.txt \
         --datasets ${HOME}/bactopia-datasets \
         --species lactobacillus \
         --coverage 100 \
         --cpus 4 \
         --min_genome_size 1000000 \
         --max_genome_size 4200000

# Use Bactopia Tools
# Create a summary of the results
bactopia tools summary --bactopia ${HOME}/bactopia --prefix lactobacillus

# Reconstruct 16S genes and create a phylogeny
bactopia tools phyloflash --phyloflash ${HOME}/bactopia-datasets/16s/138 \
                          --bactopia ${HOME}/bactopia \
                          --cpus 16 \
                          --exclude ${HOME}/bactopia-tool/summary/lactobacillus-exclude.txt

# Assign taxonmic classifications
bactopia tools gtdb --gtdb ${HOME}/bactopia-datasets/gtdb/db \
                    --bactopia ${HOME}/bactopia \
                    --cpus 48 \
                    --exclude ${HOME}/bactopia-tool/summary/lactobacillus-exclude.txt

# Determine ANI L. crispatus (GCF_003795065.1)
bactopia tools fastani --bactopia ${HOME}/bactopia \
                       --exclude ${HOME}/bactopia-tool/summary/lactobacillus-exclude.txt \
                       --accession GCF_003795065.1 \
                       --refseq_only \
                       --minFraction 0.0

# Select strains with >95% ANI
awk '{if ($3 > 95){print $0}}' ${HOME}/bactopia-tool/fastani/fastani.tsv | grep "RX" > ${HOME}/crispatus-include.txt

# Create pan-genome, and supplement with all completed L. crispatus genomes 
bactopia tools roary --bactopia ${HOME}/bactopia \
                     --cpus 20 \
                     --include ${HOME}/crispatus-include.txt \
                     --species "lactobacillus crispatus" \
                     --n
