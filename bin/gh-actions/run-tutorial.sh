#! /bin/bash
# A script to execute the Bactopia tutorial, https://bactopia.github.io/tutorial/.
# This script is only meant to be used with pytest and GitHub Actions.
set -e
set -x
CPUS=8
if [[ "${1}" == "docker" ]]; then
    PROFILE="-profile docker"
elif [[ "${1}" == "singularity" ]]; then
    PROFILE="-profile singularity"
fi

# Setup
mkdir bactopia-tutorial
cd bactopia-tutorial

# Build datasets
bactopia datasets \
    --ariba "vfdb_core,card" \
    --species "Staphylococcus aureus" \
    --include_genus \
    --limit 10 \
    --cpus ${CPUS}

# Get S. aureus datasets
git clone https://github.com/bactopia-datasets/staphylococcus-aureus.git
cp -r staphylococcus-aureus/species-specific/ datasets/
rm -rf staphylococcus-aureus/

# Single Sample Test
bactopia --accession SRX4563634 \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 30 \
         --genome_size median \
         --outdir ena-single-sample \
         --max_cpus ${CPUS} ${PROFILE}

# Multiple Samples
echo SRX4563687 > ena-accessions.txt
echo SRX4563689 >> ena-accessions.txt
bactopia --accessions ena-accessions.txt \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 30 \
         --genome_size median \
         --outdir ena-multiple-samples \
         --max_cpus ${CPUS} ${PROFILE}

# Local Samples
mkdir fastqs
cp ena-single-sample/SRX4563634/quality-control/SRX4563634*.fastq.gz fastqs/
cat fastqs/SRX4563634*.fastq.gz > fastqs/SRX4563634-SE.fastq.gz

# Single Local Samples
# Paired-end
bactopia --R1 fastqs/SRX4563634_R1.fastq.gz \
         --R2 fastqs/SRX4563634_R2.fastq.gz \
         --sample SRX4563634 \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 30 \
         --genome_size median \
         --outdir local-single-sample \
         --max_cpus ${CPUS} ${PROFILE}

# Single-end
bactopia --SE fastqs/SRX4563634-SE.fastq.gz \
         --sample SRX4563634-SE \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 30 \
         --genome_size median \
         --outdir local-single-sample \
         --max_cpus ${CPUS} ${PROFILE}

# Local Samples FOFN
bactopia prepare fastqs/ > fastqs.txt
bactopia --fastqs fastqs.txt \
         --datasets datasets/ \
         --species "Staphylococcus aureus" \
         --coverage 30 \
         --genome_size median \
         --outdir local-multiple-samples \
         --max_cpus ${CPUS} ${PROFILE}
