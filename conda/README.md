# Commands
Below is a list of the commands used to create each enviroment.

```
VERSION=1.2.1

# annotate_genome.yml
conda create -y -n bactopia-annotate_genome -c conda-forge -c bioconda prokka pigz

# antimicrobial_resistance.yml
conda create -y -n bactopia-annotate_genome -c conda-forge -c bioconda ncbi-amrfinderplus

# ariba_analysis.yml
conda create -y -n bactopia-ariba_analysis -c conda-forge -c bioconda ariba

# assemble_genome.yml
conda create -y -n bactopia-assemble_genome -c rpetit3 -c conda-forge -c bioconda shovill assembly-scan pigz

# call_variants.yml
conda create -y -n bactopia-call_variants -c conda-forge -c bioconda snippy vcf-annotator pigz

# count_31mers.yml
conda create -y -n bactopia-count_31mers -c conda-forge -c bioconda mccortex

# download_references.yml
conda create -y -n bactopia-download_references -c conda-forge -c bioconda ncbi-genome-download mash

# insertion_sequences.yml
conda create -y -n bactopia-insertion_sequences -c rpetit3 -c conda-forge -c bioconda ismapper=2.0.a

# gather_fastqs.yml
conda create -y -n bactopia-gather_fastqs -c rpetit3 -c conda-forge -c bioconda aspera-connect ena-dl

# minmers.yml
conda create -y -n bactopia-minmers -c conda-forge -c bioconda mash sourmash

# qc_reads.yml
conda create -y -n bactopia-qc_reads -c conda-forge -c bioconda bbmap fastqc fastq-scan lighter pigz

# sequence_type.yml
conda create -y -n bactopia-sequence_type -c conda-forge -c bioconda ariba blast
```

Environments were exported using `conda env export`
