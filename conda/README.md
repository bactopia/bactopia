# Commands
Below is a list of the commands used to create each enviroment.

```
# annotate_genome.yml
conda create -y -n bactopia-annotate_genome -c conda-forge -c bioconda prokka

# ariba_analysis.yml
conda create -y -n bactopia-ariba_analysis -c conda-forge -c bioconda ariba

# assemble_genome.yml
conda create -y -n bactopia-assemble_genome -c rpetit3 -c conda-forge -c bioconda shovill assembly-scan

# call_variants.yml
conda create -y -n bactopia-call_variants -c conda-forge -c snippy vcf-annotator

# count_31mers.yml
conda create -y -n bactopia-count_31mers -c conda-forge -c bioconda mccortex

# insertion_sequences.yml
conda create -y -n bactopia-insertion_sequences -c rpetit3 -c conda-forge -c bioconda ismapper=2.0.a

# minmers.yml
conda create -y -n bactopia-minmers -c conda-forge -c bioconda mash sourmash

# qc_reads.yml
conda create -y -n bactopia-qc_reads -c conda-forge -c bioconda illumina-cleanup

# sequence_type.yml
conda create -y -n bactopia-sequence_type -c conda-forge -c bioconda ariba blast
```

Environments were exported using `conda env export`
