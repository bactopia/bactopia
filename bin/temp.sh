

# Illumina-Cleanup
#
# HFLU83357_R1.fastq.gz  HFLU83357_R2.fastq.gz  logs  summary
conda create -y -n illumina-cleanup_env -c conda-forge -c bioconda illumina-cleanup
source activate illumina-cleanup_env
illumina-cleanup --fastqs fastqs.txt --outdir temp/ --coverage 100 --genome_size 1800000 --max_cpus 22 --cpus 7
source deactivate


# Contamination
mash
    - run mash screen on refseq

# Assembly
conda create -y -n assembly_env -c conda-forge -c bioconda shovill assembly-scan
conda activate assembly_env
shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize !{genome_size} --outdir . \
    --minlen 500 --cpus !{cpus} --assembler !{assembler} --noreadcorr --force
assembly-scan contigs.fa
source deactivate


# MLST
conda create -y -n mlst_env -c conda-forge -c bioconda ariba blast
conda activate mlst_env
ariba run !{ariba_mlst_ref} !{fq[0]} !{fq[1]} ariba --threads !{cpus} --verbose


mlst-blast
    - Add step to output ST as well

# kmers
conda create -y -n mccortex_env -c conda-forge -c bioconda mccortex
source activate mccortex_env
mccortex31 build -k 31 -s !{sample} -2 !{fq[0]}:!{fq[1]} -t !{cpus} -m 2000mb -q initial.ctx # Paired-End
mccortex31 build -k 31 -s !{sample} -1 !{fq[0]} -t 10 -m 2000mb -q OUTPUT.ctx # Single-End
mccortex31 clean -q -B 2 -U2 -T2 -m 2000mb  -o !{sample}.ctx initial.ctx # Clean up Cortex file (mostly remove singletons)
source deactivate

# MinHashes (21, 31, 51)
conda create -y -n minmers_env -c conda-forge -c bioconda mash sourmash
source activate minmers_env
mash sketch -o !{sample}-21 -k 21 -s 10000 -r -I !{sample} -p !{cpus} !{fq[0]} !{fq[1]}
mash sketch -o !{sample}-31 -k 31 -s 10000 -r -I !{sample} -p !{cpus} !{fq[0]} !{fq[1]}
sourmash compute --scaled 10000 -o !{sample}.sig -p !{cpus} --merge !{sample} -k 21,31,51 !{fq[0]} !{fq[1]}
# dashing currently has no conda setup
# https://github.com/dnbaker/dashing
source deactivate


# Annotation
conda create -y -n annotation_env -c conda-forge -c bioconda prokka
source activate annotation_env
prokka --cpus !{cpus} --outdir . --force --proteins !{proteins} \
    --prefix !{sample} --locustag !{sample} --centre BACTOPIA --addgenes \
    --mincontiglen 200 !{assembly}
source deactivate

# Ariba pre built DBs
ariba
    - loop through directory to determine available DBs, run all
    - at minimum, virulence (vfdb_core) and resistance (card)

# ISMapper
ismapper
    - loop through directory to determine available FASTA files, run all
