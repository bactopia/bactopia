#!/bin/bash
set -e
set -u
if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p fastqs extra
    touch fastqs/!{sample}.fastq.gz extra/!{sample}.fna.gz
else
    mkdir -p fastqs
    mkdir -p extra
    if [ "!{sample_type}" == "paired-end" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        touch extra/empty.fna.gz
    elif [ "!{sample_type}" == "single-end" ]; then
        # Single-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        touch extra/empty.fna.gz
    elif  [ "!{sample_type}" == "hybrid" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        ln -s `readlink !{extra}` extra/!{sample}.fastq.gz
    elif [ "!{sample_type}" == "sra_accession" ]; then
        # Download accession from ENA/SRA
        ASPERA=""
        ASPERA_KEY=""
        ASPERA_SPEED="--aspera_speed !{params.aspera_speed}"
        FTP_ONLY=""
        ARCHIVE=""

        # Check if ascp is available
        if [ "!{params.use_ena}" == "true" ]; then
            if [ "!{params.ftp_only}" == "false" ]; then
                if which ascp; then
                    # ascp found
                    ASCP=`which ascp`
                    ASPERA="--aspera $ASCP"
                    if readlink $ASCP; then
                        # ascp is symbolic link, need to get key
                        ASPERA_KEY="--aspera_key `readlink -f $ASCP | sed 's=bin/ascp$=etc/asperaweb_id_dsa.openssh='`"
                    fi
                fi
            else
                FTP_ONLY="--ftp_only"
            fi
            ARCHIVE="ENA"
        else
            ARCHIVE="SRA"
        fi

        # Download accession from ENA/SRA
        fastq-dl !{sample} $ARCHIVE \
            --cpus !{task.cpus} \
            --outdir fastqs/ \
            --group_by_experiment \
            --is_experiment $ASPERA $ASPERA_KEY $ASPERA_SPEED $FTP_ONLY
        touch extra/empty.fna.gz    
    elif [ "!{is_assembly}" == "true" ]; then
        if [ "!{sample_type}" == "assembly_accession" ]; then
            # Download from NCBI assembly and simulate reads
            mkdir fasta/
            ncbi-genome-download bacteria -o ./ -F fasta -p !{task.cpus} \
                                          -s !{section} -A !{sample} -r 50
            find -name "*!{sample}*.fna.gz" | xargs -I {} mv {} fasta/
            if [ "!{section}" == 'refseq' ]; then
                rename 's/(GCF_\d+).*/$1.fna.gz/' fasta/*
            else
                rename 's/(GCA_\d+).*/$1.fna.gz/' fasta/*
            fi
            zcat fasta/!{sample_name}.fna.gz > !{sample_name}-art.fna
        elif [ "!{sample_type}" == "assembly" ]; then
            if [ "!{is_compressed}" == "true" ]; then
                zcat !{extra} > !{sample_name}-art.fna
            else 
                cat !{extra} > !{sample_name}-art.fna
            fi
        fi
        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov !{fcov} \
                     -ir 0 -ir2 0 -dr 0 -dr2 0 -rs !{params.sampleseed} \
                     -na -qL 33 -qU 40 -o !{sample_name}_R \
                     --id !{sample_name} -i !{sample_name}-art.fna

        mv !{sample_name}_R1.fq fastqs/!{sample_name}_R1.fastq 
        mv !{sample_name}_R2.fq fastqs/!{sample_name}_R2.fastq
        pigz -p !{task.cpus} --fast fastqs/*.fastq
        cp !{sample_name}-art.fna extra/!{sample_name}.fna
        pigz -p !{task.cpus} --best extra/!{sample_name}.fna
    fi
fi
