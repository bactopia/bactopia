#!/bin/bash
set -e
set -u
if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p fastqs
    touch fastqs/!{sample}.fastq.gz
else
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

    mkdir -p fastqs
    if [ "!{single_end}" == "is_accession" ]; then
        # Download accession from ENA/SRA
        fastq-dl !{sample} $ARCHIVE \
            --cpus !{task.cpus} \
            --outdir fastqs/ \
            --group_by_experiment \
            --is_experiment $ASPERA $ASPERA_KEY $ASPERA_SPEED $FTP_ONLY
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    fi
fi
