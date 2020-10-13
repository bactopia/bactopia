#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p fastqs
mkdir -p extra
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
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
    if [ "!{use_ena}" == "true" ]; then
        if [ "!{ftp_only}" == "false" ]; then
            if which ascp; then
                # ascp found
                ASCP=`which ascp`
                ASPERA="--aspera $ASCP"
                if readlink $ASCP; then
                    # ascp is symbolic link, need to get key
                    ASPERA_KEY="--aspera_key `readlink -f $ASCP | sed 's=bin/ascp$=etc/asperaweb_id_dsa.openssh='`"
                fi
            else
                #ascp not found, only use ftp
                FTP_ONLY="--ftp_only"
            fi
        else
            FTP_ONLY="--ftp_only"
        fi
        ARCHIVE="ENA"
    else
        ARCHIVE="SRA"
    fi

    # fastq-dl Version
    echo "# fastq-dl Version" >> ${LOG_DIR}/!{task.process}.versions
    fastq-dl --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Download accession from ENA/SRA
    fastq-dl !{sample} $ARCHIVE \
        --cpus !{task.cpus} \
        --outdir fastqs/ \
        --group_by_experiment \
        --is_experiment $ASPERA $ASPERA_KEY $ASPERA_SPEED $FTP_ONLY > ${LOG_DIR}/fastq-dl.out 2> ${LOG_DIR}/fastq-dl.err
    touch extra/empty.fna.gz    
elif [ "!{is_assembly}" == "true" ]; then
    if [ "!{sample_type}" == "assembly_accession" ]; then
        # ncbi-genome-download Version
        echo "# ncbi-genome-download Version" >> ${LOG_DIR}/!{task.process}.versions
        ncbi-genome-download --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

        # Verify Assembly accession
        check-assembly-accession.py !{sample} > accession.txt 2> ${LOG_DIR}/check-assembly-accession.txt

        if [ -s "accession.txt" ]; then
            # Download from NCBI assembly and simulate reads
            mkdir fasta/
            ncbi-genome-download bacteria -o ./ -F fasta -p !{task.cpus} \
                                        -s !{section} -A accession.txt -r 50 !{no_cache} > ${LOG_DIR}/ncbi-genome-download.out 2> ${LOG_DIR}/ncbi-genome-download.err
            find . -name "*!{sample}*.fna.gz" | xargs -I {} mv {} fasta/
            rename 's/(GC[AF]_\d+).*/$1.fna.gz/' fasta/*
            gzip -cd fasta/!{sample}.fna.gz > !{sample}-art.fna
        else
            cp ${LOG_DIR}/check-assembly-accession.txt !{sample}-assembly-accession-error.txt
            exit
        fi
    elif [ "!{sample_type}" == "assembly" ]; then
        if [ "!{is_compressed}" == "true" ]; then
            gzip -cd !{extra} > !{sample}-art.fna
        else 
            cat !{extra} > !{sample}-art.fna
        fi
    fi
    # ART Version
    echo "# ART Version" >> ${LOG_DIR}/!{task.process}.versions
    art_illumina --help | head -n 6 | tail -n 5 >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Simulate reads from assembly, reads are 250bp without errors
    art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov !{fcov} \
                    -ir 0 -ir2 0 -dr 0 -dr2 0 -rs !{params.sampleseed} \
                    -na -qL 33 -qU 40 -o !{sample}_R \
                    --id !{sample} -i !{sample}-art.fna > ${LOG_DIR}/art.out 2> ${LOG_DIR}/art.err

    mv !{sample}_R1.fq fastqs/!{sample}_R1.fastq 
    mv !{sample}_R2.fq fastqs/!{sample}_R2.fastq
    pigz -p !{task.cpus} --fast fastqs/*.fastq
    cp !{sample}-art.fna extra/!{sample}.fna
    pigz -p !{task.cpus} --best extra/!{sample}.fna
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
