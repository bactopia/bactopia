#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
MERGED="multiple-read-sets-merged.txt"
mkdir -p fastqs
mkdir -p extra
mkdir -p ${LOG_DIR}

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

# Bactopia Version Info
echo "# Timestamp" > bactopia.versions
date --iso-8601=seconds >> bactopia.versions
echo "# Bactopia Version" >> bactopia.versions
echo "bactopia !{bactopia_version}" >> bactopia.versions
echo "# Nextflow Version" >> bactopia.versions
echo "nextflow !{nextflow_version}" >> bactopia.versions
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
if [ "!{sample_type}" == "paired-end" ]; then
    # Paired-End Reads
    ln -s `readlink !{r1[0]}` fastqs/!{sample}_R1.fastq.gz
    ln -s `readlink !{r2[0]}` fastqs/!{sample}_R2.fastq.gz
    touch extra/empty.fna.gz
elif [ "!{sample_type}" == "single-end" ]; then
    # Single-End Reads
    ln -s `readlink !{r1[0]}` fastqs/!{sample}.fastq.gz
    touch extra/empty.fna.gz
elif  [ "!{sample_type}" == "hybrid" ]; then
    # Paired-End Reads
    ln -s `readlink !{r1[0]}` fastqs/!{sample}_R1.fastq.gz
    ln -s `readlink !{r2[0]}` fastqs/!{sample}_R2.fastq.gz
    ln -s `readlink !{extra}` extra/!{sample}.fastq.gz
elif [ "!{sample_type}" == "merge-pe" ]; then 
    # Merge Paired-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "R1:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R1.fastq.gz
    echo "Merged R1:" >> ${MERGED}
    ls -l fastqs/!{sample}_R1.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

    echo "R2:" >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R2.fastq.gz
    echo "Merged R2:" >> ${MERGED}
    ls -l fastqs/!{sample}_R2.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

    touch extra/empty.fna.gz
elif [ "!{sample_type}" == "hybrid-merge-pe" ]; then 
    # Merge Paired-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "R1:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R1.fastq.gz
    echo "Merged R1:" >> ${MERGED}
    ls -l fastqs/!{sample}_R1.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

    echo "R2:" >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} >> ${MERGED}
    find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R2.fastq.gz
    echo "Merged R2:" >> ${MERGED}
    ls -l fastqs/!{sample}_R2.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

    ln -s `readlink !{extra}` extra/!{sample}.fastq.gz
elif [ "!{sample_type}" == "merge-se" ]; then 
    # Merge Single-End Reads
    echo "This sample had reads merged." > ${MERGED}
    echo "SE:" >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
    find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}.fastq.gz
    echo "Merged SE:" >> ${MERGED}
    ls -l fastqs/!{sample}.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

    touch extra/empty.fna.gz
elif [ "!{sample_type}" == "sra_accession" ]; then
    # Download accession from ENA/SRA
    FTP_ONLY="--ftp_only"
    ARCHIVE=""

    # Check if ascp is available
    if [ "!{use_ena}" == "true" ]; then
        ARCHIVE="ENA"
    else
        ARCHIVE="SRA"
    fi

    # fastq-dl Version
    echo "# fastq-dl Version" >> ${LOG_DIR}/!{task.process}.versions
    fastq-dl --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
        echo "Unable to download !{sample} from both SRA and ENA !{params.max_retry} times. This may or may 
              not be a temporary connection issue. Rather than stop the whole Bactopia run, 
              further analysis of !{sample} will be discontinued." | \
        sed 's/^\s*//' > !{sample}-fastq-download-error.txt
        exit
    else
        # Download accession from ENA/SRA
        fastq-dl !{sample} $ARCHIVE \
            --cpus !{task.cpus} \
            --outdir fastqs/ \
            --group_by_experiment \
            --is_experiment $FTP_ONLY > ${LOG_DIR}/fastq-dl.out 2> ${LOG_DIR}/fastq-dl.err
        touch extra/empty.fna.gz
    fi 
elif [ "!{is_assembly}" == "true" ]; then
    if [ "!{sample_type}" == "assembly_accession" ]; then
        # ncbi-genome-download Version
        echo "# ncbi-genome-download Version" >> ${LOG_DIR}/!{task.process}.versions
        ncbi-genome-download --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

        if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
            touch extra/empty.fna.gz
            echo "Unable to download !{sample} from NCBI Assembly !{params.max_retry} times. This may or may
                  not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                  further analysis of !{sample} will be discontinued." | \
            sed 's/^\s*//' > !{sample}-assembly-download-error.txt
            exit
        else
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
