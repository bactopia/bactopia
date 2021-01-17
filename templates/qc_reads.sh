#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p quality-control
mkdir -p ${LOG_DIR}
ERROR=0
GENOME_SIZE=`head -n 1 !{genome_size}`
TOTAL_BP=$(( !{params.coverage}*${GENOME_SIZE} ))

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# BBMap (bbduk.sh, reformat.sh) Version" >> ${LOG_DIR}/!{task.process}.versions
bbduk.sh --version 2>&1 | grep " version" >> ${LOG_DIR}/!{task.process}.versions 2>&1

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --extra !{extra} --genome_size !{genome_size} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{extra} --genome_size !{genome_size}
    fi
fi

if [ "!{params.skip_qc}" == "true" ]; then
    echo "Sequence QC was skipped for !{sample}" > quality-control/!{sample}-qc-skipped.txt
    if [[ -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` quality-control/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` quality-control/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` quality-control/!{sample}.fastq.gz
        fi
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            cp !{fq[0]} quality-control/!{sample}_R1.fastq.gz
            cp !{fq[1]} quality-control/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            cp  !{fq[0]} quality-control/!{sample}.fastq.gz
        fi
    fi
else
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        # Remove Adapters
        bbduk.sh -Xmx!{qc_ram}g \
            in=!{fq[0]} in2=!{fq[1]} \
            out=adapter-r1.fq out2=adapter-r2.fq \
            ref=!{adapters} \
            k=!{params.adapter_k} \
            ktrim=!{params.ktrim} \
            mink=!{params.mink} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            threads=!{task.cpus} \
            ftm=!{params.ftm} \
            !{qin} ordered=t \
            stats=${LOG_DIR}/bbduk-adapter.log 1> ${LOG_DIR}/bbduk-adapter.out 2> ${LOG_DIR}/bbduk-adapter.err

        # Remove PhiX
        bbduk.sh -Xmx!{qc_ram}g \
            in=adapter-r1.fq in2=adapter-r2.fq \
            out=phix-r1.fq out2=phix-r2.fq \
            ref=!{phix} \
            k=!{params.phix_k} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            qtrim=!{params.qtrim} \
            trimq=!{params.trimq} \
            minlength=!{params.minlength} \
            minavgquality=!{params.maq} \
            !{qin} qout=!{params.qout} \
            tossjunk=!{params.tossjunk} \
            threads=!{task.cpus} \
            ordered=t \
            stats=${LOG_DIR}/bbduk-phix.log 1> ${LOG_DIR}/bbduk-phix.out 2> ${LOG_DIR}/bbduk-phix.err

        # Error Correction
        if [ "!{params.skip_error_correction}" == "false" ]; then
            echo "# Lighter Version" >> ${LOG_DIR}/!{task.process}.versions
            lighter -v >> ${LOG_DIR}/!{task.process}.versions 2>&1
            lighter -od . -r phix-r1.fq -r phix-r2.fq -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t !{task.cpus} 1> ${LOG_DIR}/lighter.out 2> ${LOG_DIR}/lighter.err
        else
            echo "Skipping error correction"
            ln -s phix-r1.fq phix-r1.cor.fq
            ln -s phix-r2.fq phix-r2.cor.fq
        fi

        # Reduce Coverage
        if (( ${TOTAL_BP} > 0 )); then
            reformat.sh -Xmx!{qc_ram}g \
                in=phix-r1.cor.fq in2=phix-r2.cor.fq \
                out=subsample-r1.fq out2=subsample-r2.fq \
                samplebasestarget=${TOTAL_BP} \
                sampleseed=!{params.sampleseed} \
                overwrite=t 1> ${LOG_DIR}/reformat.out 2> ${LOG_DIR}/reformat.err
        else
            echo "Skipping coverage reduction"
            ln -s phix-r1.cor.fq subsample-r1.fq
            ln -s phix-r2.cor.fq subsample-r2.fq
        fi

        # Compress
        pigz -p !{task.cpus} -c -n subsample-r1.fq > quality-control/!{sample}_R1.fastq.gz
        pigz -p !{task.cpus} -c -n subsample-r2.fq > quality-control/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        # Remove Adapters
        bbduk.sh -Xmx!{qc_ram}g \
            in=!{fq[0]} \
            out=adapter-r1.fq \
            ref=!{adapters} \
            k=!{params.adapter_k} \
            ktrim=!{params.ktrim} \
            mink=!{params.mink} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            threads=!{task.cpus} \
            ftm=!{params.ftm} \
            ordered=t \
            stats=${LOG_DIR}/bbduk-adapter.log 1> ${LOG_DIR}/bbduk-adapter.out 2> ${LOG_DIR}/bbduk-adapter.err

        # Remove PhiX
        bbduk.sh -Xmx!{qc_ram}g \
            in=adapter-r1.fq \
            out=phix-r1.fq \
            ref=!{phix} \
            k=!{params.phix_k} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            qtrim=!{params.qtrim} \
            trimq=!{params.trimq} \
            minlength=!{params.minlength} \
            minavgquality=!{params.maq} \
            qout=!{params.qout} \
            tossjunk=!{params.tossjunk} \
            threads=!{task.cpus} \
            ordered=t \
            stats=${LOG_DIR}/bbduk-phix.log 1> ${LOG_DIR}/bbduk-phix.out 2> ${LOG_DIR}/bbduk-phix.err

        # Error Correction
        if [ "!{params.skip_error_correction}" == "false" ]; then
            echo "# Lighter Version" >> ${LOG_DIR}/!{task.process}.versions
            lighter -v >> ${LOG_DIR}/!{task.process}.versions 2>&1
            lighter -od . -r phix-r1.fq -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t !{task.cpus} 1> ${LOG_DIR}/lighter.out 2> ${LOG_DIR}/lighter.err
        else
            echo "Skipping error correction"
            ln -s phix-r1.fq phix-r1.cor.fq
        fi

        # Reduce Coverage
        if (( ${TOTAL_BP} > 0 )); then
            reformat.sh -Xmx!{qc_ram}g \
                in=phix-r1.cor.fq \
                out=subsample-r1.fq \
                samplebasestarget=${TOTAL_BP} \
                sampleseed=!{params.sampleseed} \
                overwrite=t 1> ${LOG_DIR}/reformat.out 2> ${LOG_DIR}/reformat.err
        else
            echo "Skipping coverage reduction"
            ln -s phix-r1.cor.fq subsample-r1.fq
        fi

        # Compress
        pigz -p !{task.cpus} -c -n subsample-r1.fq > quality-control/!{sample}.fastq.gz
    fi

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate FASTQ files
        rm *.fq
    fi
fi

echo "# fastq-scan Version" >> ${LOG_DIR}/!{task.process}.versions
fastq-scan -v >> ${LOG_DIR}/!{task.process}.versions 2>&1
FINAL_BP=`gzip -cd quality-control/*.gz | fastq-scan | grep "total_bp" | sed -r 's/.*:[ ]*([0-9]+),/\1/'`
if [ ${FINAL_BP} -lt "!{params.min_basepairs}" ]; then
    ERROR=1
    echo "After QC, !{sample} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
            not exceed the required minimum !{params.min_basepairs} bp. Further analysis
            is discontinued." | \
    sed 's/^\s*//' > !{sample}-low-sequence-depth-error.txt
fi

FINAL_READS=`gzip -cd quality-control/*.gz | fastq-scan | grep "read_total" | sed -r 's/.*:[ ]*([0-9]+),/\1/'`
if [ ${FINAL_READS} -lt "!{params.min_reads}" ]; then
    ERROR=1
    echo "After QC, !{sample} FASTQ(s) contain ${FINAL_READS} total reads. This does
            not exceed the required minimum !{params.min_reads} reads count. Further analysis
            is discontinued." | \
    sed 's/^\s*//' > !{sample}-low-read-count-error.txt
fi

if [ "!{is_assembly}" == "true" ]; then
    touch quality-control/reads-simulated-from-assembly.txt
fi

if [ "${ERROR}" -eq "1" ]; then
    if [ "!{single_end}" == "false" ]; then
        mv quality-control/!{sample}_R1.fastq.gz quality-control/!{sample}_R1.error-fq.gz
        mv quality-control/!{sample}_R2.fastq.gz quality-control/!{sample}_R2.error-fq.gz
    else
        mv quality-control/!{sample}.fastq.gz quality-control/!{sample}.error-fq.gz
    fi
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
