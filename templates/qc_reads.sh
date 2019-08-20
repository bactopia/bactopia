#!/bin/bash
set -e
set -u


GENOME_SIZE=`head -n 1 !{genome_size}`
TOTAL_BP=$(( !{params.coverage}*${GENOME_SIZE} ))
if [ "!{single_end}" == "false" ]; then
    # Paired-End Reads
    # Remove Adapters
    bbduk.sh -Xmx4g \
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
        ordered=t \
        stats=quality-control/logs/bbduk-adapter.log

    # Remove PhiX
    bbduk.sh -Xmx!{params.xmx} \
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
        qout=!{params.qout} \
        tossjunk=!{params.tossjunk} \
        threads=!{task.cpus} \
        ordered=t \
        stats=quality-control/logs/bbduk-phix.log

    # Error Correction
    if [ "!{params.skip_error_correction}" == "false" ]; then
        lighter -od . -r phix-r1.fq -r phix-r2.fq -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t !{task.cpus}
    else
        echo "Skipping error correction"
        ln -s phix-r1.fq phix-r1.cor.fq
        ln -s phix-r2.fq phix-r2.cor.fq
    fi

    # Reduce Coverage
    if (( ${TOTAL_BP} > 0 )); then
        reformat.sh -Xmx!{params.xmx} \
            in=phix-r1.cor.fq in2=phix-r2.cor.fq \
            out=subsample-r1.fq out2=subsample-r2.fq \
            samplebasestarget=${TOTAL_BP} \
            sampleseed=!{params.sampleseed} \
            overwrite=t
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
    bbduk.sh -Xmx4g \
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
        stats=quality-control/logs/bbduk-adapter.log

    # Remove PhiX
    bbduk.sh -Xmx!{params.xmx} \
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
        stats=quality-control/logs/bbduk-phix.log

    # Error Correction
    if [ "!{params.skip_error_correction}" == "false" ]; then
        lighter -od . -r phix-r1.fq -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t !{task.cpus}
    else
        echo "Skipping error correction"
        ln -s phix-r1.fq phix-r1.cor.fq
    fi

    # Reduce Coverage
    if (( ${TOTAL_BP} > 0 )); then
        reformat.sh -Xmx!{params.xmx} \
            in=phix-r1.cor.fq \
            out=subsample-r1.fq \
            samplebasestarget=${TOTAL_BP} \
            sampleseed=!{params.sampleseed} \
            overwrite=t
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

FINAL_BP=`zcat quality-control/*.gz | fastq-scan | grep "total_bp" | sed -r 's/.*:([0-9]+),/\1/'`
if [ ${FINAL_BP} -lt "!{params.min_basepairs}" ]; then
    rm -f quality-control/*.fastq.gz
    echo "After QC, !{sample} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
          not exceed the required minimum !{params.min_basepairs} bp. Further analysis
          is discontinued." | \
    sed 's/^\s*//' > quality-control/low-sequence-depth-after-qc-error.txt
fi

FINAL_READS=`zcat quality-control/*.gz | fastq-scan | grep "read_total" | sed -r 's/.*:([0-9]+),/\1/'`
if [ ${FINAL_READS} -lt "!{params.min_reads}" ]; then
    rm -f quality-control/*.fastq.gz
    echo "After QC, !{sample} FASTQ(s) contain ${FINAL_READS} total reads. This does
          not exceed the required minimum !{params.min_reads} reads count. Further analysis
          is discontinued." | \
    sed 's/^\s*//' > quality-control/low-read-count-after-qc-error.txt
fi
