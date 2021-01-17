#!/bin/bash
set -e
set -u
OUTPUT="!{sample}-genome-size.txt"
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --extra !{extra} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{extra}
    fi
fi

if [ "!{genome_size}" == "null" ]; then
    # Use mash to estimate the genome size, if a genome size cannot be
    # estimated set the genome size to 0
    echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
    mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    if [ "!{single_end}" == "false" ]; then
        mash sketch -o test -k 31 -m 3 -r !{fq[0]} !{fq[1]} 2>&1 | \
            grep "Estimated genome size:" | \
            awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
    else
        mash sketch -o test -k 31 -m 3 !{fq[0]} 2>&1 | \
            grep "Estimated genome size:" | \
            awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
    fi
    rm -rf test.msh
    ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`

    if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ]; then
        # Probably high coverage, try increasing number of kmer copies to 10
        if [ "!{single_end}" == "false" ]; then
            mash sketch -o test -k 31 -m 10 -r !{fq[0]} !{fq[1]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        else
            mash sketch -o test -k 31 -m 10 !{fq[0]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        fi
        rm -rf test.msh
    elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
        # Probably low coverage, try decreasing the number of kmer copies to 1
        if [ "!{single_end}" == "false" ]; then
            mash sketch -o test -k 31 -m 1 -r !{fq[0]} !{fq[1]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        else
            mash sketch -o test -k 31 -m 1 !{fq[0]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        fi
        rm -rf test.msh
    fi

    ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`
    if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ]; then
        rm ${OUTPUT}
        echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
                allowed genome size (!{params.max_genome_size} bp). If this is unexpected, please
                investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
                of !{sample} will be discontinued." | \
        sed 's/^\s*//' > !{sample}-genome-size-error.txt
    elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
        rm ${OUTPUT}
        echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
                allowed genome size (!{params.min_genome_size} bp). If this is unexpected, please
                investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
                of !{sample} will be discontinued." | \
        sed 's/^\s*//' > !{sample}-genome-size-error.txt
    fi
else
    # Use the genome size given by the user. (Should be >= 0)
    echo "!{genome_size}" > ${OUTPUT}
fi

# pass along FASTQs
mkdir -p fastqs
if [[ -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
    fi
else
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
        cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        cp  !{fq[0]} fastqs/!{sample}.fastq.gz
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
