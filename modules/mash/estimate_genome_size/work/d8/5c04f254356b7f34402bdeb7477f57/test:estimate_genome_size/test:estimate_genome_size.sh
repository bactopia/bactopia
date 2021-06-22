#!/bin/bash
set -e
set -u
OUTPUT="SRR2838702-genome-size.txt"
LOG_DIR="test:estimate_genome_size"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/test:estimate_genome_size.versions
date --iso-8601=seconds >> ${LOG_DIR}/test:estimate_genome_size.versions

# Verify AWS files were staged
if [[ ! -L "input.1" ]]; then
    if [ "false" == "true" ]; then
        check-staging.py --fq1 input.1 --extra input.2 --is_single
    else
        check-staging.py --fq1 input.1 --fq2 null --extra input.2
    fi
fi

if [ "1" == "null" ]; then
    # Use mash to estimate the genome size, if a genome size cannot be
    # estimated set the genome size to 0
    echo "# Mash Version" >> ${LOG_DIR}/test:estimate_genome_size.versions
    mash --version >> ${LOG_DIR}/test:estimate_genome_size.versions 2>&1
    if [ "false" == "false" ]; then
        mash sketch -o test -k 31 -m 3 -r input.1 null 2>&1 | \
            grep "Estimated genome size:" | \
            awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
    else
        mash sketch -o test -k 31 -m 3 input.1 2>&1 | \
            grep "Estimated genome size:" | \
            awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
    fi
    rm -rf test.msh
    ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`

    if [ ${ESTIMATED_GENOME_SIZE} -gt "18040666" ]; then
        # Probably high coverage, try increasing number of kmer copies to 10
        if [ "false" == "false" ]; then
            mash sketch -o test -k 31 -m 10 -r input.1 null 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        else
            mash sketch -o test -k 31 -m 10 input.1 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        fi
        rm -rf test.msh
    elif [ ${ESTIMATED_GENOME_SIZE} -lt "100000" ]; then
        # Probably low coverage, try decreasing the number of kmer copies to 1
        if [ "false" == "false" ]; then
            mash sketch -o test -k 31 -m 1 -r input.1 null 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        else
            mash sketch -o test -k 31 -m 1 input.1 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        fi
        rm -rf test.msh
    fi

    ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`
    if [ ${ESTIMATED_GENOME_SIZE} -gt "18040666" ]; then
        rm ${OUTPUT}
        echo "SRR2838702 estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
                allowed genome size (18040666 bp). If this is unexpected, please
                investigate SRR2838702 to determine a cause (e.g. metagenomic, contaminants, etc...).
                Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
                of SRR2838702 will be discontinued." | \
        sed 's/^\s*//' > SRR2838702-genome-size-error.txt
    elif [ ${ESTIMATED_GENOME_SIZE} -lt "100000" ]; then
        rm ${OUTPUT}
        echo "SRR2838702 estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
                allowed genome size (100000 bp). If this is unexpected, please
                investigate SRR2838702 to determine a cause (e.g. metagenomic, contaminants, etc...).
                Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
                of SRR2838702 will be discontinued." | \
        sed 's/^\s*//' > SRR2838702-genome-size-error.txt
    fi
else
    # Use the genome size given by the user. (Should be >= 0)
    echo "1" > ${OUTPUT}
fi

# pass along FASTQs
mkdir -p fastqs
if [[ -L "input.1" ]]; then
    if [ "false" == "false" ]; then
        # Paired-End Reads
        ln -s `readlink input.1` fastqs/SRR2838702_R1.fastq.gz
        ln -s `readlink null` fastqs/SRR2838702_R2.fastq.gz
    else
        # Single-End Reads
        ln -s `readlink input.1` fastqs/SRR2838702.fastq.gz
    fi
else
    if [ "false" == "false" ]; then
        # Paired-End Reads
        cp input.1 fastqs/SRR2838702_R1.fastq.gz
        cp null fastqs/SRR2838702_R2.fastq.gz
    else
        # Single-End Reads
        cp  input.1 fastqs/SRR2838702.fastq.gz
    fi
fi


if [ "false" == "false" ]; then 
    cp .command.err ${LOG_DIR}/test:estimate_genome_size.err
    cp .command.out ${LOG_DIR}/test:estimate_genome_size.out
    cp .command.sh ${LOG_DIR}/test:estimate_genome_size.sh || :
    cp .command.trace ${LOG_DIR}/test:estimate_genome_size.trace || :
else
    rm -rf ${LOG_DIR}/
fi
