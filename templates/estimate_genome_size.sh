#!/bin/bash
set -e
set -u

OUTPUT="!{sample}-genome-size.txt"
LOG_DIR="!{task.process}"
if [ "!{params.dry_run}" == "true" ]; then
    touch ${OUTPUT}
else
    mkdir -p ${LOG_DIR}
    touch ${LOG_DIR}/!{task.process}.versions
    if [ "!{params.genome_size}" == "null" ]; then
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
    elif [ "!{params.genome_size}" == "min" ]; then
        # Use the minimum genome size based on completed RefSeq genomes
        echo "!{species_genome_size.min}" > ${OUTPUT}
    elif [ "!{params.genome_size}" == "median" ]; then
        # Use the median genome size based on completed RefSeq genomes
        echo "!{species_genome_size.median}" > ${OUTPUT}
    elif [ "!{params.genome_size}" == "mean" ]; then
        # Use the mean genome size based on completed RefSeq genomes
        echo "!{species_genome_size.mean}" > ${OUTPUT}
    elif [ "!{params.genome_size}" == "max" ]; then
        # Use the maximum genome size based on completed RefSeq genomes
        echo "!{species_genome_size.max}" > ${OUTPUT}
    else
        # Use the genome size given by the user. (Should be >= 0)
        echo "!{params.genome_size}" > ${OUTPUT}
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.run ${LOG_DIR}/!{task.process}.run
        cp .command.sh ${LOG_DIR}/!{task.process}.sh
        cp .command.trace ${LOG_DIR}/!{task.process}.trace
    else
        rm -rf ${LOG_DIR}/
    fi
fi
