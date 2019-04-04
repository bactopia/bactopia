#!/bin/bash
set -e
set -u
OUTPUT="genome-size.txt"

if [ "!{params.genome_size}" == "null" ]; then
    # Use mash to estimate the genome size, if a genome size cannot be
    # estimated set the genome size to 0
    mash sketch -o test -k 31 -m 3 !{fq[0]} 2>&1 | \
        grep "Estimated genome size:" | \
        awk '{if($4){printf("%d\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
    rm test.msh
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
