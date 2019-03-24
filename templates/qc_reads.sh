#!/bin/bash
set -e
set -u

parse_params() {
   if [ "$2" == "null" ]; then
        echo ""
    else
        echo "$1 $2"
    fi
}

OUTDIR=$(parse_params "--outdir" !{params.outdir})
ADAPTERS=$(parse_params "--adapters" !{params.adapters})
PHIX=$(parse_params "--phix" !{params.phix})
GENOME_SIZE=`head -n 1 !{genome_size_file}`

printf "sample\tr1\tr2\n!{sample}\t!{fq[0]}\t!{fq2}\n" > temp-fastqs.txt

illumina-cleanup --fastqs temp-fastqs.txt \
    --coverage !{params.coverage} \
    --genome_size ${GENOME_SIZE} \
    --max_cpus !{params.max_cpus} \
    --cpus !{cpus} \
    --adapter_k !{params.adapter_k} \
    --phix_k !{params.phix_k} \
    --ktrim !{params.ktrim} \
    --mink !{params.mink} \
    --hdist !{params.hdist} \
    --tpe !{params.tpe} \
    --tbo !{params.tbo} \
    --qtrim !{params.qtrim} \
    --trimq !{params.trimq} \
    --maq !{params.maq} \
    --minlength !{params.minlength} \
    --ftm !{params.ftm} \
    --tossjunk !{params.tossjunk} \
    --qout !{params.qout} \
    --xmx !{params.xmx} \
    --maxcor !{params.maxcor} \
    --sampleseed !{params.sampleseed} \
    $OUTDIR \
    $ADAPTERS \
    $PHIX

mv !{sample}/ quality-control/
