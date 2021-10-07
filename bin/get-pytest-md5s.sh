#! /bin/bash

WORK_DIR=$(basename $1)

find $1 -type f | sort | xargs -I {} md5sum {} | awk '{print "    - path: "$2"\n      md5sum: "$1}' | \
    sed "s/.run$/.run\n      contains: \['NEXTFLOW TASK', '\$NXF_ENTRY'\]/;s/.trace$/.trace\n      contains: \['nextflow.trace'\]/" | \
    sed "s=${1}=${WORK_DIR}/="
