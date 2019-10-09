#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    touch !{sample}.sig !{sample}-k31.msh !{sample}-k21.msh
else
    zcat !{fastq} | mash sketch -o !{sample}-k21 -k 21 -s !{params.mash_sketch} -r -I !{sample} -
    zcat !{fastq} | mash sketch -o !{sample}-k31 -k 31 -s !{params.mash_sketch} -r -I !{sample} -
    sourmash compute --scaled !{params.sourmash_scale} -o !{sample}.sig -p !{task.cpus} \
                     --track-abundance --merge !{sample} -k 21,31,51 !{fastq}
fi
