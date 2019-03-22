#!/bin/bash
set -e
set -u

zcat !{fastq} | mash sketch -o !{sample}-21 -k 21 -s !{params.mash_sketch} -r -I !{sample} -
zcat !{fastq} | mash sketch -o !{sample}-31 -k 31 -s !{params.mash_sketch} -r -I !{sample} -
sourmash compute --scaled !{params.sourmash_scale} -o !{sample}.sig -p !{task.cpus} --merge !{sample} -k 21,31,51 !{fastq}
