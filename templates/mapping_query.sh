#!/bin/bash
set -e
set -u

avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\1/'`
bwa index !{query}
if [ "${avg_len}" -gt "70" ]; then
    bwa mem -M -t !{task.cpus} !{query} !{fq} > bwa.sam
else
    if [ "!{single_end}" == "true" ]; then
        bwa aln -f bwa.sai -t !{task.cpus} !{query} !{fq[0]}
        bwa samse -n 9999 !{query} bwa.sai !{fq[0]} > bwa.sam
    else
        bwa aln -f r1.sai -t !{task.cpus} !{query} !{fq[0]}
        bwa aln -f r2.sai -t !{task.cpus} !{query} !{fq[1]}
        bwa sampe -n 9999 !{query} r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam
    fi
fi
samtools view -bS bwa.sam | samtools sort -o !{query_name}.bam -

# Write per-base coverage
genomeCoverageBed -ibam !{query_name}.bam -d | \
    pigz --best -n -c -p !{task.cpus} - > !{query_name}.cov.gz

mkdir -p mapping/!{query_name}
mv !{query_name}.bam !{query_name}.cov.gz mapping/!{query_name}
rm bwa.sam
