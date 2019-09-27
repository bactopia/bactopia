#!/bin/bash
set -e
set -u

avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\1/'`
bwa index !{query}
if [ "${avg_len}" -gt "70" ]; then
    bwa mem -M -t !{task.cpus} !{bwa_mem_opts} !{query} !{fq} > bwa.sam
else
    if [ "!{single_end}" == "true" ]; then
        bwa aln -f bwa.sai -t !{task.cpus} !{bwa_aln_opts} !{query} !{fq[0]}
        bwa samse -n !{params.bwa_n} !{bwa_samse_opts} !{query} bwa.sai !{fq[0]} > bwa.sam
    else
        bwa aln -f r1.sai -t !{task.cpus} !{bwa_aln_opts} !{query} !{fq[0]}
        bwa aln -f r2.sai -t !{task.cpus} !{bwa_aln_opts} !{query} !{fq[1]}
        bwa sampe -n !{params.bwa_n} !{bwa_sampe_opts} !{query} r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam
    fi
fi
samtools view !{f_value} -bS bwa.sam | samtools sort -o !{query_name}.bam -

# Write per-base coverage
genomeCoverageBed -ibam !{query_name}.bam -d > !{query_name}.coverage.txt
if [[ !{params.compress} == "true" ]]; then
    pigz --best -n -p !{task.cpus} !{query_name}.coverage.txt
fi

mkdir -p mapping/!{query_name}
mv !{query_name}.bam !{query_name}.cov.gz mapping/!{query_name}
rm bwa.sam
