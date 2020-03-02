#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p mapping
    touch mapping/dry-run.coverage.txt
else
    avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\1/'`
    ls *.fasta | xargs -I {} grep -H "^>" {} | awk '{print $1}' | sed 's/.fasta:>/\t/' > mapping.txt
    cat *.fasta > multifasta.fa
    bwa index multifasta.fa
    if [ "${avg_len}" -gt "70" ]; then
        bwa mem -M -t !{task.cpus} !{bwa_mem_opts} multifasta.fa !{fq} > bwa.sam
    else
        if [ "!{single_end}" == "true" ]; then
            bwa aln -f bwa.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa  !{fq[0]}
            bwa samse -n !{params.bwa_n} !{bwa_samse_opts} multifasta.fa  bwa.sai !{fq[0]} > bwa.sam
        else
            bwa aln -f r1.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa  !{fq[0]}
            bwa aln -f r2.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa  !{fq[1]}
            bwa sampe -n !{params.bwa_n} !{bwa_sampe_opts} multifasta.fa  r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam
        fi
    fi
    # Write per-base coverage
    samtools view -bS bwa.sam | samtools sort -o cov.bam -
    genomeCoverageBed -ibam cov.bam -d > cov.txt
    split-coverages.py mapping.txt cov.txt --outdir mapping
    awk '{print "samtools view -b cov.bam "$2" > mapping/"$1".bam"}' mapping.txt | xargs -I {} sh -c '{}'
    if [[ !{params.compress} == "true" ]]; then
        pigz --best -n -p !{task.cpus} mapping/*.txt
    fi
fi
