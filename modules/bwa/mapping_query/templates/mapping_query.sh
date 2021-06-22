#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --extra !{query} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{query}
    fi
fi

avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\1/'`
ls *.fasta | xargs -I {} grep -H "^>" {} | awk '{print $1}' | sed 's/.fasta:>/\t/' > mapping.txt
cat *.fasta > multifasta.fa

echo "# bwa Version" >> ${LOG_DIR}/!{task.process}.versions
bwa 2>&1 | grep "Version" >> ${LOG_DIR}/!{task.process}.versions 2>&1
bwa index multifasta.fa > ${LOG_DIR}/bwa-index.out 2> ${LOG_DIR}/bwa-index.err
if [ "${avg_len}" -gt "70" ]; then
    bwa mem -M -t !{task.cpus} !{bwa_mem_opts} multifasta.fa !{fq} > bwa.sam
else
    if [ "!{single_end}" == "true" ]; then
        bwa aln -f bwa.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > ${LOG_DIR}/bwa-aln.out 2> ${LOG_DIR}/bwa-aln.err
        bwa samse -n !{params.bwa_n} !{bwa_samse_opts} multifasta.fa  bwa.sai !{fq[0]} > bwa.sam 2> ${LOG_DIR}/bwa-samse.err
    else
        bwa aln -f r1.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > ${LOG_DIR}/bwa-aln.out 2> ${LOG_DIR}/bwa-aln.err
        bwa aln -f r2.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[1]} >> ${LOG_DIR}/bwa-aln.out 2>> ${LOG_DIR}/bwa-aln.err
        bwa sampe -n !{params.bwa_n} !{bwa_sampe_opts} multifasta.fa  r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam 2> ${LOG_DIR}/bwa-sampe.err
    fi
fi
# Write per-base coverage
echo "# samtools Version" >> ${LOG_DIR}/!{task.process}.versions
samtools 2>&1 | grep "Version" >> ${LOG_DIR}/!{task.process}.versions 2>&1
samtools view -bS bwa.sam | samtools sort -o cov.bam - > ${LOG_DIR}/samtools.out 2> ${LOG_DIR}/samtools.err

echo "# bedtools Version" >> ${LOG_DIR}/!{task.process}.versions
bedtools --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
genomeCoverageBed -ibam cov.bam -d > cov.txt 2> ${LOG_DIR}/genomeCoverageBed.err
split-coverages.py mapping.txt cov.txt --outdir mapping

if [[ !{params.compress} == "true" ]]; then
    pigz --best -n -p !{task.cpus} mapping/*.txt
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
