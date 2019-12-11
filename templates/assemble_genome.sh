#!/bin/bash
set -e
set -u
OUTDIR=assembly
if [ "!{params.dry_run}" == "true" ]; then
    touch !{sample}.fna shovill.dry_run.txt
else
    GENOME_SIZE=`head -n 1 !{genome_size}`
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize ${GENOME_SIZE} \
            --outdir ${OUTDIR} \
            --force \
            --minlen !{params.min_contig_len} \
            --mincov !{params.min_contig_cov} \
            --namefmt "!{params.contig_namefmt}" \
            --keepfiles \
            --cpus !{task.cpus} \
            --ram !{shovill_ram} \
            --assembler !{params.assembler} \
            --noreadcorr !{opts} !{kmers} !{nostitch} !{nocorr}
    else
        # Single-End Reads
        shovill-se --se !{fq[0]} --depth 0 --gsize ${GENOME_SIZE} \
            --outdir ${OUTDIR} \
            --force \
            --minlen !{params.min_contig_len} \
            --mincov !{params.min_contig_cov} \
            --namefmt "!{params.contig_namefmt}" \
            --keepfiles \
            --cpus !{task.cpus} \
            --ram !{shovill_ram} \
            --assembler !{params.assembler} !{opts} !{kmers} !{nocorr}
    fi

    TOTAL_CONTIGS=`grep -c "^>" ${OUTDIR}/contigs.fa || true`
    if [ "${TOTAL_CONTIGS}" -gt "0" ]; then
        sed -i -r  's/^>(contig[0-9]+)(.*)/>gnl|\1|!{sample}\2/' ${OUTDIR}/contigs.fa
        mv ${OUTDIR}/contigs.fa ${OUTDIR}/!{sample}.fna
        assembly-scan ${OUTDIR}/!{sample}.fna > ${OUTDIR}/!{sample}.fna.json
        TOTAL_CONTIG_SIZE=`grep "total_contig_length" ${OUTDIR}/!{sample}.fna.json | sed -r 's/.*: ([0-9]+)/\1/'`
        if [ ${TOTAL_CONTIG_SIZE} -lt "!{params.min_genome_size}" ]; then
            mv ${OUTDIR}/!{sample}.fna ${OUTDIR}/!{sample}-error.fna
            mv ${OUTDIR}/!{sample}.fna.json ${OUTDIR}/!{sample}-error.fna.json
            echo "!{sample} assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                  size (!{params.min_genome_size} bp). If this is unexpected, please investigate !{sample} to
                  determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                  Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                  based analysis of !{sample} will be discontinued." | \
            sed 's/^\s*//' > !{sample}-assembly-error.txt
        fi

        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/!{sample}.fna
        fi
    else
        echo "!{sample} assembled successfully, but 0 contigs were formed. Please investigate
              !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
              outcome. Further assembly-based analysis of !{sample} will be discontinued." | \
        sed 's/^\s*//' > !{sample}-assembly-error.txt
    fi

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate files
        rm -fv ${OUTDIR}/shovill.bam* ${OUTDIR}/flash.extendedFrags* ${OUTDIR}/flash.notCombined* ${OUTDIR}/skesa.fasta.*
    fi
fi
