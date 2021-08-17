nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "assemble_genome"

process ASSEMBLE_GENOME {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${sample}"
    label "max_cpu_75"
    label "assemble_genome"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "assembly/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${sample}-assembly-error.txt"

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra), path(genome_size)

    output:
    path "assembly/*"
    path "${sample}-assembly-error.txt" optional true
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("assembly/${sample}.{fna,fna.gz}"),emit: SEQUENCE_TYPE, optional:true
    tuple val(sample), val(single_end), path("assembly/${sample}.{fna,fna.gz}"), emit: MAKE_BLASTDB, optional: true
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("assembly/${sample}.{fna,fna.gz}"), path("total_contigs_*"),emit: ANNOTATION, optional:true
    tuple val(sample), path("assembly/${sample}.{fna,fna.gz}"), path(genome_size),emit: ASSEMBLY_QC, optional: true
    path "${PROCESS_NAME}/*" optional true

    shell:
    contig_namefmt = params.contig_namefmt ? params.contig_namefmt : "${sample}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0]
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    no_miniasm = params.no_miniasm ? "--no_miniasm" : ""
    no_rotate = params.no_rotate ? "--no_rotate" : ""
    no_pilon = params.no_pilon ? "--no_pilon" : ""
    keep = params.keep_all_files ? "--keep 3" : "--keep 1"
    use_original_assembly = null
    if (sample_type.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    '''
    OUTDIR=assembly
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    GENOME_SIZE=`head -n 1 !{genome_size}`
    if [[ "!{sample_type}" == "hybrid"  || "!{params.assembler}" == "unicycler" ]]; then
        EXTRA=""
        if [ "!{sample_type}" == "hybrid" ]; then
            EXTRA="-l !{extra}"
        fi
        echo "# unicycler Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        unicycler --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        unicycler -1 !{fq[0]} -2 !{fq[1]} \
            ${EXTRA} -o ${OUTDIR} \
            --no_correct \
            --min_fasta_length !{params.min_contig_len} \
            --threads !{task.cpus} \
            !{keep} --mode !{params.unicycler_mode} \
            !{no_miniasm} !{no_rotate} !{no_pilon} --min_polish_size !{params.min_polish_size} \
            --min_component_size !{params.min_component_size} \
            --min_dead_end_size !{params.min_dead_end_size} > ${LOG_DIR}/unicycler.out 2> ${LOG_DIR}/unicycler.err
        sed -r 's/^>([0-9]+)(.*)/>!{sample}_\\1\\2/' ${OUTDIR}/assembly.fasta > ${OUTDIR}/!{sample}.fna
        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/*.gfa
            pigz -n --best -p !{task.cpus} ${OUTDIR}/*.fasta
        fi
    elif [ "!{use_original_assembly}" == "true" ]; then
        mkdir ${OUTDIR}
        gzip -cd !{extra} > ${OUTDIR}/!{sample}.fna
    else
        echo "# shovill Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        shovill --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        shovill --check >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

        if [ "!{params.assembler}" == "spades" ]; then
            echo "# SPAdes Version (this assembler was used)" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
            spades.py --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        elif [ "!{params.assembler}" == "skesa" ]; then
            echo "# SKESA Version (this assembler was used)" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
            skesa --version 2>&1 | tail -n 1 >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        elif [ "!{params.assembler}" == "velvet" ]; then
            echo "# Velvet Version (this assembler was used)" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
            velvetg | grep "^Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        else
            echo "# MEGAHIT Version (this assembler was used)" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
            megahit --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        fi

        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            shovill --R1 !{fq[0]} --R2 !{fq[1]} --depth 0 --gsize ${GENOME_SIZE} \
                --outdir ${OUTDIR} \
                --force \
                --minlen !{params.min_contig_len} \
                --mincov !{params.min_contig_cov} \
                --namefmt "!{contig_namefmt}" \
                --keepfiles \
                --cpus !{task.cpus} \
                --ram !{shovill_ram} \
                --assembler !{params.assembler} \
                --noreadcorr !{opts} !{kmers} !{nostitch} !{nocorr} > ${LOG_DIR}/shovill.out 2> ${LOG_DIR}/shovill.err
        else
            # Single-End Reads
            shovill-se --se !{fq[0]} --depth 0 --gsize ${GENOME_SIZE} \
                --outdir ${OUTDIR} \
                --force \
                --minlen !{params.min_contig_len} \
                --mincov !{params.min_contig_cov} \
                --namefmt "!{contig_namefmt}" \
                --keepfiles \
                --cpus !{task.cpus} \
                --ram !{shovill_ram} \
                --assembler !{params.assembler} !{opts} !{kmers} !{nocorr} > ${LOG_DIR}/shovill.out 2> ${LOG_DIR}/shovill.err
        fi
        mv ${OUTDIR}/contigs.fa ${OUTDIR}/!{sample}.fna
        
        if [ "!{params.keep_all_files}" == "false" ]; then
            # Remove intermediate files
            rm -fv ${OUTDIR}/shovill.bam* ${OUTDIR}/flash.extendedFrags* ${OUTDIR}/flash.notCombined* ${OUTDIR}/skesa.fasta.* ${OUTDIR}/*.fq.gz 
        fi
    fi

    TOTAL_CONTIGS=`grep -c "^>" ${OUTDIR}/!{sample}.fna || true`
    touch "total_contigs_${TOTAL_CONTIGS}"
    if [ "${TOTAL_CONTIGS}" -gt "0" ]; then
        assembly-scan ${OUTDIR}/!{sample}.fna > ${OUTDIR}/!{sample}.fna.json 2> ${LOG_DIR}/assembly-scan.err
        TOTAL_CONTIG_SIZE=`grep "total_contig_length" ${OUTDIR}/!{sample}.fna.json | sed -r 's/.*: ([0-9]+)/\1/'`
        if [ ${TOTAL_CONTIG_SIZE} -lt "!{params.min_genome_size}" ]; then
            mv ${OUTDIR}/!{sample}.fna ${OUTDIR}/!{sample}-error.fna
            mv ${OUTDIR}/!{sample}.fna.json ${OUTDIR}/!{sample}-error.fna.json
            echo "!{sample} assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                    size (!{params.min_genome_size} bp). If this is unexpected, please investigate !{sample} to
                    determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-assembly-error.txt
        fi

        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/!{sample}.fna
        fi
    else
        echo "!{sample} assembled successfully, but 0 contigs were formed. Please investigate
                !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
                outcome. Further assembly-based analysis of !{sample} will be discontinued." | \
        sed 's/^\\s*//' > !{sample}-assembly-error.txt
    fi

    # pass the FASTQs along
    mkdir -p fastqs
    if [[ -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
            cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            cp  !{fq[0]} fastqs/!{sample}.fastq.gz
        fi
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{PROCESS_NAME}.err
        cp .command.out ${LOG_DIR}/!{PROCESS_NAME}.out
        cp .command.sh ${LOG_DIR}/!{PROCESS_NAME}.sh || :
        cp .command.trace ${LOG_DIR}/!{PROCESS_NAME}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir assembly
    mkdir fastqs
    mkdir ${PROCESS_NAME}
    touch total_contigs_${sample}
    touch ${sample}-assembly-error.txt
    touch fastqs/${sample}.fastq.gz
    touch assembly/${sample}
    touch assembly/${sample}.fna
    touch assembly/${sample}.fna.gz
    touch ${PROCESS_NAME}/${sample}
    """
}
