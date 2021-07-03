nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)

process BLAST {
    /*
    Query gene FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"
    label "max_cpus"
    label "blast"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blast/${query}/*.{json,json.gz}"

    input:
    tuple val(sample), path(blastdb)
    path(query)

    output:
    path("blast/${query}/*.{json,json.gz}")
    file "${task.process}/*" optional true

    shell:

    '''
    LOG_DIR="!{task.process}"
    OUTDIR=blast/!{query}
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    if [ "!{query}" == "proteins" ]; then
        echo "# tblastn Version" >> ${LOG_DIR}/!{task.process}.versions
        tblastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    else
        echo "# blastn Version" >> ${LOG_DIR}/!{task.process}.versions
        blastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    fi
    echo "# Parallel Version" >> ${LOG_DIR}/!{task.process}.versions
    parallel --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    mkdir -p ${OUTDIR}

    for fasta in !{query}/*; do
        type=`readlink -f ${fasta}`
        name=$(basename "${fasta%.*}")
        mkdir -p temp_json
        if [ "!{query}" == "genes" ]; then
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            blastn -db blastdb/!{sample} \
                -outfmt 15 \
                -evalue 1 \
                -perc_identity !{params.perc_identity} \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - \
                -out temp_json/${name}_{#}.json
        elif [ "!{query}" == "primers" ]; then
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            blastn -db blastdb/!{sample} \
                -outfmt 15 \
                -task blastn \
                -dust no \
                -word_size 7 \
                -perc_identity !{params.perc_identity} \
                -evalue 1 \
                -query - \
                -out temp_json/${name}_{#}.json
        else
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            tblastn -db blastdb/!{sample} \
                -outfmt 15 \
                -evalue 0.0001 \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - \
                -out temp_json/${name}_{#}.json
        fi

        merge-blast-json.py temp_json > ${OUTDIR}/${name}.json
        rm -rf temp_json

        if [[ !{params.compress} == "true" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.json
        fi
    done

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir ${task.process}
    mkdir blast/${query}
    touch ${task.process}/${sample}
    touch blast/${query}/${sample}.json
    touch blast/${query}/${sample}.json.gz
    """
}

process MAKE_BLASTDB {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"
    label "blast"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    tuple val(sample), path(fasta)

    output:
    path("blastdb/*")
    tuple val(sample), path("blastdb/*"), emit: BLAST_DB, optional:true
    file "${task.process}/*" optional true

    shell:
    '''
    LOG_DIR="!{task.process}"
    mkdir blastdb
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    echo "# makeblastdb Version" >> ${LOG_DIR}/!{task.process}.versions
    makeblastdb -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Verify AWS files were staged
    if [[ ! -L "!{fasta}" ]]; then
        check-staging.py --assembly !{fasta}
    fi

    if [[ !{params.compress} == "true" ]]; then
        gzip -cd !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
    else
        cat !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir blastdb
    mkdir ${task.process}
    touch blastdb/${sample}
    touch ${task.process}/${sample}
    """
}