nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "annotate_genome"

process ANTIMICROBIAL_RESISTANCE {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample}"
    label "max_cpus"
    label "antimicrobial_resistance"

    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "logs/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${amrdir}/*"

    input:
    tuple val(sample), path(genes), path(proteins)
    each path(amrdb)

    output:
    path "${amrdir}/*"
    path "logs/*" optional true

    shell:
    amrdir = "antimicrobial-resistance"
    plus = params.amr_plus ? "--plus" : ""
    report_common = params.amr_report_common ? "--report_common" : ""
    organism_gene = ""
    organism_protein = ""
    if (params.amr_organism) {
        organism_gene = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-gene-point-mutations.txt"
        organism_protein = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-protein-point-mutations.txt"
    }
    '''
    LOG_DIR="logs/!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    if [[ !{params.compress} == "true" ]]; then
        gzip -cd !{genes} > !{sample}.ffn
        gzip -cd !{proteins} > !{sample}.faa
    fi
    echo !{amrdb}
    ls -lha !{amrdb}

    tar -xzvf !{amrdb}
    mkdir !{amrdir}

    # amrfinder Version
    echo "# amrfinder Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    amrfinder --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    amrfinder -n !{sample}.ffn \
            -d amrfinderdb/ \
            -o !{amrdir}/!{sample}-gene-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_gene} !{plus} !{report_common} > ${LOG_DIR}/amrfinder-gene.out 2> ${LOG_DIR}/amrfinder-gene.err

    amrfinder -p !{sample}.faa \
            -d amrfinderdb/ \
            -o !{amrdir}/!{sample}-protein-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_protein} !{plus} !{report_common} > ${LOG_DIR}/amrfinder-protein.out 2> ${LOG_DIR}/amrfinder-protein.err

    if [[ !{params.compress} == "true" ]]; then
        rm !{sample}.faa !{sample}.ffn
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
    amrdir = "antimicrobial-resistance"
    """
    mkdir ${amrdir}
    mkdir logs
    touch ${amrdir}/${sample}
    touch logs/${sample}
    """
}
