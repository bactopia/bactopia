nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "antimicrobial_resistance"

process ANTIMICROBIAL_RESISTANCE {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${meta.id}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(meta), path(genes), path(proteins)
    each path(amrdb)

    output:
    tuple val(meta), path("*{gene,protein}-{point-mutations,report}.txt"), emit: results
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    plus = params.amr_plus ? "--plus" : ""
    report_common = params.amr_report_common ? "--report_common" : ""
    organism_gene = ""
    organism_protein = ""
    if (params.amr_organism) {
        organism_gene = "-O ${params.amr_organism} --point_mut_all ${meta.id}-gene-point-mutations.txt"
        organism_protein = "-O ${params.amr_organism} --point_mut_all ${meta.id}-protein-point-mutations.txt"
    }
    '''
    if [[ !{params.skip_compression} == "false" ]]; then
        # Files passed to other modules
        gunzip -c !{meta.id}.faa.gz > !{meta.id}.faa
        gunzip -c !{meta.id}.ffn.gz > !{meta.id}.ffn
    fi

    tar -xzvf !{amrdb}
    amrfinder -n !{meta.id}.ffn \
            -d amrfinderdb/ \
            -o !{meta.id}-gene-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_gene} !{plus} !{report_common} > amrfinder-gene.stdout.txt 2> amrfinder-gene.stderr.txt

    amrfinder -p !{meta.id}.faa \
            -d amrfinderdb/ \
            -o !{meta.id}-protein-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_protein} !{plus} !{report_common} > amrfinder-protein.stdout.txt 2> amrfinder-protein.stderr.txt

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    antimicrobial_resistance:
        amrfinder:  $(echo $(amrfinder --version 2>&1))
    END_VERSIONS
    '''

    stub:
    """
    mkdir ${PUBLISH_DIR}
    touch ${PUBLISH_DIR}/${meta.id}
    """
}
