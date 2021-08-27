nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "antimicrobial_resistance"

process ANTIMICROBIAL_RESISTANCE {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.overwrite,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), path(genes), path(proteins)
    each path(amrdb)

    output:
    tuple val(sample), path("*{gene,protein}-{point-mutations,report}.txt"), emit: results
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    plus = params.amr_plus ? "--plus" : ""
    report_common = params.amr_report_common ? "--report_common" : ""
    organism_gene = ""
    organism_protein = ""
    if (params.amr_organism) {
        organism_gene = "-O ${params.amr_organism} --point_mut_all ${sample}-gene-point-mutations.txt"
        organism_protein = "-O ${params.amr_organism} --point_mut_all ${sample}-protein-point-mutations.txt"
    }
    '''
    if [[ !{params.skip_compression} == "false" ]]; then
        # Files passed to other modules
        gunzip -c !{sample}.faa.gz > !{sample}.faa
        gunzip -c !{sample}.ffn.gz > !{sample}.ffn
    fi

    tar -xzvf !{amrdb}
    amrfinder -n !{sample}.ffn \
            -d amrfinderdb/ \
            -o !{sample}-gene-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_gene} !{plus} !{report_common} > amrfinder-gene.stdout.txt 2> amrfinder-gene.stderr.txt

    amrfinder -p !{sample}.faa \
            -d amrfinderdb/ \
            -o !{sample}-protein-report.txt \
            --ident_min !{params.amr_ident_min} \
            --coverage_min !{params.amr_coverage_min} \
            --translation_table !{params.amr_translation_table} \
            --threads !{task.cpus} !{organism_protein} !{plus} !{report_common} > amrfinder-protein.stdout.txt 2> amrfinder-protein.stderr.txt

    # Capture versions
    amrfinder --version >> amrfinder.version.txt 2>&1
    '''

    stub:
    """
    mkdir ${PUBLISH_DIR}
    touch ${PUBLISH_DIR}/${sample}
    """
}
