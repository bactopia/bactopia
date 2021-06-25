nextflow.enable.dsl = 2

process ANTIMICROBIAL_RESISTANCE {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "logs/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${amrdir}/*"

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
    template "antimicrobial_resistance.sh"

    stub:
    amrdir = "antimicrobial-resistance"
    """
    mkdir ${amrdir}
    mkdir logs
    touch ${amrdir}/${sample}
    touch logs/${sample}
    """
}

//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        path(params.genes),
        path(params.proteins)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.amrdb)
        )
    antimicrobial_resistance(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect())
}
