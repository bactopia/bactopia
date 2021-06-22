nextflow.enable.dsl = 2

process ANNOTATE_GENOME {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "annotation/${sample}*"

    input:
    tuple val(sample), val(single_end), file(fq), file(fasta), file(total_contigs)
    file prokka_proteins
    file prodigal_tf

    output:
    file "annotation/${sample}*"
    tuple val(sample), file("annotation/${sample}.{ffn,ffn.gz}"),emit: PLASMID_BLAST,optional: true
    tuple val(sample),
        file("annotation/${sample}.{ffn,ffn.gz}"),
        file("annotation/${sample}.{faa,faa.gz}"),emit: ANTIMICROBIAL_RESISTANCE, optional: true
    file "${task.process}/*" optional true

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    contig_count = total_contigs.getName().replace('total_contigs_', '')
    genus = "Genus"
    species = "species"
    proteins = ""
    if (prokka_proteins.getName() != 'EMPTY_PROTEINS') {
        proteins = "--proteins ${prokka_proteins}"
        if (SPECIES.contains("-")) {
            genus = SPECIES.split('-')[0].capitalize()
            species = SPECIES.split('-')[1]
        } else {
            genus = SPECIES.capitalize()
            species = "spp."
        }
    }

    prodigal = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prodigal = "--prodigaltf ${prodigal_tf}"
    }

    compliant = params.compliant ? "--compliant" : ""
    locustag = "--locustag ${sample}"
    renamed = false
    // Contig ID must <= 37 characters
    if ("gnl|${params.centre}|${sample}_${contig_count}".length() > 37) {
        locustag = ""
        compliant = "--compliant"
        renamed = true
    }
    addgenes = params.nogenes ? "" : "--addgenes"
    addmrna = params.addmrna ? "--addmrna" : ""
    rawproduct = params.rawproduct ? "--rawproduct" : ""
    cdsrnaolap = params.cdsrnaolap ? "--cdsrnaolap" : ""
    norrna = params.norrna ? "--norrna" : ""
    notrna = params.notrna ? "--notrna" : ""
    rnammer = params.rnammer ? "--rnammer" : ""
    rfam = params.rnammer ? "--rfam" : ""
    template "annotate_genome.sh"

    stub:
    """
    mkdir annotation
    mkdir ${task.process}
    touch annotation/${sample}
    touch annotation/${sample}.ffn
    touch annotation/${sample}.ffn.gz
    touch annotation/${sample}.faa
    touch annotation/${sample}.faa.gz
    touch "${task.process}/${sample}"
    """
}


//###############
//Module testing
//###############

workflow test{
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        file(params.fq),
        file(params.fasta),
        file(params.total_contigs)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        file(params.prokka_proteins)
        )
    TEST_PARAMS_CH3 = Channel.of(
        file(params.prodigal_tf)
        )

    annotate_genome(TEST_PARAMS_CH,TEST_PARAMS_CH2,TEST_PARAMS_CH3)
}
