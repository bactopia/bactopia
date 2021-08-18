nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "annotate_genome"

process ANNOTATE_GENOME {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.overwrite,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), val(single_end), file(fq), file(fasta), file(total_contigs)
    file prokka_proteins
    file prodigal_tf

    output:
    tuple val(sample), path("${sample}.faa"), emit: faa
    tuple val(sample), path("${sample}.ffn"), emit: ffn
    path "results/*", emit: results
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
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
    '''
    if [ "!{renamed}" == "true" ]; then
        echo "Original sample name (!{sample}) not used due to creating a contig ID >37 characters"
    fi

    prokka --outdir results \
        --force \
        --prefix '!{sample}' \
        --genus '!{genus}' \
        --species '!{species}' \
        --evalue '!{params.prokka_evalue}' \
        --coverage !{params.prokka_coverage} \
        --cpus !{task.cpus} \
        --centre '!{params.centre}' \
        --mincontiglen !{params.min_contig_len} \
        !{locustag} \
        !{prodigal} \
        !{addgenes} \
        !{compliant} \
        !{proteins} \
        !{rawproduct} \
        !{cdsrnaolap} \
        !{addmrna} \
        !{norrna} \
        !{notrna} \
        !{rnammer} \
        !{rfam} \
        !{fasta} > prokka.stdout.txt 2> prokka.stderr.txt

    # Files passed to other modules
    ln -s results/!{sample}.faa
    ln -s results/!{sample}.ffn

    # Capture version
    prokka --version >> prokka.version.txt 2>&1
    '''

    stub:
    """
    mkdir annotation
    touch annotation/${sample}
    touch annotation/${sample}.ffn
    touch annotation/${sample}.faa
    """
}
