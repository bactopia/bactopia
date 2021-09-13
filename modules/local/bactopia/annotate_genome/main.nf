nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "annotate_genome"

process ANNOTATE_GENOME {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), path(genome_size), path(fasta), path(total_contigs)
    path prokka_proteins
    path prodigal_tf

    output:
    tuple val(sample), path("${sample}.{ffn,ffn.gz}"), path("${sample}.{faa,faa.gz}"), emit: annotations
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
        proteins_name = prokka_proteins.getName()
        if(proteins_name.contains("-")) {
            genus = proteins_name.split('-')[0].capitalize()
            species = proteins_name.split('-')[1]
        } else {
            genus = proteins_name.capitalize()
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

    if [[ "!{params.skip_compression}" == "false" ]]; then
        gunzip -c !{sample}.fna.gz > !{sample}.fna
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
        !{sample}.fna > prokka.stdout.txt 2> prokka.stderr.txt

    if [[ "!{params.skip_compression}" == "false" ]]; then
        find results/ -type f | \
            grep -v -E "\\.err$|\\.log$|\\.tsv$|\\.txt$"| \
            xargs -I {} pigz -n --best -p !{task.cpus} {}

        # Files passed to other modules
        ln -s results/!{sample}.faa.gz !{sample}.faa.gz
        ln -s results/!{sample}.ffn.gz !{sample}.ffn.gz
    else
        # Files passed to other modules
        ln -s results/!{sample}.faa !{sample}.faa
        ln -s results/!{sample}.ffn !{sample}.ffn
    fi

    # Capture version
    prokka --version > prokka.version.txt 2>&1
    '''

    stub:
    """
    mkdir annotation
    touch annotation/${sample}
    touch annotation/${sample}.ffn
    touch annotation/${sample}.faa
    """
}
