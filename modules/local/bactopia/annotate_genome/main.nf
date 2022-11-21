nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'annotate_genome')

process ANNOTATE_GENOME {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${meta.id}"
    label 'annotate_genome'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(genome_size), path(fasta), path(total_contigs), path(prokka_proteins), path(prodigal_tf)

    output:
    tuple val(meta), path("${meta.id}.{ffn,ffn.gz}"), path("${meta.id}.{faa,faa.gz}"), emit: annotations
    path "results/*", emit: results
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

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
    locustag = "--locustag ${meta.id}"
    renamed = false
    // Contig ID must <= 37 characters
    if ("gnl|${params.centre}|${meta.id}_${contig_count}".length() > 37) {
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
        echo "Original sample name (!{meta.id}) not used due to creating a contig ID >37 characters"
    fi

    if [[ "!{params.skip_compression}" == "false" ]]; then
        gunzip -c !{fasta} > !{meta.id}.fna
    fi

    prokka --outdir results \
        --force \
        --prefix '!{meta.id}' \
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
        !{meta.id}.fna
    mv results/!{meta.id}.err ./
    mv results/!{meta.id}.log ./

    if [[ "!{params.skip_compression}" == "false" ]]; then
        find results/ -type f | \
            grep -v -E "\\.err$|\\.log$|\\.tsv$|\\.txt$"| \
            xargs -I {} pigz -n --best -p !{task.cpus} {}

        # Files passed to other modules
        ln -s results/!{meta.id}.faa.gz !{meta.id}.faa.gz
        ln -s results/!{meta.id}.ffn.gz !{meta.id}.ffn.gz
    else
        # Files passed to other modules
        ln -s results/!{meta.id}.faa !{meta.id}.faa
        ln -s results/!{meta.id}.ffn !{meta.id}.ffn
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        prokka:  $(echo $(prokka --version 2>&1) | sed 's/prokka //')
    END_VERSIONS
    '''
}
