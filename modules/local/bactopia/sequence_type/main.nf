nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'sequence_type')

process SEQUENCE_TYPE {
    /* Determine MLST types using BLAST */
    tag "${meta.id} - ${schema}"
    label "max_cpus_1"
    label "sequence_type"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir: schema) }

    input:
    tuple val(meta), path(assembly), path(fq)
    each path(dataset)

    output:
    path "results/*", emit: results
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    dataset_tarball = dataset.getName()
    schema = dataset_tarball.replace('.tar.gz', '')
    blast_opts = params.skip_compression ? "" : "--compressed"
    '''
    tar -xzvf !{dataset_tarball}

    # Run BLAST
    OUTDIR="results/!{schema}"
    mkdir -p ${OUTDIR}/blast
    mlst-blast.py !{assembly} !{schema}/blastdb ${OUTDIR}/blast/!{meta.id}-blast.json --cpu !{task.cpus} !{blast_opts}

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        blastn: $(echo $(blastn -version 2>&1) | sed 's/^.*blastn: //;s/ .*$//')
    END_VERSIONS
    '''
}
