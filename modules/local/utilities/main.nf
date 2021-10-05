nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

process MERGE_TABLES {
    /* Merge multiple TSVs into a single TSV */
    beforeScript 'ulimit -Ss unlimited'
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(tables).collect()

    output:
    tuple val(meta), path("${meta.id}.${extension}"), emit: table
    path "*.version.txt", emit: version

    shell:
    extension = meta.containsKey('extension') ? meta.extension : 'tsv'
    '''
    head -n1 !{tables}[0] > !{meta.id}.!{extension}
    tail -n+2 !{tables} >> !{meta.id}.!{extension}

    echo \$(head --version 2>&1) | grep "^head" | sed 's/head (GNU coreutils) //' > head.version.txt
    echo \$(tail --version 2>&1) | grep "^head" | sed 's/tail (GNU coreutils) //' > tail.version.txt
    '''
}
