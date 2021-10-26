// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process HICAP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::hicap=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hicap:1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hicap:1.0.3--py_0"
    }

    input:
    tuple val(meta), path(fasta)
    path database_dir
    path model_fp

    output:
    tuple val(meta), path("*.gbk"), emit: gbk, optional: true
    tuple val(meta), path("*.svg"), emit: svg, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    hicap \\
        --query_fp $fasta_name \\
        $database_args \\
        $model_args \\
        $options.args \\
        --threads $task.cpus \\
        --debug \\
        -o ./

    if [ ! -f ${meta.id}.tsv ]; then
        echo "isolate<TAB>predicted_serotype<TAB>attributes<TAB>genes_identified<TAB>locus_location<TAB>region_I_genes<TAB>region_II_genes<TAB>region_III_genes<TAB>IS1016_hits" | sed 's/<TAB>/\t/g' > ${meta.id}.tsv
        echo "${meta.id}<TAB>cap_not_found<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-" | sed 's/<TAB>/\t/g' >> ${meta.id}.tsv
    else
        sed -i 's/#isolate/isolate/' ${meta.id}.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    hicap:
        hicap: \$( echo \$( hicap --version 2>&1 ) | sed 's/^.*hicap //' )
    END_VERSIONS
    """
}
