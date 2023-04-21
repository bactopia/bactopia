// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'prokka')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::blast=2.12.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BLAST {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    tuple val(meta), path(blastdb, stageAs: 'blastdb/*')
    path fasta

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml"                   , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def cat_cmd = fasta.getName().endsWith(".gz") ? "zcat" : "cat"
    def blast_cmd = params.blast_type ==  "proteins" ? "tblastn" : "blastn"
    def added_opts = params.blast_type == "primers" ? '-task blastn -dust no -word_size 7' : ''
    """
    mkdir temp_out
    ${cat_cmd} ${query} | \
        sed -e 's/<[^>]*>//g' |
        parallel --gnu --plain -j ${task.cpus} --recstart '>' -N 1 --pipe \
            ${blast_cmd} -db blastdb/${meta.id} \
                ${options.args} ${added_opts} \
                -query - \
                -out temp_out/\${name}_{#}.txt
    ls temp_out | head -n 1 | xargs -I {} head -n 1 temp_out/{} > ${prefix}.txt
    ls temp_out | xargs -I {} tail -n +2 temp_out/{} >> ${prefix}.txt
    rm -rf temp_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(echo $(blastn -version 2>&1) | sed 's/^.*blastn: //;s/ .*\$//')
        parallel: \$(echo $(parallel --version 2>&1) | sed 's/^GNU parallel //;s/ .*\$//')
        pigz: \$(echo $(pigz --version 2>&1) | sed 's/pigz //')
        tblastn: \$(echo $(tblastn -version 2>&1) | sed 's/^.*tblastn: //;s/ .*\$//') 
    END_VERSIONS
    """
}
