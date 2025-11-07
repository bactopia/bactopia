process AMRFINDERPLUS_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(genes), path(proteins), path(gff)
    each path(db)

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    tuple val(meta), path("*.{log,err}")            , emit: logs, optional: true
    tuple val(meta), path(".command.begin")         , emit: nf_begin
    tuple val(meta), path(".command.err")           , emit: nf_err
    tuple val(meta), path(".command.log")           , emit: nf_log
    tuple val(meta), path(".command.out")           , emit: nf_out
    tuple val(meta), path(".command.run")           , emit: nf_run
    tuple val(meta), path(".command.sh")            , emit: nf_sh
    tuple val(meta), path(".command.trace")         , emit: nf_trace
    tuple val(meta), path("versions.yml")           , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // WF specific parameters
    def fna_is_compressed = genes.getName().endsWith(".gz") ? true : false
    def faa_is_compressed = proteins.getName().endsWith(".gz") ? true : false
    def gff_is_compressed = gff.getName().endsWith(".gz") ? true : false
    organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fna_name = genes.getName().replace(".gz", "")
    faa_name = proteins.getName().replace(".gz", "")
    gff_name = gff.getName().replace(".gz", "")
    annotation_format = gff_name.endsWith(".gff") ? "prokka" : "bakta"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$fna_is_compressed" == "true" ]; then
        gzip -c -d $genes > $fna_name
    fi

    if [ "$faa_is_compressed" == "true" ]; then
        gzip -c -d $proteins > $faa_name
    fi

    if [ "$gff_is_compressed" == "true" ]; then
        gzip -c -d $gff > $gff_name
    fi

    # Extract database
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        AMRFINDER_DB=\$(find database/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    else
        AMRFINDER_DB=\$(find $db/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    fi

    # Full AMRFinderPlus search combining results
    amrfinder \\
        --nucleotide $fna_name \\
        --protein $faa_name \\
        --gff $gff_name \\
        --annotation_format $annotation_format \\
        $organism_param \\
        $task.ext.args \\
        --database \$AMRFINDER_DB \\
        --threads $task.cpus \\
        --name $prefix > ${prefix}.tsv

    # Clean up
    DB_VERSION=\$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    rm -rf database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$DB_VERSION
    END_VERSIONS
    """
}
