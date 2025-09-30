process DATASETS {
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path("${bactopia_version}/amrfinderplus.tar.gz")   , emit: amrfinderplus_db
    path("${bactopia_version}/mlst.tar.gz")            , emit: mlst_db
    path("mash-refseq88.k21.msh.xz")                   , emit: mash_db
    path("gtdb-rs207.genomic-reps.dna.k31.lca.json.gz"), emit: sourmash_db
    path("*.{log,err}")  , emit: logs, optional: true
    path ".command.begin", emit: nf_begin
    path ".command.err"  , emit: nf_err
    path ".command.log"  , emit: nf_log
    path ".command.out"  , emit: nf_out
    path ".command.run"  , emit: nf_run
    path ".command.sh"   , emit: nf_sh
    path ".command.trace", emit: nf_trace
    path "versions.yml"  , emit: versions

    script:
    bactopia_version = task.ext.amrfinder_url.replace("https://datasets.bactopia.com/datasets/","").replace("/amrfinderplus.tar.gz","")
    """
    mkdir -p ${bactopia_version}
    wget -O ${bactopia_version}/amrfinderplus.tar.gz ${task.ext.amrfinder_url}
    wget -O ${bactopia_version}/mlst.tar.gz ${task.ext.mlst_url}
    wget -O gtdb-rs207.genomic-reps.dna.k31.lca.json.gz ${task.ext.sourmash_url}
    wget -O mash-refseq88.k21.msh.xz ${task.ext.mash_url}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version | head -n 1 | sed 's/^GNU Wget //'))
    END_VERSIONS
    """
}
