nextflow.preview.types = true

process DATASETS {
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    amrfinderplus_db = file(("${bactopia_version}/amrfinderplus.tar.gz"))
    mlst_db          = file(("${bactopia_version}/mlst.tar.gz"))
    mash_db          = file(("mash-refseq88.k21.msh.xz"))
    sourmash_db      = file(("gtdb-rs207.genomic-reps.dna.k31.lca.json.gz"))

    script:
    bactopia_version = task.ext.amrfinder_url.replace("https://datasets.bactopia.com/datasets/", "").replace("/amrfinderplus.tar.gz", "")
    """
    mkdir -p ${bactopia_version}
    wget -O ${bactopia_version}/amrfinderplus.tar.gz ${task.ext.amrfinder_url}
    wget -O ${bactopia_version}/mlst.tar.gz ${task.ext.mlst_url}
    wget -O gtdb-rs207.genomic-reps.dna.k31.lca.json.gz ${task.ext.sourmash_url}
    wget -O mash-refseq88.k21.msh.xz ${task.ext.mash_url}
    """
}
