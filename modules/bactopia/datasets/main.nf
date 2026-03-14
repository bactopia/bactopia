/**
 * Download pre-compiled datasets required by Bactopia.
 *
 * Fetches the core datasets (AMR, MLST, Mash, Sourmash) hosted by the Bactopia project.
 * These are used to populate the local cache for offline use.
 *
 * @status stable
 * @keywords download, database, setup, amr, mlst, minhash, sourmash, gtdb
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download
 * @citation bactopia, amrfinderplus, mlst, mash, sourmash
 *
 * @note Internet Required
 * This process requires an active internet connection to fetch files from `datasets.bactopia.com`.
 *
 * @output record(amrfinderplus_db, mlst_db, mash_db, sourmash_db)
 * - `amrfinderplus_db`: A compressed tarball of the [AMRFinderPlus](https://github.com/ncbi/amr) database
 * - `mlst_db`: A compressed tarball of the [PubMLST](https://pubmlst.org/) schemes
 * - `mash_db`: Pre-computed [Mash](https://github.com/marbl/Mash) sketches (RefSeq)
 * - `sourmash_db`: Pre-computed [Sourmash](https://github.com/sourmash-bio/sourmash) signatures (GTDB)
 */
nextflow.preview.types = true

process DATASETS {
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    amrfinderplus_db = file(("${bactopia_version}/amrfinderplus.tar.gz"))
    mlst_db          = file(("mlst.tar.gz"))
    mash_db          = file(("mash-refseq88.k21.msh.xz"))
    sourmash_db      = file(("gtdb-rs207.genomic-reps.dna.k31.lca.json.gz"))

    script:
    bactopia_version = task.ext.amrfinder_url.replace("https://datasets.bactopia.com/datasets/", "").replace("/amrfinderplus.tar.gz", "")
    """
    mkdir -p ${bactopia_version}
    wget -O ${bactopia_version}/amrfinderplus.tar.gz ${task.ext.amrfinder_url}
    wget -O mlst.tar.gz ${task.ext.mlst_url}
    wget -O gtdb-rs207.genomic-reps.dna.k31.lca.json.gz ${task.ext.sourmash_url}
    wget -O mash-refseq88.k21.msh.xz ${task.ext.mash_url}
    """
}
