/**
 * Download and provide pre-compiled datasets required by Bactopia.
 *
 * This subworkflow wraps the DATASETS module and extracts individual database
 * paths as separate channel emissions for downstream consumption.
 *
 * @status stable
 * @keywords download, database, setup, amr, mlst, minhash, sourmash, gtdb
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download
 *
 * @output amrfinderplus_db  Path to the AMRFinderPlus database tarball
 * @output mlst_db           Path to the PubMLST database tarball
 * @output mash_db           Path to the Mash RefSeq sketch
 * @output sourmash_db       Path to the Sourmash GTDB signatures
 */
// bactopia-lint: ignore S003,S005,S010
nextflow.preview.types = true

include { DATASETS as DATASETS_MODULE } from '../../../modules/bactopia/datasets/main'

workflow DATASETS {

    main:
    DATASETS_MODULE()

    emit:
    amrfinderplus_db = DATASETS_MODULE.out.map { r -> r.amrfinderplus_db }
    mlst_db          = DATASETS_MODULE.out.map { r -> r.mlst_db }
    mash_db          = DATASETS_MODULE.out.map { r -> r.mash_db }
    sourmash_db      = DATASETS_MODULE.out.map { r -> r.sourmash_db }
}
