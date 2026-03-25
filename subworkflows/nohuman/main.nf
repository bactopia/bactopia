/**
 * Remove human reads from sequencing data using nohuman.
 *
 * This subworkflow uses [nohuman](https://github.com/mbhall88/nohuman) to identify and remove
 * human reads from FASTQ files using a Kraken2 database built from Human Pangenome Reference
 * Consortium (HPRC) genomes. It optionally downloads the database if not already available.
 *
 * @status stable
 * @keywords human, contamination, decontamination, scrubbing, reads, nohuman, kraken2
 * @tags complexity:moderate input-type:single output-type:single features:database-dependent,conditional-logic,resource-download
 * @citation nohuman
 *
 * @modules nohuman_download, nohuman_run
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end forward)
 * - `r2`: Illumina R2 reads (paired-end reverse)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to nohuman database directory or tarball (ignored if download_nohuman is true)
 *
 * @input download_nohuman
 * Boolean flag to download the database instead of using the provided path
 *
 * @output sample_outputs
 * - `scrubbed`: FASTQ files with human reads removed
 * - `scrub_report`: Kraken2 classification report (optional)
 */
nextflow.preview.types = true

include { NOHUMAN_DOWNLOAD                } from '../../modules/nohuman/download/main'
include { NOHUMAN_RUN as NOHUMAN_MODULE   } from '../../modules/nohuman/run/main'

workflow NOHUMAN {
    take:
    reads: Channel<Record>
    database: Path
    download_nohuman: Boolean

    main:
    if (download_nohuman) {
        NOHUMAN_DOWNLOAD()
        NOHUMAN_MODULE(reads, NOHUMAN_DOWNLOAD.out.db)
    } else {
        NOHUMAN_MODULE(reads, database)
    }

    emit:
    sample_outputs = NOHUMAN_MODULE.out
}
