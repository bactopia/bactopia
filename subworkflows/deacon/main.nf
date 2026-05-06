/**
 * Remove host reads from sequencing data using deacon.
 *
 * This subworkflow uses [deacon](https://github.com/bede/deacon) to identify and remove
 * host reads from FASTQ files using SIMD-accelerated minimizer comparison against a
 * pre-built reference index (default: panhuman-1). It optionally downloads the index
 * if not already available.
 *
 * @status stable
 * @keywords host, contamination, decontamination, depletion, filtering, minimizer, reads, deacon
 * @tags complexity:moderate input-type:single output-type:single features:database-dependent,conditional-logic,resource-download
 * @citation deacon
 *
 * @modules deacon_fetch, deacon_filter
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to deacon minimizer index file (.idx) (ignored if download_deacon is true)
 *
 * @input download_deacon
 * Boolean flag to download the index instead of using the provided path
 *
 * @output sample_outputs
 * - `special_meta`: Simplified metadata record for downstream report joining
 * - `r1?`: Filtered paired-end forward reads
 * - `r2?`: Filtered paired-end reverse reads
 * - `se?`: Filtered single-end reads
 * - `lr?`: Filtered long reads
 * - `scrub_report`: Read filtering statistics report
 *
 * @output run_outputs
 */
nextflow.enable.types = true

include { DEACON_FETCH                   } from '../../modules/deacon/fetch/main'
include { DEACON_FILTER as DEACON_MODULE } from '../../modules/deacon/filter/main'

workflow DEACON {
    take:
    reads: Channel<Record>
    database: Path?
    download_deacon: Boolean

    main:
    ch_deacon = channel.empty()
    if (download_deacon) {
        ch_deacon_fetch = DEACON_FETCH()
        ch_deacon = DEACON_MODULE(reads, ch_deacon_fetch.map { r -> r.db })
    } else {
        ch_deacon = DEACON_MODULE(reads, database)
    }

    emit:
    // Published outputs
    sample_outputs = ch_deacon
    run_outputs = channel.empty()
}
