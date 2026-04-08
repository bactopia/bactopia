/**
 * Remove human contamination from sequencing reads for SRA submission.
 *
 * This subworkflow uses the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify
 * and remove human reads from sequencing data. It first initializes a human reference database
 * and then scrubs the input reads to ensure they meet SRA submission requirements.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords contamination, human, scrub, sra, sequencing, fastq
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation srahumanscrubber
 *
 * @modules srahumanscrubber_initdb, srahumanscrubber_scrub
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 *   - `special_meta`: Groovy Map with name for downstream aggregation
 *   - `scrubbed`: Scrubbed FASTQ files with human reads removed
 *   - `scrubbed_extra`: Placeholder files for pipeline compatibility
 *   - `scrub_report`: Report of scrubbing statistics
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB  } from '../../modules/srahumanscrubber/scrub/main'
include { gatherCsvtk             } from 'plugin/nf-bactopia'

workflow SRAHUMANSCRUBBER {
    take:
    reads: Channel<Record>

    main:
    ch_srahumanscrubber_initdb = SRAHUMANSCRUBBER_INITDB()
    ch_srahumanscrubber_scrub = SRAHUMANSCRUBBER_SCRUB(reads, ch_srahumanscrubber_initdb.map{ r -> r.db })

    emit:
    // Published outputs
    sample_outputs = ch_srahumanscrubber_scrub
    run_outputs = channel.empty()
}
