/**
 * Remove human contamination from sequencing reads for SRA submission.
 *
 * This subworkflow uses the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify
 * and remove human reads from sequencing data. It first initializes a human reference database
 * and then scrubs the input reads to ensure they meet SRA submission requirements.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords contamination, human, scrub, sra, sequencing, fastq
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation srahumanscrubber
 *
 * @modules srahumanscrubber_initdb, srahumanscrubber_scrub
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 *   - `special_meta`: Groovy Map with ID for downstream aggregation
 *   - `scrubbed`: Scrubbed FASTQ files with human reads removed
 *   - `scrubbed_extra`: Placeholder files for pipeline compatibility
 *   - `scrub_report`: Report of scrubbing statistics
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB  } from '../../modules/srahumanscrubber/scrub/main'
include { gather                  } from 'plugin/nf-bactopia'

workflow SRAHUMANSCRUBBER {
    take:
    reads: Channel<Record>

    main:
    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)

    emit:
    sample_outputs = SRAHUMANSCRUBBER_SCRUB.out
}
