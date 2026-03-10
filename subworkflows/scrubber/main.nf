/**
 * Remove contaminant sequences from metagenomic data.
 *
 * This subworkflow removes human and other contaminant sequences from metagenomic reads using either
 * the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) or a Kraken2-based approach (k2scrubber)
 * with the HPRC human database. It provides flexible contamination removal with detailed reporting
 * and aggregates results across multiple samples.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, decontamination, human removal, read filtering
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic, aggregation
 * @citation kraken2, srahumanscrubber
 *
 * @subworkflows srahumanscrubber, k2scrubber
 * @modules csvtk_concat
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input use_srascrubber
 * Boolean flag to choose between SRA Human Scrubber (true) or k2scrubber (false) for decontamination.
 *
 * @output sample_outputs  Per-sample module record from the selected scrubbing tool
 * @output scrubbed        Clean metagenomic reads after contaminant removal (unified from SRA/K2)
 * @output scrubbed_extra  Additional cleaned reads from extended filtering
 * @output special_tsv     Contamination reports keyed by special_meta for downstream joining
 * @output run_outputs     Merged contamination reports across all samples from CSVTK_CONCAT
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { K2SCRUBBER       } from '../k2scrubber/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gather           } from 'plugin/nf-bactopia'

workflow SCRUBBER {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    use_srascrubber: Boolean

    main:
    ch_sample_outputs = channel.empty()
    ch_scrub_report = channel.empty()
    ch_special_report = channel.empty()
    ch_scrubbed = channel.empty()
    ch_scrubbed_extra = channel.empty()

    if (use_srascrubber) {
        SRAHUMANSCRUBBER(reads)
        ch_sample_outputs = SRAHUMANSCRUBBER.out.sample_outputs
        ch_scrub_report = SRAHUMANSCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.scrub_report) }
        ch_special_report = SRAHUMANSCRUBBER.out.sample_outputs.map { r -> tuple(r.special_meta, r.scrub_report) }
        ch_scrubbed = SRAHUMANSCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.scrubbed) }
        ch_scrubbed_extra = SRAHUMANSCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.scrubbed_extra) }
    } else {
        K2SCRUBBER(reads)
        ch_sample_outputs = K2SCRUBBER.out.sample_outputs
        ch_scrub_report = K2SCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.scrub_report) }
        ch_special_report = K2SCRUBBER.out.sample_outputs.map { r -> tuple(r.special_meta, r.scrub_report) }
        ch_scrubbed = K2SCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.unclassified) }
        ch_scrubbed_extra = K2SCRUBBER.out.sample_outputs.map { r -> tuple(r.meta, r.unclassified_extra) }
    }

    CSVTK_CONCAT(gather(ch_scrub_report, 'scrubber'), 'tsv', 'tsv')

    emit:
    sample_outputs = ch_sample_outputs
    scrubbed = ch_scrubbed
    scrubbed_extra = ch_scrubbed_extra
    special_tsv = ch_special_report
    run_outputs = CSVTK_CONCAT.out
}
