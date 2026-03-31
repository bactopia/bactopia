/**
 * Remove contaminant sequences from metagenomic data.
 *
 * This subworkflow removes human and other contaminant sequences from metagenomic reads using either
 * the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) or [nohuman](https://github.com/mbhall88/nohuman)
 * with the HPRC human database. It provides flexible contamination removal with detailed reporting
 * and aggregates results across multiple samples.
 *
 * @status stable
 * @keywords metagenomics, decontamination, human removal, read filtering
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic,aggregation
 * @citation nohuman, srahumanscrubber
 *
 * @subworkflows srahumanscrubber, nohuman
 * @modules csvtk_concat
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input use_srascrubber
 * Boolean flag to choose between SRA Human Scrubber (true) or nohuman (false) for decontamination.
 *
 * @input nohuman_db
 * Path to nohuman database directory or tarball (used when use_srascrubber is false)
 *
 * @input download_nohuman
 * Boolean flag to download the nohuman database instead of using the provided path
 *
 * @input nohuman_save_as_tarball
 * Boolean flag to save downloaded nohuman database as tarball
 *
 * @output sample_outputs  Per-sample module record from the selected scrubbing tool
 * @output scrubbed        Clean metagenomic reads after contaminant removal
 * @output scrubbed_extra  Additional cleaned reads for extended filtering
 * @output special_tsv     Contamination reports keyed by special_meta for downstream joining
 * @output run_outputs     Merged contamination reports across all samples from CSVTK_CONCAT
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { NOHUMAN          } from '../nohuman/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gatherCsvtk      } from 'plugin/nf-bactopia'
include { filterWithData   } from 'plugin/nf-bactopia'

workflow SCRUBBER {
    take:
    reads: Channel<Record>
    use_srascrubber: Boolean
    nohuman_db: Path?
    download_nohuman: Boolean
    nohuman_save_as_tarball: Boolean

    main:
    ch_sample_outputs = channel.empty()
    ch_special_report = channel.empty()

    if (use_srascrubber) {
        SRAHUMANSCRUBBER(reads)
        ch_sample_outputs = SRAHUMANSCRUBBER.out.sample_outputs
        ch_special_report = SRAHUMANSCRUBBER.out.sample_outputs.map { r -> record(special_meta: r.special_meta, scrub_report: r.scrub_report) }
    } else {
        NOHUMAN(reads, nohuman_db, download_nohuman, nohuman_save_as_tarball)
        ch_sample_outputs = NOHUMAN.out.sample_outputs
        ch_special_report = NOHUMAN.out.sample_outputs.map { r -> record(special_meta: r.special_meta, scrub_report: r.scrub_report) }
    }

    CSVTK_CONCAT(gatherCsvtk(ch_sample_outputs, 'scrub_report', [name: 'scrubber']), 'tsv', 'tsv')

    emit:
    // Downstream inputs
    scrubbed = filterWithData(ch_sample_outputs, ['r1', 'r2', 'se', 'lr'])
    scrubbed_extra = filterWithData(ch_sample_outputs, ['r1', 'r2', 'se', 'lr', 'fna'])
    special_tsv = ch_special_report
    // Published outputs
    sample_outputs = ch_sample_outputs
    run_outputs = CSVTK_CONCAT.out
}
