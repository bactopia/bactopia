/**
 * Remove contaminant sequences from metagenomic data.
 *
 * This subworkflow removes human and other contaminant sequences from metagenomic reads using
 * [deacon](https://github.com/bede/deacon) (default), [nohuman](https://github.com/mbhall88/nohuman),
 * or the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber). It provides flexible
 * contamination removal with detailed reporting and aggregates results across multiple samples.
 *
 * @status stable
 * @keywords metagenomics, decontamination, human removal, read filtering
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic,aggregation
 * @citation deacon, kraken2, srahumanscrubber
 *
 * @modules csvtk_concat
 * @subworkflows deacon, srahumanscrubber, nohuman
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input use_srascrubber
 * Boolean flag to use SRA Human Scrubber for decontamination
 *
 * @input use_nohuman
 * Boolean flag to use nohuman for decontamination
 *
 * @input nohuman_db
 * Path to nohuman database directory or tarball (used when use_nohuman is true)
 *
 * @input download_nohuman
 * Boolean flag to download the nohuman database instead of using the provided path
 *
 * @input nohuman_save_as_tarball
 * Boolean flag to save downloaded nohuman database as tarball
 *
 * @input deacon_db
 * Path to deacon minimizer index file (.idx) (used when deacon is selected)
 *
 * @input download_deacon
 * Boolean flag to download the deacon index instead of using the provided path
 *
 * @output sample_outputs
 * - `special_meta`: Simplified metadata record for downstream report joining
 * - `r1?`: Scrubbed paired-end forward reads
 * - `r2?`: Scrubbed paired-end reverse reads
 * - `se?`: Scrubbed single-end reads
 * - `lr?`: Scrubbed long reads
 * - `scrub_report`: Contamination removal statistics report
 *
 * @output run_outputs
 * - `csv`: Aggregated contamination reports across all samples
 *
 * @output scrubbed
 * - `r1?`: Scrubbed paired-end forward reads
 * - `r2?`: Scrubbed paired-end reverse reads
 * - `se?`: Scrubbed single-end reads
 * - `lr?`: Scrubbed long reads
 *
 * @output scrubbed_extra
 * - `r1?`: Scrubbed paired-end forward reads
 * - `r2?`: Scrubbed paired-end reverse reads
 * - `se?`: Scrubbed single-end reads
 * - `lr?`: Scrubbed long reads
 * - `fna?`: Assembly file (passed through)
 *
 * @output special_tsv
 * - `special_meta`: Simplified metadata record for downstream report joining
 * - `scrub_report`: Contamination removal statistics report
 */
nextflow.enable.types = true

include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { NOHUMAN          } from '../nohuman/main'
include { DEACON           } from '../deacon/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { gatherCsvtk      } from 'plugin/nf-bactopia'
include { filterWithData   } from 'plugin/nf-bactopia'

workflow SCRUBBER {
    take:
    reads: Channel<Record>
    use_srascrubber: Boolean
    use_nohuman: Boolean
    nohuman_db: Path?
    download_nohuman: Boolean
    nohuman_save_as_tarball: Boolean
    deacon_db: Path?
    download_deacon: Boolean

    main:
    ch_sample_outputs = channel.empty()
    ch_special_report = channel.empty()
    ch_reads = filterWithData(reads, ['r1', 'r2', 'se', 'lr'])

    if (use_srascrubber) {
        ch_srahumanscrubber = SRAHUMANSCRUBBER(ch_reads)
        ch_sample_outputs = ch_srahumanscrubber.sample_outputs
        ch_special_report = ch_srahumanscrubber.sample_outputs.map { r -> record(special_meta: r.special_meta, scrub_report: r.scrub_report) }
    } else if (use_nohuman) {
        ch_nohuman = NOHUMAN(ch_reads, nohuman_db, download_nohuman, nohuman_save_as_tarball)
        ch_sample_outputs = ch_nohuman.sample_outputs
        ch_special_report = ch_nohuman.sample_outputs.map { r -> record(special_meta: r.special_meta, scrub_report: r.scrub_report) }
    } else {
        ch_deacon = DEACON(ch_reads, deacon_db, download_deacon)
        ch_sample_outputs = ch_deacon.sample_outputs
        ch_special_report = ch_deacon.sample_outputs.map { r -> record(special_meta: r.special_meta, scrub_report: r.scrub_report) }
    }

    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_sample_outputs, 'scrub_report', [name: 'scrubber']), 'tsv', 'tsv')

    emit:
    // Downstream inputs
    scrubbed = filterWithData(ch_sample_outputs, ['r1', 'r2', 'se', 'lr'])
    scrubbed_extra = filterWithData(ch_sample_outputs, ['r1', 'r2', 'se', 'lr', 'fna'])
    special_tsv = ch_special_report
    // Published outputs
    sample_outputs = ch_sample_outputs
    run_outputs = ch_csvtk_concat
}
