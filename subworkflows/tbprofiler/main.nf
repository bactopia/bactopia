/**
 * Profiling tool for Mycobacterium tuberculosis to detect resistance and strain type.
 *
 * This subworkflow performs comprehensive profiling of Mycobacterium tuberculosis
 * from sequencing reads using [TBProfiler](https://github.com/jodyphelan/TBProfiler).
 * The tool detects drug resistance mutations, determines lineage and strain type,
 * and provides detailed variant calling results. It combines individual sample
 * results with population-level analysis for surveillance and epidemiological studies.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords Mycobacterium, tuberculosis, drug resistance, lineage, variants
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation tbprofiler
 *
 * @modules tbprofiler_profile, tbprofiler_collate
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 * - `bam`: Aligned BAM file
 * - `csv`: Results in CSV format
 * - `json`: Compressed JSON results file
 * - `txt`: Results in text format
 * - `vcf`: Compressed VCF file with variants
 *
 * @output run_outputs
 * - `csv`: Main collated results in CSV format
 * - `variants_csv`: Collated variants in CSV format
 * - `variants_txt`: Collated variants in text format
 * - `itol`: iTOL formatted files for visualization
 */
nextflow.preview.types = true

include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'
include { gatherCsvtk             } from 'plugin/nf-bactopia'

workflow TBPROFILER {
    take:
    reads: Channel<Record>

    main:
    TBPROFILER_PROFILE(reads)
    TBPROFILER_COLLATE(gatherCsvtk(TBPROFILER_PROFILE.out, 'json', [name: 'tbprofiler']))

    emit:
    // Published outputs
    sample_outputs = TBPROFILER_PROFILE.out
    run_outputs = TBPROFILER_COLLATE.out
}
