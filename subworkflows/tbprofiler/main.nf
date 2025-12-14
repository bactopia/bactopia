/**
 * Profiling tool for Mycobacterium tuberculosis to detect resistance and strain type.
 *
 * This subworkflow performs comprehensive profiling of Mycobacterium tuberculosis
 * from sequencing reads using [TBProfiler](https://github.com/jodyphelan/TBProfiler).
 * The tool detects drug resistance mutations, determines lineage and strain type,
 * and provides detailed variant calling results. It combines individual sample
 * results with population-level analysis for surveillance and epidemiological studies.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords Mycobacterium, tuberculosis, drug resistance, lineage, variants
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation tbprofiler
 *
 * @modules tbprofiler_profile, tbprofiler_collate
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output csv          TBProfiler resistance and lineage results in CSV format
 * @output json         Detailed analysis results with variant information in JSON format
 * @output txt          Summary report of resistance mutations and lineage
 * @output bam          Aligned reads in BAM format for visualization
 * @output vcf          Variant calls in VCF format
 * @output merged_csv   Combined resistance and lineage results from all samples
 * @output variants_csv Collated variant information across all samples
 * @output variants_txt Summary of detected resistance mutations
 * @output itol         iTOL input file for phylogenetic tree annotation
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { TBPROFILER_PROFILE } from '../../modules/tbprofiler/profile/main'
include { TBPROFILER_COLLATE } from '../../modules/tbprofiler/collate/main'
include { flattenPaths       } from 'plugin/nf-bactopia'
include { gather             } from 'plugin/nf-bactopia'

workflow TBPROFILER {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>

    main:
    TBPROFILER_PROFILE(reads)
    TBPROFILER_COLLATE(gather(TBPROFILER_PROFILE.out.json, 'tbprofiler', 'json'))

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_PROFILE.out.csv
    json: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_PROFILE.out.json
    txt: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_PROFILE.out.txt
    bam: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_PROFILE.out.bam
    vcf: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_PROFILE.out.vcf
    merged_csv: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_COLLATE.out.csv
    variants_csv: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_COLLATE.out.variants_csv
    variants_txt: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_COLLATE.out.variants_txt
    itol: Channel<Tuple<Map, Set<Path>>> = TBPROFILER_COLLATE.out.itol

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.csv,
        TBPROFILER_PROFILE.out.json,
        TBPROFILER_PROFILE.out.txt,
        TBPROFILER_PROFILE.out.bam,
        TBPROFILER_PROFILE.out.vcf,
        TBPROFILER_COLLATE.out.csv,
        TBPROFILER_COLLATE.out.variants_csv,
        TBPROFILER_COLLATE.out.variants_txt,
        TBPROFILER_COLLATE.out.itol
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.logs,
        TBPROFILER_COLLATE.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.nf_logs,
        TBPROFILER_COLLATE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        TBPROFILER_PROFILE.out.versions,
        TBPROFILER_COLLATE.out.versions
    ])
}
