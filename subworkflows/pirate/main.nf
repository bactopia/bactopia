/**
 * Build a pangenome from GFF3 annotations using PIRATE.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [PIRATE](http://github.com/SionBayliss/PIRATE).
 * PIRATE is a scalable pangenome toolbox that clusters orthologous genes at multiple identity thresholds.
 * It is particularly useful for highly diverse datasets as it can handle divergent gene families
 * and provides flexible clustering options for different analytical needs.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation pirate
 *
 * @modules pirate as pirate_module
 *
 * @input tuple(meta, gff)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `gff`: Set of GFF3 annotation files representing the genomic annotations for each sample
 *
 * @output aln      Core-genome alignment file containing genes present across all input genomes
 * @output csv      Gene presence/absence matrix in CSV format showing which genes are present in each genome
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PIRATE {
    take:
    gff: Channel<Tuple<Map, Set<Path>>>

    main:
    PIRATE_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Set<Path>>> = PIRATE_MODULE.out.aln
    csv: Channel<Tuple<Map, Set<Path>>> = PIRATE_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PIRATE_MODULE.out.supplemental,
        PIRATE_MODULE.out.aln,
        PIRATE_MODULE.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PIRATE_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PIRATE_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PIRATE_MODULE.out.versions])
}
