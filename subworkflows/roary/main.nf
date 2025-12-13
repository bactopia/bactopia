/**
 * Build a pangenome from GFF3 annotations using Roary.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [Roary](https://github.com/sanger-pathogens/Roary).
 * Roary is a rapid pangenome pipeline that processes large numbers of annotated genomes to produce
 * gene presence/absence matrices and core-genome alignments. It is particularly optimized for
 * bacterial datasets and can handle hundreds of genomes efficiently.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation roary
 *
 * @modules roary as roary_module
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

include { ROARY as ROARY_MODULE } from '../../modules/roary/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow ROARY {
    take:
    gff : Channel<Tuple<Map, Set<Path>>>

    main:
    ROARY_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Set<Path>>> = ROARY_MODULE.out.aln
    csv: Channel<Tuple<Map, Set<Path>>> = ROARY_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ROARY_MODULE.out.supplemental,
        ROARY_MODULE.out.aln,
        ROARY_MODULE.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ROARY_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([ROARY_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ROARY_MODULE.out.versions])
}
