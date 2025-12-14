/**
 * Predict serotypes of Shigella from reads or assemblies.
 *
 * This subworkflow uses [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper) to predict
 * serotypes of *Shigella* strains from either Illumina/Nanopore reads or assembled genomes.
 * It analyzes antigen-encoding genes to determine the serotype classification of each isolate.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords shigella, serotype, typing, prediction, antigen genes
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation shigatyper
 *
 * @modules shigatyper, csvtk_concat
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output tsv         Per-sample TSV files containing Shigella serotype predictions
 * @output hits        Per-sample TSV files containing detailed gene hit information
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SHIGATYPER as SHIGATYPER_MODULE } from '../../modules/shigatyper/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { flattenPaths                    } from 'plugin/nf-bactopia'
include { gather                          } from 'plugin/nf-bactopia'

workflow SHIGATYPER {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>

    main:
    SHIGATYPER_MODULE(reads)
    CSVTK_CONCAT(gather(SHIGATYPER_MODULE.out.tsv, 'shigatyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SHIGATYPER_MODULE.out.tsv
    hits: Channel<Tuple<Map, Set<Path>>> = SHIGATYPER_MODULE.out.hits
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.tsv,
        SHIGATYPER_MODULE.out.hits,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGATYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
