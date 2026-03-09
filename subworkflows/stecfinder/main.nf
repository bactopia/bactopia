/**
 * Identify and serotype Shiga toxin-producing E. coli (STEC) from assemblies.
 *
 * This subworkflow uses [STECFinder](https://github.com/LanLab/STECFinder) to identify
 * and serotype Shiga toxin-producing *E. coli* (STEC) strains using genomic cluster-specific
 * markers. It screens assemblies for virulence genes and serotype markers to classify
 * STEC isolates into their known serotypes.
 *
 * @status stable
 * @keywords escherichia coli, stec, serotype, virulence genes, shiga toxin
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation stecfinder
 *
 * @modules stecfinder, csvtk_concat
 *
 * @input tuple(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @output tsv         Per-sample TSV files containing STEC identification and serotyping results
 * @output merged_tsv  Consolidated TSV file containing STEC results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { STECFINDER as STECFINDER_MODULE } from '../../modules/stecfinder/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gather                          } from 'plugin/nf-bactopia'

workflow STECFINDER {
    take:
    seqs: Channel<Tuple<Map, Path, Path?, Path?, Path?, Path?>>

    main:
    STECFINDER_MODULE(seqs)
    CSVTK_CONCAT(gather(STECFINDER_MODULE.out, 'stecfinder', field: 'tsv'), 'tsv', 'tsv')
    emit:
    sample_outputs = STECFINDER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
