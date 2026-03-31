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
 * @input record(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @output sample_outputs
 * - `tsv`: TSV file with STEC gene markers results
 *
 * @output run_outputs
 * - `csv`: Merged STEC results from all samples
 */
nextflow.preview.types = true

include { STECFINDER as STECFINDER_MODULE } from '../../modules/stecfinder/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                     } from 'plugin/nf-bactopia'

workflow STECFINDER {
    take:
    seqs: Channel<Record>

    main:
    STECFINDER_MODULE(seqs)
    CSVTK_CONCAT(gatherCsvtk(STECFINDER_MODULE.out, 'tsv', [name: 'stecfinder']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = STECFINDER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
