/**
 * Gene Allele Mutation Microbial Assessment.
 *
 * This subworkflow performs rapid identification, classification, and annotation of
 * translated gene matches from sequencing data using [GAMMA](https://github.com/rastanton/GAMMA).
 * The tool screens input sequences against a protein database to identify gene
 * variants, mutations, and allele types, providing detailed annotation and classification.
 *
 * @status stable
 * @keywords gene, allele, mutation, variant, antimicrobial resistance
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation gamma
 *
 * @modules csvtk_concat, gamma
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for gene allele identification
 *
 * @input db
 * Protein database file for sequence comparison (required)
 *
 * @output sample_outputs
 * - `gamma`: Main GAMMA output file containing annotated gene matches
 * - `psl`: Raw alignment details in PSL format
 * - `gff`: Gene matches in GFF3 format
 * - `fasta`: Extracted nucleotide sequences of the matched genes
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { GAMMA as GAMMA_MODULE } from '../../modules/gamma/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { gather                } from 'plugin/nf-bactopia'

workflow GAMMA {
    take:
    assembly: Channel<Record>
    db: Path

    main:
    GAMMA_MODULE(assembly, db)
    CSVTK_CONCAT(gather(GAMMA_MODULE.out, 'gamma', [name: 'gamma']), 'tsv', 'tsv')

    emit:
    sample_outputs = GAMMA_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
