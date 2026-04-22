/**
 * Pan-genome wide association studies.
 *
 * This subworkflow performs genome-wide association studies (GWAS) on pan-genome data
 * using [Scoary](https://github.com/AdmiralenOla/Scoary). The tool identifies genes
 * associated with binary traits such as pathogenicity, host specificity, or antibiotic
 * resistance. It calculates statistical associations between gene presence/absence
 * and phenotypic traits across multiple bacterial isolates.
 *
 * @status stable
 * @keywords GWAS, association, pan-genome, traits, statistical
 * @tags complexity:simple input-type:multiple output-type:single features:conditional-input
 * @citation scoary
 *
 * @modules scoary
 *
 * @input record(meta, csv)
 * - `meta`: Groovy Record containing sample information
 * - `csv`: Gene presence/absence matrix from pan-genome analysis in CSV format
 *
 * @input traits
 * Trait file containing binary phenotypic characteristics for each isolate (optional)
 *
 * @output sample_outputs
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'

workflow SCOARY {
    take:
    csv: Channel<Record>
    traits: Path?

    main:
    ch_scoary = SCOARY_MODULE(csv, traits)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_scoary
}
