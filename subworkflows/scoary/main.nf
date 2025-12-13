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
 * @input tuple(meta, csv)
 * - `meta`: Groovy Map containing sample information
 * - `csv`: Gene presence/absence matrix from pan-genome analysis in CSV format
 *
 * @input traits
 * Trait file containing binary phenotypic characteristics for each isolate (optional)
 *
 * @output csv         Scoary GWAS results with statistical associations between genes and traits
 * @output results     Aggregated results channel containing all output files
 * @output logs        Execution logs from the Scoary analysis
 * @output nf_logs     Nextflow execution scripts and logs for debugging
 * @output versions    Software version information
 */
nextflow.preview.types = true

include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCOARY {
    take:
    csv: Channel<Tuple<Map, Set<Path>>>
    traits: Path?

    main:    
    SCOARY_MODULE(csv, traits)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Set<Path>>> = SCOARY_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([SCOARY_MODULE.out.csv])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SCOARY_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SCOARY_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SCOARY_MODULE.out.versions])
}
