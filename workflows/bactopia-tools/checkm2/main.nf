#!/usr/bin/env nextflow
/**
 * Machine learning-based assessment of microbial genome assembly quality.
 *
 * This Bactopia Tool uses [CheckM2](https://github.com/chklovski/CheckM2) to assess the quality
 * of microbial genomes recovered from isolates, single cells, and metagenomes using
 * advanced machine learning approaches.
 *
 * @status stable
 * @keywords assembly quality, microbial genomes, machine learning, completeness, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,assembly-qa,machine-learning
 * @citation checkm2
 *
 * @subworkflows bactopiatool_init, checkm2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input checkm2_db
 * Path to CheckM2 database
 *
 * @input download_checkm2
 * Download CheckM2 database if not found locally
 *
 * @section Quality Assessment
 * @publish diamond_output/**    Directory with intermediate results from CheckM2 processing
 * @publish protein_files/**     Directory containing protein files used for analysis
 * @publish quality_report.tsv   Output file with completeness statistics
 *
 * @section Merged Results
 * @publish checkm2.tsv          Merged TSV file with CheckM2 results from all samples
 *
 * @section Execution Logs
 * @publish logs/**              Tool execution logs
 * @publish logs/nf-*            Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml         Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    checkm2_db        : Path
    download_checkm2  : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { CHECKM2           } from '../../../subworkflows/checkm2/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    CHECKM2(
        BACTOPIATOOL_INIT.out.assembly,
        params.checkm2_db,
        params.download_checkm2
    )

    ch_sample_nf_logs = CHECKM2.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = CHECKM2.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = CHECKM2.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = CHECKM2.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results      >> "${r.meta.output_dir}/"
            r.supplemental >> "${r.meta.output_dir}/"
            r.logs         >> "${r.meta.logs_dir}/"
            r.versions     >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
