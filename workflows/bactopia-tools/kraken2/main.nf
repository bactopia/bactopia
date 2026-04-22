#!/usr/bin/env nextflow
/**
 * Taxonomic classification of metagenomic sequence reads.
 *
 * This Bactopia Tool uses [Kraken2](https://github.com/DerrickWood/kraken2) to assign taxonomic
 * classifications to metagenomic sequence reads. It creates reports of taxonomic assignments
 * and generates interactive Krona plots for visualization of classification results.
 *
 * @status stable
 * @keywords taxonomy, metagenomics, classification, fastq, kraken2, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,classification,visualization
 * @citation kraken2
 *
 * @subworkflows utils_bactopia-tools, kraken2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input kraken2_db
 * Path to Kraken2 database for taxonomic classification
 *
 * @section Per-Sample Results
 * @publish *.report          Kraken2 taxonomic classification report
 * @publish *.kraken2         Kraken2 classification results
 * @publish *.html            Krona interactive HTML report
 *
 * @section Merged Results
 * @publish merged-report     Merged Kraken2 reports from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    kraken2_db : Value<Path>
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { KRAKEN2             } from '../../../subworkflows/kraken2/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_kraken2 = KRAKEN2(ch_bactopiatool.reads, params.kraken2_db)

    publish:
    // Per-sample
    sample_outputs = ch_kraken2.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_kraken2.sample_outputs)
    // Run-level
    run_outputs = ch_kraken2.run_outputs
    run_nf_logs = collectNextflowLogs(ch_kraken2.run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
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
