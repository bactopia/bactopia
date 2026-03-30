#!/usr/bin/env nextflow
/**
 * Estimate taxonomic abundance of metagenomic samples.
 *
 * This Bactopia Tool uses [Bracken](https://github.com/jenniferlu717/Bracken) to estimate
 * taxonomic abundance from Kraken2 results. It also runs [Kraken2](https://ccb.jhu.edu/software/kraken2/)
 * for taxonomic classification and generates [Krona](https://github.com/marbl/Krona) interactive charts.
 *
 * @status stable
 * @keywords metagenomics, classification, abundance, kraken2, bracken, krona, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,abundance-profiling,database-dependent
 * @citation bracken, kraken2, krona
 *
 * @subworkflows bactopiatool_init, bracken
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input kraken2_db
 * Path to Kraken2 database for taxonomic classification
 *
 * @section Classification Results
 * @publish *.kraken2.report.txt   Kraken2 classification report
 * @publish *.bracken.report.txt   Bracken abundance estimates
 * @publish *.krona.html          Krona interactive visualization
 *
 * @section Summary Reports
 * @publish bracken-summary.tsv    Summary of classification results across all samples
 * @publish bracken-matrix.tsv     Abundance matrix for downstream analysis
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    kraken2_db : Path
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BRACKEN             } from '../../../subworkflows/bracken/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BRACKEN(BACTOPIATOOL_INIT.out.reads, params.kraken2_db)

    publish:
    // Per-sample
    sample_outputs = BRACKEN.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(BRACKEN.out.sample_outputs)
    // Run-level
    run_outputs = BRACKEN.out.run_outputs
    run_nf_logs = collectNextflowLogs(BRACKEN.out.run_outputs)
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
