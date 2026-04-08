#!/usr/bin/env nextflow
/**
 * Automatic Multi-Locus Sequence Type (MLST) calling from assembled contigs.
 *
 * This Bactopia Tool uses [mlst](https://github.com/tseemann/mlst) to scan genome assemblies
 * and determine the sequence type based on [PubMLST](https://pubmlst.org/) schemes. The workflow
 * automatically detects the appropriate MLST scheme for each organism and provides
 * standardized sequence type assignments.
 *
 * @status stable
 * @keywords typing, MLST, sequence type, alleles, PubMLST
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, mlst, pubmlst
 *
 * @subworkflows bactopiatool_init, mlst
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                      Tab-delimited file with MLST results including scheme, ST, and allele profiles
 *
 * @section Merged Results
 * @publish mlst.tsv                   Merged TSV file containing MLST results from all samples
 *
 * @section Execution Logs
 * @publish logs/mlst/*                Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    mlst_db : Value<Path?>
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { MLST                } from '../../../subworkflows/mlst/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_mlst = MLST(ch_bactopiatool.assembly, params.mlst_db)

    publish:
    // Per-sample
    sample_outputs = ch_mlst.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_mlst.sample_outputs)
    // Run-level
    run_outputs = ch_mlst.run_outputs
    run_nf_logs = collectNextflowLogs(ch_mlst.run_outputs)
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
