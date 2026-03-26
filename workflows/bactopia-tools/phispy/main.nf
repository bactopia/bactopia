#!/usr/bin/env nextflow
/**
 * Prediction of prophages in bacterial and archaeal genomes.
 *
 * This Bactopia Tool uses [PhiSpy](https://github.com/linsalrob/PhiSpy) to identify prophages
 * from bacterial and archaeal genomes using machine learning approaches.
 *
 * @status stable
 * @keywords prophage, phage, annotation, bacterial, archaeal, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,phage-detection
 * @citation phispy
 *
 * @subworkflows bactopiatool_init, phispy
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Tab-delimited file containing prophage predictions
 *
 * @section Merged Results
 * @publish phispy.tsv       Merged TSV file containing prophage predictions from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PHISPY            } from '../../../subworkflows/phispy/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PHISPY(BACTOPIATOOL_INIT.out.gffs)
    ch_sample_nf_logs = PHISPY.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = PHISPY.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = PHISPY.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = PHISPY.out.run_outputs
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
