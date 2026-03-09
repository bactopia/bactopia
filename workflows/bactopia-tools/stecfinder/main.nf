#!/usr/bin/env nextflow
/**
 * Serotype identification of Shiga toxin-producing E. coli.
 *
 * This Bactopia Tool uses [STECFinder](https://github.com/LanLab/STECFinder) to identify
 * the serotype of Shiga toxin-producing _E. coli_ (STEC) from sequencing data.
 * STECFinder determines the serotype as well as the O-antigen and H-antigens.
 *
 * @status stable
 * @keywords stec, serotype, e coli, shiga toxin, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,serotyping
 * @citation stecfinder
 *
 * @subworkflows bactopiatool_init, stecfinder
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Serotype identification results
 *
 * @section Merged Results
 * @publish stecfinder.tsv   Merged TSV file containing STECFinder results from all samples
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
include { STECFINDER        } from '../../../subworkflows/stecfinder/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    STECFINDER(BACTOPIATOOL_INIT.out.assembly_reads)
    ch_sample_nf_logs = STECFINDER.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = STECFINDER.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = STECFINDER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = STECFINDER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
