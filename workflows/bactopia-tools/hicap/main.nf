#!/usr/bin/env nextflow
/**
 * Identify cap locus serotype and structure in Haemophilus influenzae assemblies.
 *
 * This Bactopia Tool uses [hicap](https://github.com/scwatts/hicap) with assemblies for
 * _in silico_ typing of the _Haemophilus influenzae_ capsular locus.
 *
 * @status stable
 * @keywords haemophilus influenzae, serotyping, capsular locus, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation hicap
 *
 * @subworkflows bactopiatool_init, hicap
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input hicap_database_dir
 * Path to hicap database directory
 *
 * @input hicap_model_fp
 * Path to hicap model file
 *
 * @section Per-Sample Results
 * @publish *.summary       Summary of serotype prediction
 * @publish *.gff           Annotated capsular locus in GFF format
 *
 * @section Merged Results
 * @publish hicap.tsv        Merged TSV file containing hicap results from all samples
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
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    hicap_database_dir : Path
    hicap_model_fp     : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { HICAP             } from '../../../subworkflows/hicap/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    HICAP(
        BACTOPIATOOL_INIT.out.assembly,
        params.hicap_database_dir,
        params.hicap_model_fp
    )

    ch_sample_nf_logs = HICAP.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = HICAP.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = HICAP.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = HICAP.out.run_outputs
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
