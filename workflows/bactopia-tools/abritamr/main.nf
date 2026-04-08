#!/usr/bin/env nextflow
/**
 * A NATA accredited tool for reporting the presence of antimicrobial resistance genes.
 *
 * This Bactopia Tool uses [abriTAMR](https://github.com/MDU-PHL/abritamr) to identify
 * antimicrobial resistance genes in bacterial genomes. It runs AMRFinderPlus on each
 * sample and collates the results into functional classes, producing detailed reports
 * on resistance genes, partial matches, and virulence factors. It is accredited by NATA
 * for use in reporting presence of reportable AMR genes in Victoria, Australia.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, virulence, amr, nata, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation abritamr, amrfinderplus, csvtk
 *
 * @subworkflows bactopiatool_init, abritamr
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.abritamr.txt           Tab-delimited file combining non-empty summary files from abriTAMR
 * @publish *.amrfinder.out          Raw output from AMRFinderPlus (per sequence)
 * @publish *.summary_matches.txt    Tab-delimited file with AMR gene matches per sequence
 * @publish *.summary_partials.txt   Tab-delimited file with partial AMR gene matches
 * @publish *.summary_virulence.txt  Tab-delimited file with virulence gene classifications
 *
 * @section Merged Results
 * @publish abritamr.tsv             Merged TSV file containing AMR summaries from all samples
 *
 * @section Execution Logs
 * @publish logs/abritamr/*          Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml             Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { ABRITAMR            } from '../../../subworkflows/abritamr/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_abritamr = ABRITAMR(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_abritamr.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_abritamr.sample_outputs)
    // Run-level
    run_outputs = ch_abritamr.run_outputs
    run_nf_logs = collectNextflowLogs(ch_abritamr.run_outputs)
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
