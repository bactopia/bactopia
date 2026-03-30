#!/usr/bin/env nextflow
/**
 * Serotyping of Streptococcus pneumoniae from Illumina paired-end reads.
 *
 * This Bactopia Tool uses [Seroba](https://github.com/sanger-pathogens/seroba) to predict the
 * serotype of *Streptococcus pneumoniae* samples from raw sequencing reads. Seroba employs
 * a k-mer based approach to identify and type pneumococcal capsules, determining both the
 * serotype and serogroup based on the presence of specific capsular polysaccharide synthesis
 * (cps) locus sequences. The tool is specifically designed for Illumina paired-end reads and
 * provides accurate serotype predictions essential for epidemiological surveillance and vaccine
 * development studies.
 *
 * @status stable
 * @keywords serotyping, Streptococcus pneumoniae, capsule, cps locus, vaccine, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, seroba
 *
 * @subworkflows bactopiatool_init, seroba
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                              Tab-delimited file containing predicted serotype for each sample
 * @publish *detailed_serogroup_info.txt      Detailed information about serotype prediction and coverage
 *
 * @section Merged Results
 * @publish seroba.tsv                         Merged TSV file containing serotype predictions from all samples
 *
 * @section Execution Logs
 * @publish logs/**                            Tool execution logs including Seroba output
 * @publish logs/nf-*                          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                       Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SEROBA              } from '../../../subworkflows/seroba/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SEROBA(BACTOPIATOOL_INIT.out.reads)

    publish:
    // Per-sample
    sample_outputs = SEROBA.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(SEROBA.out.sample_outputs)
    // Run-level
    run_outputs = SEROBA.out.run_outputs
    run_nf_logs = collectNextflowLogs(SEROBA.out.run_outputs)
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
