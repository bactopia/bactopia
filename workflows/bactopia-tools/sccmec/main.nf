#!/usr/bin/env nextflow
/**
 * Typing of SCCmec cassettes in Staphylococcus aureus assemblies.
 *
 * This Bactopia Tool uses [sccmec](https://github.com/rpetit3/sccmec) to identify and type
 * Staphylococcal Cassette Chromosome mec (SCCmec) elements in *Staphylococcus aureus* assemblies.
 * SCCmec cassettes are mobile genetic elements that carry the mecA gene and other methicillin
 * resistance determinants. The tool performs BLAST searches against target-specific sequences
 * and full cassette references to determine SCCmec types and subtypes, providing detailed
 * reports of BLAST hits and type predictions for epidemiological surveillance of MRSA.
 *
 * @status stable
 * @keywords resistance, Staphylococcus aureus, MRSA, SCCmec, typing, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, sccmec
 *
 * @subworkflows bactopiatool_init, sccmec
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                           Summary report of predicted SCCmec type for each sample
 * @publish *.targets.blastn.tsv            Tab-delimited file of all target-specific BLAST hits
 * @publish *.targets.details.tsv           Detailed breakdown of type predictions based on target hits
 * @publish *.regions.blastn.tsv            Tab-delimited file of all full cassette BLAST hits
 * @publish *.regions.details.tsv           Detailed breakdown of type predictions based on full cassettes
 *
 * @section Merged Results
 * @publish sccmec.tsv                      Merged TSV file containing SCCmec typing results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                         Tool execution logs including BLAST output
 * @publish logs/nf-*                       Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                    Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SCCMEC              } from '../../../subworkflows/sccmec/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SCCMEC(BACTOPIATOOL_INIT.out.assembly)

    publish:
    // Per-sample
    sample_outputs = SCCMEC.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(SCCMEC.out.sample_outputs)
    // Run-level
    run_outputs = SCCMEC.out.run_outputs
    run_nf_logs = collectNextflowLogs(SCCMEC.out.run_outputs)
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
