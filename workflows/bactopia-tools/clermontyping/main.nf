#!/usr/bin/env nextflow
/**
 * In silico phylotyping of Escherichia genus.
 *
 * This Bactopia Tool uses [ClermonTyping](https://github.com/happykhan/ClermonTyping)
 * to conduct _in silico_ prediction of phylotype for _Escherichia_ genomes. It uses
 * genome assemblies to assign them to _E. albertii_, _E. fergusonii_, _Escherichia_
 * clades I–V, _E. coli sensu stricto_ as well as to the main _E. coli_ phylogroups.
 *
 * @status stable
 * @keywords e coli, phylotyping, phylogroups, clermontyping, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation clermontyping
 *
 * @subworkflows utils_bactopia-tools, clermontyping
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.blast.xml          BLAST XML file with ClermonTyping analysis results
 * @publish *.html               HTML file with ClermonTyping analysis results
 * @publish *.mash.tsv           TSV file with Mash distances
 * @publish *.phylogroups.txt    TSV file with final phylogroup assignments
 *
 * @section Merged Results
 * @publish clermontyping.tsv    Merged TSV file with ClermonTyping results from all samples
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
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { CLERMONTYPING       } from '../../../subworkflows/clermontyping/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_clermontyping = CLERMONTYPING(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_clermontyping.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_clermontyping.sample_outputs)
    // Run-level
    run_outputs = ch_clermontyping.run_outputs
    run_nf_logs = collectNextflowLogs(ch_clermontyping.run_outputs)
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
