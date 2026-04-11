#!/usr/bin/env nextflow
/**
 * Comprehensive typing of Neisseria meningitidis.
 *
 * This Bactopia Tool uses [meningotype](https://github.com/MDU-PHL/meningotype)
 * for _in silico_ typing of _Neisseria meningitidis_ genomes. It uses assembly contigs
 * to determine the serotype, MLST, finetyping (_porA_, _fetA_, _porB_), and
 * Bexsero antigen sequence typing (BAST) (_fHbp_, _NHBA_, _NadA_, _PorA_).
 *
 * @status stable
 * @keywords neisseria meningitidis, serotyping, mlst, finetyping, fasta, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing,mlst
 * @citation meningotype
 *
 * @subworkflows utils_bactopia-tools, meningotype
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt                Comprehensive typing report
 * @publish *-allele.tsv         Allele typing results
 * @publish *-mlst.tsv           MLST typing results
 *
 * @section Merged Results
 * @publish meningotype.tsv       Merged TSV file containing meningotype results from all samples
 *
 * @section Execution Logs
 * @publish logs/**               Tool execution logs
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { MENINGOTYPE         } from '../../../subworkflows/meningotype/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_meningotype = MENINGOTYPE(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_meningotype.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_meningotype.sample_outputs)
    // Run-level
    run_outputs = ch_meningotype.run_outputs
    run_nf_logs = collectNextflowLogs(ch_meningotype.run_outputs)
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
