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
 * @subworkflows bactopiatool_init, meningotype
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

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MENINGOTYPE       } from '../../../subworkflows/meningotype/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MENINGOTYPE(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = MENINGOTYPE.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MENINGOTYPE.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    publish:
    sample_outputs = MENINGOTYPE.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MENINGOTYPE.out.run_outputs
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
