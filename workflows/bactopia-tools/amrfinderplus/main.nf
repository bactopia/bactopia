#!/usr/bin/env nextflow
/**
 * Bactopia Tool: Amrfinderplus.
 *
 * Identify antimicrobial resistance genes and point mutations in bacterial genomes.
 * This Bactopia Tool uses [AMRFinder+](https://github.com/ncbi/amr) to screen assemblies and proteins
 * for antimicrobial resistance genes, virulence genes, and resistance-associated point mutations.
 * It identifies acquired AMR genes and some point mutations in protein or assembled nucleotide sequences.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance, virulence, genes, proteins, mutations
 * @tags complexity:moderate input-type:directory output-type:multiple features:database-dependent,amr-detection
 * @citation amrfinderplus
 *
 * @subworkflows bactopiatool_init, amrfinderplus
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @input amrfinderplus_db
 * Path to AMRFinder+ database for AMR gene detection
 *
 * @section Per-Sample Results
 * @publish *.tsv                    AMR gene detection results in TSV format
 * @publish *-mutations.tsv          Point mutations associated with antimicrobial resistance
 *
 * @section Merged Results
 * @publish amrfinderplus.tsv        Combined AMR detection results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                  Tool execution logs
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml              Software version information
   */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    amrfinderplus_db : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { AMRFINDERPLUS       } from '../../../subworkflows/amrfinderplus/main'
include { DATASETS            } from '../../../subworkflows/bactopia/datasets/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()

    if (params.amrfinderplus_db) {
        ch_amrfinderplus = AMRFINDERPLUS(
            ch_bactopiatool.assembly_proteins_gff,
            params.amrfinderplus_db
        )
    } else {
        ch_datasets = DATASETS()
        ch_amrfinderplus = AMRFINDERPLUS(
            ch_bactopiatool.assembly_proteins_gff,
            ch_datasets.amrfinderplus_db
        )
    }

    publish:
    // Per-sample
    sample_outputs = ch_amrfinderplus.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_amrfinderplus.sample_outputs)
    // Run-level
    run_outputs = ch_amrfinderplus.run_outputs
    run_nf_logs = collectNextflowLogs(ch_amrfinderplus.run_outputs)
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
