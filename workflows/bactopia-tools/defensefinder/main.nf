#!/usr/bin/env nextflow
/**
 * Systematic identification of anti-phage defense systems.
 *
 * This Bactopia Tool uses [DefenseFinder](https://github.com/mdmparis/defense-finder)
 * to systematically search for and identify all known anti-phage defense systems
 * in bacterial genomes using HMM-based protein domain detection.
 *
 * @status stable
 * @keywords anti-phage, defense systems, hmm, protein domains, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,hmm-search
 *
 * @subworkflows bactopiatool_init, defensefinder
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.prt                          FASTA file containing all proteins found in defense systems
 * @publish *.prt.idx                      Index file for the proteins file
 * @publish *defense_finder_genes.tsv      TSV file with each gene found in defense systems
 * @publish *defense_finder_hmmer.tsv      TSV file with each HMM hit
 * @publish *defense_finder_systems.tsv    TSV file with information about each system found
 * @publish *.macsydata.tar.gz             Raw MACSyFinder output file (requires --defensefinder_preserveraw)
 *
 * @section Merged Results
 * @publish defensefinder-genes.tsv        Merged TSV of all genes found in defense systems
 * @publish defensefinder-hmmer.tsv        Merged TSV of all HMM hits
 * @publish defensefinder-systems.tsv      Merged TSV of all information about systems found
 *
 * @section Execution Logs
 * @publish logs/**                        Tool execution logs
 * @publish logs/nf-*                      Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                   Software version information
 *
 * @citation defensefinder
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { DEFENSEFINDER       } from '../../../subworkflows/defensefinder/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    DEFENSEFINDER(BACTOPIATOOL_INIT.out.proteins)

    publish:
    // Per-sample
    sample_outputs = DEFENSEFINDER.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(DEFENSEFINDER.out.sample_outputs)
    // Run-level
    run_outputs = DEFENSEFINDER.out.run_outputs
    run_nf_logs = collectNextflowLogs(DEFENSEFINDER.out.run_outputs)
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
