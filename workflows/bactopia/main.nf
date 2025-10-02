#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Core
include { BACTOPIA_INIT   } from '../../subworkflows/utils/bactopia'
include { AMRFINDERPLUS   } from '../../subworkflows/amrfinderplus/main'
include { ASSEMBLER       } from '../../subworkflows/bactopia/assembler/main'
include { DATASETS        } from '../../modules/bactopia/datasets/main'
include { GATHER          } from '../../subworkflows/bactopia/gather/main'
include { SKETCHER        } from '../../subworkflows/bactopia/sketcher/main'
include { MLST            } from '../../subworkflows/mlst/main'
include { QC              } from '../../subworkflows/bactopia/qc/main'
include { paramsHelp      } from 'plugin/nf-bactopia'
include { workflowSummary } from 'plugin/nf-bactopia'

// Annotation wih Bakta or Prokka
include { BAKTA           } from '../../subworkflows/bakta/main'
include { PROKKA          } from '../../subworkflows/prokka/main'

// Merlin
include { MERLIN          } from '../../subworkflows/merlin/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Check if help is requested
    if (params.help || params.help_all) {
        log.info paramsHelp()
        exit 0
    }


    // Initialize and execute the workflow
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()
    BACTOPIA_INIT()

    // Core steps
    DATASETS()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)
    ch_results = ch_results.mix(GATHER.out.results)
    ch_logs = ch_logs.mix(GATHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
    ch_versions = ch_versions.mix(GATHER.out.versions)

    // QC samples
    QC(
        GATHER.out.raw_fastq,
        params.adapters ? file(params.adapters) : [],
        params.phix ? file(params.phix) : []
    )
    ch_results = ch_results.mix(QC.out.results)
    ch_logs = ch_logs.mix(QC.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QC.out.nf_logs)
    ch_versions = ch_versions.mix(QC.out.versions)

    // Assemble genomes
    ASSEMBLER(QC.out.fastq)
    ch_results = ch_results.mix(ASSEMBLER.out.results)
    ch_logs = ch_logs.mix(ASSEMBLER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ASSEMBLER.out.nf_logs)
    ch_versions = ch_versions.mix(ASSEMBLER.out.versions)

    // Sketch and query
    SKETCHER(ASSEMBLER.out.fna, DATASETS.out.mash_db, DATASETS.out.sourmash_db)
    ch_results = ch_results.mix(SKETCHER.out.results)
    ch_logs = ch_logs.mix(SKETCHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SKETCHER.out.nf_logs)
    ch_versions = ch_versions.mix(SKETCHER.out.versions)

    // Annotate samples
    ch_annotations = Channel.empty()
    if (params.use_bakta) {
        BAKTA(
            ASSEMBLER.out.fna,
            params.bakta_db ? file(params.bakta_db) : [],
            params.download_bakta,
            params.bakta_save_as_tarball,
            params.proteins ? file(params.proteins) : [],
            params.prodigal_tf ? file(params.prodigal_tf) : [],
            params.replicons ? file(params.replicons) : []
        )
        ch_results = ch_results.mix(BAKTA.out.results)
        ch_logs = ch_logs.mix(BAKTA.out.logs)
        ch_nf_logs = ch_nf_logs.mix(BAKTA.out.nf_logs)
        ch_versions = ch_versions.mix(BAKTA.out.versions)
        ch_annotations = ch_annotations.mix(BAKTA.out.annotations)
    } else {
        PROKKA(
            ASSEMBLER.out.fna,
            params.proteins ? file(params.proteins) : [],
            params.prodigal_tf ? file(params.prodigal_tf) : []
        )
        ch_results = ch_results.mix(PROKKA.out.results)
        ch_logs = ch_logs.mix(PROKKA.out.logs)
        ch_nf_logs = ch_nf_logs.mix(PROKKA.out.nf_logs)
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        ch_annotations = ch_annotations.mix(PROKKA.out.annotations)
    }

    // AMR
    AMRFINDERPLUS(ch_annotations, DATASETS.out.amrfinderplus_db)
    ch_results = ch_results.mix(AMRFINDERPLUS.out.results)
    ch_logs = ch_logs.mix(AMRFINDERPLUS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(AMRFINDERPLUS.out.nf_logs)
    ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions)

    // MLST
    MLST(ASSEMBLER.out.fna, DATASETS.out.mlst_db)
    ch_results = ch_results.mix(MLST.out.results)
    ch_logs = ch_logs.mix(MLST.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MLST.out.nf_logs)
    ch_versions = ch_versions.mix(MLST.out.versions)

    if (params.ask_merlin) {
        MERLIN(
            ASSEMBLER.out.fna.join(QC.out.fastq_only, by:[0]),
            DATASETS.out.mash_db,
            // emmtyper
            params.emmtyper_blastdb ? file(params.emmtyper_blastdb, checkIfExists: true) : [],
            // hicap
            params.database_dir ? file(params.database_dir, checkIfExists: true) : [],
            params.model_fp ? file(params.model_fp, checkIfExists: true) : [],
            // staphtyper
            params.repeats ? file(params.repeats, checkIfExists: true) : [],
            params.repeat_order ? file(params.repeat_order, checkIfExists: true) : []
        )
        ch_results = ch_results.mix(MERLIN.out.results)
        ch_logs = ch_logs.mix(MERLIN.out.logs)
        ch_nf_logs = ch_nf_logs.mix(MERLIN.out.nf_logs)
        ch_versions = ch_versions.mix(MERLIN.out.versions)
    }

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = ch_results
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
