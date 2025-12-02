#!/usr/bin/env nextflow
nextflow.preview.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW PARAMETERS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params {
    rundir   : String

    // Tool-specific parameters
    adapters              : Path?
    phix                  : Path?
    use_bakta             : Boolean
    bakta_db              : Path?
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path?
    bakta_prodigal_tf     : Path?
    bakta_replicons       : Path?
    prokka_proteins       : Path?
    prokka_prodigal_tf    : Path?
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

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

// Annotation wih Bakta or Prokka
include { BAKTA           } from '../../subworkflows/bakta/main'
include { PROKKA          } from '../../subworkflows/prokka/main'

// Merlin
include { STAPHTYPER      } from '../../subworkflows/staphtyper/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Initialize and execute the workflow
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>
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
        params.adapters,
        params.phix
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
    ch_annotations = channel.empty()
    if (params.use_bakta) {
        BAKTA(
            ASSEMBLER.out.fna,
            params.bakta_db,
            params.download_bakta,
            params.bakta_save_as_tarball,
            params.bakta_proteins,
            params.bakta_prodigal_tf,
            params.bakta_replicons
        )
        ch_results = ch_results.mix(BAKTA.out.results)
        ch_logs = ch_logs.mix(BAKTA.out.logs)
        ch_nf_logs = ch_nf_logs.mix(BAKTA.out.nf_logs)
        ch_versions = ch_versions.mix(BAKTA.out.versions)
        ch_annotations = ch_annotations.mix(BAKTA.out.annotations)
    } else {
        PROKKA(
            ASSEMBLER.out.fna,
            params.prokka_proteins,
            params.prokka_prodigal_tf
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

    // Staphtyper
    STAPHTYPER(
        ASSEMBLER.out.fna,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )
    ch_results = ch_results.mix(STAPHTYPER.out.results)
    ch_logs = ch_logs.mix(STAPHTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(STAPHTYPER.out.nf_logs)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
