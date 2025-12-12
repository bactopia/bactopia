#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia.
 *
 * Comprehensive bacterial analysis pipeline for complete genomic characterization.
 * This workflow performs end-to-end analysis including quality control, assembly,
 * annotation, antimicrobial resistance detection, MLST typing, and optional
 * pathogen-specific analysis through Merlin.
 *
 * @status stable
 * @keywords bacteria, assembly, annotation, AMR, MLST, genomics
 *
 * @subworkflows bactopia_init, amrfinderplus, assembler, datasets, gather, sketcher,
 * @subworkflows mlst, qc, bakta, prokka, merlin
 *
 * @input rundir
 * Directory containing raw sequencing reads
 *
 * @input adapters
 * Path to adapter sequences file for removal during QC
 *
 * @input phix
 * Path to PhiX sequences for contamination removal during QC
 *
 * @input use_bakta
 * Use Bakta for genome annotation instead of Prokka
 *
 * @input bakta_db
 * Path to Bakta database for annotation
 *
 * @input download_bakta
 * Download Bakta database if not provided
 *
 * @input bakta_save_as_tarball
 * Save Bakta database as tarball for reuse
 *
 * @input bakta_proteins
 * Path to trusted protein sequences for Bakta annotation
 *
 * @input bakta_prodigal_tf
 * Path to Prodigal training file for Bakta
 *
 * @input bakta_replicons
 * Path to replicon sequences for Bakta
 *
 * @input prokka_proteins
 * Path to protein sequences for Prokka annotation
 *
 * @input prokka_prodigal_tf
 * Path to Prodigal training file for Prokka
 *
 * @input emmtyper_blastdb
 * Path to emmtyper BLAST database for Merlin
 *
 * @input hicap_database_dir
 * Path to HiCap database directory for Merlin
 *
 * @input hicap_model_fp
 * Path to HiCap model file for Merlin
 *
 * @input ask_merlin
 * Enable Merlin pathogen-specific analysis
 *
 * @input spatyper_repeats
 * Path to Spatyper repeats database for Merlin
 *
 * @input spatyper_repeat_order
 * Path to Spatyper repeat order file for Merlin
 *
 * @section Quality Control
 * @publish fastqc/*   FastQC quality control reports
 * @publish multiqc/*   MultiQC aggregated quality reports
 *
 * @section Assembly
 * @publish *.fasta   Assembled genome sequences
 * @publish assembly-stats.txt Assembly quality metrics
 * @publish quast.html   QUAST assembly quality report
 *
 * @section Annotation
 * @publish *.gff    Genome annotation in GFF3 format
 * @publish *.gbk    Genome annotation in GenBank format
 * @publish *.faa    Protein sequences
 * @publish *.fna    Nucleotide sequences
 * @publish *.tsv    Annotation summary tables
 *
 * @section Typing
 * @publish mlst.txt   MLST sequence type results
 *
 * @section Antimicrobial Resistance
 * @publish amrfinderplus.tsv AMR gene detection results
 * @publish amrfinderplus.mutation.tsv AMR mutation results
 *
 * @section Comparative Analysis
 * @publish mash-dist.tsv Mash distance matrix
 * @publish sketch.msh Mash sketch files
 * @publish sourmash.sig Sourmash signatures
 *
 * @section Pathogen-Specific Analysis
 * @note Only created if --ask_merlin is enabled
 * @publish merlin/     Merlin pathogen-specific analysis results
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml Software version information
 */

params {
    rundir   : String

    // Tool-specific parameters
    // QC
    adapters              : Path?
    phix                  : Path?

    // Annotation
    use_bakta             : Boolean
    bakta_db              : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
    bakta_prodigal_tf     : Path? = "${projectDir}/data/empty/EMPTY_PRODIGAL_TF" // TODO: Remove when Path? is fixed
    bakta_replicons       : Path? = "${projectDir}/data/empty/EMPTY_REPLICONS"   // TODO: Remove when Path? is fixed
    prokka_proteins       : Path
    prokka_prodigal_tf    : Path? = "${projectDir}/data/empty/EMPTY_PRODIGAL_TF" // TODO: Remove when Path? is fixed

    // Merlin
    emmtyper_blastdb      : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    hicap_database_dir    : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    hicap_model_fp        : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
    ask_merlin            : Boolean
    spatyper_repeats      : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    spatyper_repeat_order : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
}

// Core
include { BACTOPIA_INIT   } from './subworkflows/utils/bactopia/main'
include { AMRFINDERPLUS   } from './subworkflows/amrfinderplus/main'
include { ASSEMBLER       } from './subworkflows/bactopia/assembler/main'
include { DATASETS        } from './modules/bactopia/datasets/main'
include { GATHER          } from './subworkflows/bactopia/gather/main'
include { SKETCHER        } from './subworkflows/bactopia/sketcher/main'
include { MLST            } from './subworkflows/mlst/main'
include { QC              } from './subworkflows/bactopia/qc/main'

// Annotation wih Bakta or Prokka
include { BAKTA           } from './subworkflows/bakta/main'
include { PROKKA          } from './subworkflows/prokka/main'

// Merlin
include { MERLIN          } from './subworkflows/merlin/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
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
    ch_annotations = channel.empty() as Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>>
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

    // Merlin
    if (params.ask_merlin) {
        MERLIN(
            ASSEMBLER.out.fna_fq,
            DATASETS.out.mash_db,
            // emmtyper
            params.emmtyper_blastdb,
            // hicap
            params.hicap_database_dir,
            params.hicap_model_fp,
            // staphtyper
            params.spatyper_repeats,
            params.spatyper_repeat_order
        )
        ch_results = ch_results.mix(MERLIN.out.results)
        ch_logs = ch_logs.mix(MERLIN.out.logs)
        ch_nf_logs = ch_nf_logs.mix(MERLIN.out.nf_logs)
        ch_versions = ch_versions.mix(MERLIN.out.versions)
    }

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
