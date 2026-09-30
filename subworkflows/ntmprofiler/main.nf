/**
 * Identify NTM species and predict antimicrobial resistance from sequencing reads.
 *
 * This subworkflow runs [NTM-Profiler](https://github.com/jodyphelan/NTM-Profiler)
 * for each sample and collates the machine-readable results across the run.
 *
 * @status stable
 * @keywords mycobacterium, ntm, species identification, antimicrobial resistance, surveillance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation,database-dependent
 * @citation ntmprofiler
 *
 * @modules ntmprofiler_update, ntmprofiler_profile, ntmprofiler_collate
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db?
 * Optional directory or compressed tarball containing a custom NTM-Profiler database.
 * When omitted, the default database is prepared once and reused by all tasks.
 *
 * @output sample_outputs
 * - `csv`: Per-sample results in CSV format
 * - `json`: Per-sample compressed machine-readable results
 * - `txt`: Per-sample human-readable profiling report
 *
 * @output run_outputs
 * - `csv`: Collated species, lineage, and resistance results
 * - `variants_csv`: Collated variant calls
 * - `snp_dist_json?`: SNP-distance graph when a distance database is present
 */
nextflow.enable.types = true

include { NTMPROFILER_PROFILE } from '../../modules/ntmprofiler/profile/main'
include { NTMPROFILER_COLLATE } from '../../modules/ntmprofiler/collate/main'
include { NTMPROFILER_UPDATE  } from '../../modules/ntmprofiler/update/main'
include { gather              } from 'plugin/nf-bactopia'

workflow NTMPROFILER {
    take:
    reads: Channel<Record>
    db: Path?

    main:
    if (db == null) {
        ch_ntmprofiler_db = NTMPROFILER_UPDATE().map { r -> r.db }
    } else {
        ch_ntmprofiler_db = db
    }

    ch_ntmprofiler_profile = NTMPROFILER_PROFILE(reads, ch_ntmprofiler_db)
    ch_ntmprofiler_collate = NTMPROFILER_COLLATE(
        gather(ch_ntmprofiler_profile, 'json', [name: 'ntmprofiler']),
        ch_ntmprofiler_db
    )

    emit:
    // Published outputs
    sample_outputs = ch_ntmprofiler_profile
    run_outputs = ch_ntmprofiler_collate
}
