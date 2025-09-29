//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//
include { PIRATE   } from '../pirate/main'
include { ROARY    } from '../roary/main'
include { PANAROO  } from '../panaroo/main'
include { SNPDISTS } from '../snpdists/main'

workflow PANGENOME {
    take:
    gff // channel: [ val(meta), [ gff ] ]
    use_pirate
    use_roary

    main:
    ch_aln = Channel.empty()
    ch_csv = Channel.empty()
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()

    // Choose pangenome tool based on params
    if (use_pirate) {
        PIRATE(gff)
        ch_aln = PIRATE.out.aln
        ch_csv = PIRATE.out.csv
        ch_results = PIRATE.out.results
        ch_logs = PIRATE.out.logs
        ch_nf_logs = PIRATE.out.nf_logs
        ch_versions = PIRATE.out.versions
    } else if (use_roary) {
        ROARY(gff)
        ch_aln = ROARY.out.aln
        ch_csv = ROARY.out.csv
        ch_results = ROARY.out.results
        ch_logs = ROARY.out.logs
        ch_nf_logs = ROARY.out.nf_logs
        ch_versions = ROARY.out.versions
    } else {
        PANAROO(gff)
        ch_aln = PANAROO.out.aln
        ch_csv = PANAROO.out.csv
        ch_results = PANAROO.out.results
        ch_logs = PANAROO.out.logs
        ch_nf_logs = PANAROO.out.nf_logs
        ch_versions = PANAROO.out.versions
    }

    // Per-sample SNP distances
    ch_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome.distance", process_name: "snpdists"], aln]}.set{ ch_unmasked_aln }
    SNPDISTS(ch_unmasked_aln)
    ch_results = ch_results.mix(SNPDISTS.out.results)
    ch_logs = ch_logs.mix(SNPDISTS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNPDISTS.out.nf_logs)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    emit:
    // Individual outputs
    aln = ch_aln
    csv = ch_csv

    // Generic aggregate outputs
    results = ch_results
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions
}
