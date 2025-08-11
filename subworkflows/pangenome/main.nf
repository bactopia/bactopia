//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//
include { PIRATE } from '../pirate/main'
include { ROARY } from '../roary/main'
include { PANAROO } from '../panaroo/main'
include { CLONALFRAMEML } from '../clonalframeml/main'
include { IQTREE as FINAL_TREE } from '../iqtree/main'
include { SNPDISTS as SNPDISTS_UNMASKED } from '../../modules/snpdists/main'
include { SNPDISTS as SNPDISTS_MASKED } from '../../modules/snpdists/main'
include { SCOARY } from '../scoary/main'

workflow PANGENOME {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ch_needs_prokka = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Choose pangenome tool based on params
    if (params.use_pirate) {
        PIRATE(gff)
        ch_pg_versions = PIRATE.out.versions
        ch_pg_logs = PIRATE.out.logs
        ch_pg_nf_logs = PIRATE.out.nf_logs
        ch_pg_aln = PIRATE.out.aln
        ch_pg_csv = PIRATE.out.csv
    } else if (params.use_roary) {
        ROARY(gff)
        ch_pg_versions = ROARY.out.versions
        ch_pg_logs = ROARY.out.logs
        ch_pg_nf_logs = ROARY.out.nf_logs
        ch_pg_aln = ROARY.out.aln
        ch_pg_csv = ROARY.out.csv
    } else {
        PANAROO(gff)
        ch_pg_versions = PANAROO.out.versions
        ch_pg_logs = PANAROO.out.logs
        ch_pg_nf_logs = PANAROO.out.nf_logs
        ch_pg_aln = PANAROO.out.aln
        ch_pg_csv = PANAROO.out.csv
    }
    ch_versions = ch_versions.mix(ch_pg_versions)
    ch_logs = ch_logs.mix(ch_pg_logs)
    ch_nf_logs = ch_nf_logs.mix(ch_pg_nf_logs)

    // Per-sample SNP distances
    SNPDISTS_UNMASKED(ch_pg_aln)
    ch_versions = ch_versions.mix(SNPDISTS_UNMASKED.out.versions)
    ch_logs = ch_logs.mix(SNPDISTS_UNMASKED.out.logs)
    ch_nf_logs = ch_nf_logs.mix(
        SNPDISTS_UNMASKED.out.nf_begin,
        SNPDISTS_UNMASKED.out.nf_err,
        SNPDISTS_UNMASKED.out.nf_log,
        SNPDISTS_UNMASKED.out.nf_out,
        SNPDISTS_UNMASKED.out.nf_run,
        SNPDISTS_UNMASKED.out.nf_sh,
        SNPDISTS_UNMASKED.out.nf_trace
    )
    
    // Identify Recombination
    if (!params.skip_recombination) {
        // Run ClonalFrameML
        CLONALFRAMEML(ch_pg_aln)
        ch_versions = ch_versions.mix(CLONALFRAMEML.out.versions)
        ch_logs = ch_logs.mix(CLONALFRAMEML.out.logs)
        ch_nf_logs = ch_nf_logs.mix(CLONALFRAMEML.out.nf_logs)
    }

    // Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        if (params.skip_recombination) {
            FINAL_TREE(ch_pg_aln)
        } else {
            FINAL_TREE(CLONALFRAMEML.out.masked_aln)
            SNPDISTS_MASKED(CLONALFRAMEML.out.masked_aln)
            ch_logs = ch_logs.mix(SNPDISTS_MASKED.out.logs)
            ch_nf_logs = ch_nf_logs.mix(
                SNPDISTS_MASKED.out.nf_begin,
                SNPDISTS_MASKED.out.nf_err,
                SNPDISTS_MASKED.out.nf_log,
                SNPDISTS_MASKED.out.nf_out,
                SNPDISTS_MASKED.out.nf_run,
                SNPDISTS_MASKED.out.nf_sh,
                SNPDISTS_MASKED.out.nf_trace
            )
        }
        ch_versions = ch_versions.mix(FINAL_TREE.out.versions)
        ch_logs = ch_logs.mix(FINAL_TREE.out.logs)
        ch_nf_logs = ch_nf_logs.mix(FINAL_TREE.out.nf_logs)
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY(ch_pg_csv)
        ch_versions = ch_versions.mix(SCOARY.out.versions)
        ch_logs = ch_logs.mix(SCOARY.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SCOARY.out.nf_logs)
    }

    emit:
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
