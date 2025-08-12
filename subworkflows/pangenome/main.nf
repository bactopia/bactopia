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

    ch_logs = Channel.empty()
    // Choose pangenome tool based on params
    if (params.use_pirate) {
        PIRATE(gff)

    } else if (params.use_roary) {
        ROARY(gff)

    } else {
        PANAROO(gff)

    }
    ch_versions = ch_versions.mix(ch_pg_versions)
    ch_logs = ch_logs.mix(ch_pg_logs)
    // Per-sample SNP distances
    SNPDISTS_UNMASKED(ch_pg_aln)
    ch_versions = ch_versions.mix(SNPDISTS_UNMASKED.out.versions)
    ch_logs = ch_logs.mix(SNPDISTS_UNMASKED.out.logs)
    // Identify Recombination
    if (!params.skip_recombination) {
        // Run ClonalFrameML
        CLONALFRAMEML(ch_pg_aln)
        ch_versions = ch_versions.mix(CLONALFRAMEML.out.versions)
        ch_logs = ch_logs.mix(CLONALFRAMEML.out.logs)
    }

    // Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        if (params.skip_recombination) {
            FINAL_TREE(ch_pg_aln)
        } else {
            FINAL_TREE(CLONALFRAMEML.out.masked_aln)
            SNPDISTS_MASKED(CLONALFRAMEML.out.masked_aln)
            ch_logs = ch_logs.mix(SNPDISTS_MASKED.out.logs)
        }
        ch_versions = ch_versions.mix(FINAL_TREE.out.versions)
        ch_logs = ch_logs.mix(FINAL_TREE.out.logs)
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY(ch_pg_csv)
        ch_versions = ch_versions.mix(SCOARY.out.versions)
        ch_logs = ch_logs.mix(SCOARY.out.logs)
    }

    emit:
    ch_needs_prokka = Channel.empty()
    ch_pg_aln = PANAROO.out.aln
    ch_pg_csv = PANAROO.out.csv
    ch_pg_logs = PANAROO.out.logs
    ch_pg_nf_logs = PANAROO.out.nf_logs
    ch_pg_versions = PANAROO.out.versions
    logs = ch_logs
    nf_logs = CLONALFRAMEML.out.nf_begin.mix(
        CLONALFRAMEML.out.nf_err,
        CLONALFRAMEML.out.nf_log,
        CLONALFRAMEML.out.nf_out,
        CLONALFRAMEML.out.nf_run,
        CLONALFRAMEML.out.nf_sh,
        CLONALFRAMEML.out.nf_trace,
        FINAL_TREE.out.nf_begin,
        FINAL_TREE.out.nf_err,
        FINAL_TREE.out.nf_log,
        FINAL_TREE.out.nf_out,
        FINAL_TREE.out.nf_run,
        FINAL_TREE.out.nf_sh,
        FINAL_TREE.out.nf_trace,
        SCOARY.out.nf_begin,
        SCOARY.out.nf_err,
        SCOARY.out.nf_log,
        SCOARY.out.nf_out,
        SCOARY.out.nf_run,
        SCOARY.out.nf_sh,
        SCOARY.out.nf_trace,
        SNPDISTS_MASKED.out.nf_begin,
        SNPDISTS_MASKED.out.nf_err,
        SNPDISTS_MASKED.out.nf_log,
        SNPDISTS_MASKED.out.nf_out,
        SNPDISTS_MASKED.out.nf_run,
        SNPDISTS_MASKED.out.nf_sh,
        SNPDISTS_MASKED.out.nf_trace,
        SNPDISTS_UNMASKED.out.nf_begin,
        SNPDISTS_UNMASKED.out.nf_err,
        SNPDISTS_UNMASKED.out.nf_log,
        SNPDISTS_UNMASKED.out.nf_out,
        SNPDISTS_UNMASKED.out.nf_run,
        SNPDISTS_UNMASKED.out.nf_sh,
        SNPDISTS_UNMASKED.out.nf_trace
    )
    versions = ch_versions
}
