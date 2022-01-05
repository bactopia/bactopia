//
// bakta - Rapid annotation of bacterial genomes and plasmids
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'bakta')
options.is_module = params.wf == 'bakta' ? true : false
options.args =  = [
    params.skip_trna ? "--skip-trna" : "",
    params.skip_tmrna ? "--skip-tmrna" : "",
    params.skip_rrna ? "--skip-rrna" : "",
    params.skip_ncrna ? "--skip-ncrna" : "",
    params.skip_ncrna_region ? "--skip-ncrna-region" : "",
    params.skip_crispr ? "--skip-crispr" : "",
    params.skip_cds ? "--skip-cds" : "",
    params.skip_sorf ? "--skip-sorf" : "",
    params.skip_gap ? "--skip-gap" : "",
    params.skip_ori ? "--skip-ori" : "",
    params.compliant ? "--compliant" : "",
    params.keep_contig_headers ? "--keep-contig-headers" : "",
    "--min-contig-length ${params.min_contig_length}",
    "${params.bakta_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE_DIR = params.bakta_db ? file(params.bakta_db) : []
PROTEINS = params.proteins ? file(params.proteins) : []
PRODIGAL_TF = params.prodigal_tf ? file(params.prodigal_tf) : []
REPLICONS = params.replicons ? file(params.replicons) : []

include { BAKTA as BAKTA_MODULE } from '../../../modules/nf-core/modules/bakta/main' addParams( options: options )

workflow BAKTA {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    BAKTA_MODULE(fasta, DATABASE_DIR, PROTEINS, PRODIGAL_TF, REPLICONS)
    ch_versions = ch_versions.mix(BAKTA_MODULE.out.versions.first())

    emit:
    embl = BAKTA_MODULE.out.embl
    faa = BAKTA_MODULE.out.faa
    ffn = BAKTA_MODULE.out.ffn
    fna = BAKTA_MODULE.out.fna
    gbff = BAKTA_MODULE.out.gbff
    gff = BAKTA_MODULE.out.gff
    hypotheticals_tsv = BAKTA_MODULE.out.hypotheticals_tsv
    hypotheticals_faa = BAKTA_MODULE.out.hypotheticals_faa
    tsv = BAKTA_MODULE.out.tsv
    versions = ch_versions // channel: [ versions.yml ]
}
