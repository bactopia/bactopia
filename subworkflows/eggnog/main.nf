//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
nextflow.preview.types = true

include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../modules/eggnog/mapper/main'
include { flattenPaths    } from 'plugin/nf-bactopia'
include { gather          } from 'plugin/nf-bactopia'

workflow EGGNOG {
    take:
    faa: Channel<Tuple<Map, Set<Path>>>
    database: Path
    download_eggnog: Boolean

    main:
    if (download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(faa, database)
    }

    emit:
    // Individual outputs
    hits: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.hits
    seed_orthologs: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.seed_orthologs
    annotations: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.annotations
    xlsx: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.xlsx
    orthologs: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.orthologs
    genepred: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.genepred
    gff: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.gff
    no_anno: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.no_anno
    pfam: Channel<Tuple<Map, Path>> = EGGNOG_MAPPER.out.pfam

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        EGGNOG_MAPPER.out.hits,
        EGGNOG_MAPPER.out.seed_orthologs,
        EGGNOG_MAPPER.out.annotations,
        EGGNOG_MAPPER.out.xlsx,
        EGGNOG_MAPPER.out.orthologs,
        EGGNOG_MAPPER.out.genepred,
        EGGNOG_MAPPER.out.gff,
        EGGNOG_MAPPER.out.no_anno,
        EGGNOG_MAPPER.out.pfam
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([EGGNOG_MAPPER.out.versions])
}
