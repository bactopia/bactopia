//
// pneumocat - Assign capsular type to Streptococcus pneumoniae from sequence reads
//
nextflow.preview.types = true

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow PNEUMOCAT {
    take:
    fastq: Channel<Tuple<Map, Set<Path>>>

    main:
    PNEUMOCAT_MODULE(fastq)

    emit:
    // Individual outputs
    xml: Channel<Tuple<Map, Path>> = PNEUMOCAT_MODULE.out.xml
    txt: Channel<Tuple<Map, Path>> = PNEUMOCAT_MODULE.out.txt

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PNEUMOCAT_MODULE.out.xml,
        PNEUMOCAT_MODULE.out.txt
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.versions])
}
