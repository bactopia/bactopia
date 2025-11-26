//
// merlin - MinmER assisted species-specific bactopia tool seLectIoN
//

include { MERLINDIST    } from '../mashdist/main';
include { CLERMONTYPING } from '../clermontyping/main';
include { ECTYPER       } from '../ectyper/main';
include { EMMTYPER      } from '../emmtyper/main';
include { HICAP         } from '../hicap/main';
include { HPSUISSERO    } from '../hpsuissero/main';
include { GENOTYPHI     } from '../genotyphi/main';
include { KLEBORATE     } from '../kleborate/main';
include { LEGSTA        } from '../legsta/main';
include { LISSERO       } from '../lissero/main';
include { MENINGOTYPE   } from '../meningotype/main';
include { NGMASTER      } from '../ngmaster/main';
include { PASTY         } from '../pasty/main';
include { PBPTYPER      } from '../pbptyper/main';
include { SEQSERO2      } from '../seqsero2/main';
include { SEROBA        } from '../seroba/main';
include { SHIGAPASS     } from '../shigapass/main';
include { SHIGATYPER    } from '../shigatyper/main';
include { SHIGEIFINDER  } from '../shigeifinder/main';
include { SISTR         } from '../sistr/main';
include { SSUISSERO     } from '../ssuissero/main';
include { STAPHTYPER    } from '../staphtyper/main';
include { STECFINDER    } from '../stecfinder/main';
include { TBPROFILER    } from '../tbprofiler/main';
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow MERLIN {
    take:
    assembly // channel: [ val(meta), [ assembly ] ]
    mash_db // channel: [ mash_db ]
    emmtyper_blastdb
    hicap_database_dir
    hicap_model_fp
    staphtyper_repeats
    staphtyper_repeat_order

    main:
    // ID potential species
    MERLINDIST(assembly, mash_db)

    // Escherichia/Shigella
    ch_escherichia = MERLINDIST.out.escherichia.map{_meta, _assembly, _found -> [_meta, _assembly]}
    ch_escherichia_fq = MERLINDIST.out.escherichia_fq.map{_meta, _reads, _found -> [_meta, _reads]}
    ch_escherichia_fna_fq = MERLINDIST.out.escherichia_fna_fq.map{_meta, _assembly, _reads, _found -> [_meta, _assembly, _reads]}
    CLERMONTYPING(ch_escherichia)
    ECTYPER(ch_escherichia)
    SHIGAPASS(ch_escherichia)
    SHIGATYPER(ch_escherichia_fq)
    SHIGEIFINDER(ch_escherichia)
    STECFINDER(ch_escherichia_fna_fq)

    // Haemophilus
    ch_haemophilus = MERLINDIST.out.haemophilus.map{_meta, _assembly, _found -> [_meta, _assembly]}
    HICAP(ch_haemophilus, hicap_database_dir, hicap_model_fp)
    HPSUISSERO(ch_haemophilus)

    // Klebsiella
    ch_klebsiella = MERLINDIST.out.klebsiella.map{_meta, _assembly, _found -> [_meta, _assembly]}
    KLEBORATE(ch_klebsiella)

    // Legionella
    ch_legionella = MERLINDIST.out.legionella.map{_meta, _assembly, _found -> [_meta, _assembly]}
    LEGSTA(ch_legionella)

    // Listeria
    ch_listeria = MERLINDIST.out.listeria.map{_meta, _assembly, _found -> [_meta, _assembly]}
    LISSERO(ch_listeria)

    // Mycobacterium
    ch_mycobacterium_fq = MERLINDIST.out.mycobacterium_fq.map{_meta, _reads, _found -> [_meta, _reads]}
    TBPROFILER(ch_mycobacterium_fq)

    // Neisseria
    ch_neisseria = MERLINDIST.out.neisseria.map{_meta, _assembly, _found -> [_meta, _assembly]}
    MENINGOTYPE(ch_neisseria)
    NGMASTER(ch_neisseria)

    // Pseudomonas
    ch_pseudomonas = MERLINDIST.out.pseudomonas.map{_meta, _assembly, _found -> [_meta, _assembly]}
    PASTY(ch_pseudomonas)

    // Salmonella
    ch_salmonella = MERLINDIST.out.salmonella.map{_meta, _assembly, _found -> [_meta, _assembly]}
    ch_salmonella_fq = MERLINDIST.out.salmonella_fq.map{_meta, _reads, _found -> [_meta, _reads]}
    GENOTYPHI(ch_salmonella_fq)
    SEQSERO2(ch_salmonella)
    SISTR(ch_salmonella)

    // Staphylococcus
    ch_staphylococcus = MERLINDIST.out.staphylococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}
    STAPHTYPER(ch_staphylococcus, staphtyper_repeats, staphtyper_repeat_order)

    // Streptococcus
    ch_streptococcus = MERLINDIST.out.streptococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}
    ch_streptococcus_fq = MERLINDIST.out.streptococcus_fq.map{_meta, _reads, _found -> [_meta, _reads]}
    EMMTYPER(ch_streptococcus, emmtyper_blastdb)
    PBPTYPER(ch_streptococcus)
    SEROBA(ch_streptococcus_fq)
    SSUISSERO(ch_streptococcus)

    emit:
    results = MERLINDIST.out.results.mix(
        CLERMONTYPING.out.results,
        ECTYPER.out.results,
        EMMTYPER.out.results,
        GENOTYPHI.out.results,
        HICAP.out.results,
        HPSUISSERO.out.results,
        KLEBORATE.out.results,
        LEGSTA.out.results,
        LISSERO.out.results,
        MENINGOTYPE.out.results,
        NGMASTER.out.results,
        PASTY.out.results,
        PBPTYPER.out.results,
        SEQSERO2.out.results,
        SEROBA.out.results,
        SHIGAPASS.out.results,
        SHIGATYPER.out.results,
        SHIGEIFINDER.out.results,
        STECFINDER.out.results,
        SISTR.out.results,
        SSUISSERO.out.results,
        STAPHTYPER.out.results,
        TBPROFILER.out.results
    )
    logs = MERLINDIST.out.logs.mix(
        CLERMONTYPING.out.logs,
        ECTYPER.out.logs,
        EMMTYPER.out.logs,
        GENOTYPHI.out.logs,
        HICAP.out.logs,
        HPSUISSERO.out.logs,
        KLEBORATE.out.logs,
        LEGSTA.out.logs,
        LISSERO.out.logs,
        MENINGOTYPE.out.logs,
        NGMASTER.out.logs,
        PASTY.out.logs,
        PBPTYPER.out.logs,
        SEQSERO2.out.logs,
        SEROBA.out.logs,
        SHIGAPASS.out.logs,
        SHIGATYPER.out.logs,
        SHIGEIFINDER.out.logs,
        STECFINDER.out.logs,
        SISTR.out.logs,
        SSUISSERO.out.logs,
        STAPHTYPER.out.logs,
        TBPROFILER.out.logs
    )
    nf_logs = MERLINDIST.out.nf_logs.mix(
        CLERMONTYPING.out.nf_logs,
        ECTYPER.out.nf_logs,
        EMMTYPER.out.nf_logs,
        GENOTYPHI.out.nf_logs,
        HICAP.out.nf_logs,
        HPSUISSERO.out.nf_logs,
        KLEBORATE.out.nf_logs,
        LEGSTA.out.nf_logs,
        LISSERO.out.nf_logs,
        MENINGOTYPE.out.nf_logs,
        MERLINDIST.out.nf_logs,
        NGMASTER.out.nf_logs,
        PASTY.out.nf_logs,
        PBPTYPER.out.nf_logs,
        SEQSERO2.out.nf_logs,
        SEROBA.out.nf_logs,
        SHIGAPASS.out.nf_logs,
        SHIGATYPER.out.nf_logs,
        SHIGEIFINDER.out.nf_logs,
        STECFINDER.out.nf_logs,
        SISTR.out.nf_logs,
        SSUISSERO.out.nf_logs,
        STAPHTYPER.out.nf_logs,
        TBPROFILER.out.nf_logs
    )
    versions = MERLINDIST.out.versions.mix(
        CLERMONTYPING.out.versions,
        ECTYPER.out.versions,
        EMMTYPER.out.versions,
        GENOTYPHI.out.versions,
        HICAP.out.versions,
        HPSUISSERO.out.versions,
        KLEBORATE.out.versions,
        LEGSTA.out.versions,
        LISSERO.out.versions,
        MENINGOTYPE.out.versions,
        NGMASTER.out.versions,
        PASTY.out.versions,
        PBPTYPER.out.versions,
        SEQSERO2.out.versions,
        SEROBA.out.versions,
        SHIGAPASS.out.versions,
        SHIGATYPER.out.versions,
        SHIGEIFINDER.out.versions,
        STECFINDER.out.versions,
        SISTR.out.versions,
        SSUISSERO.out.versions,
        STAPHTYPER.out.versions,
        TBPROFILER.out.versions
    )
}
