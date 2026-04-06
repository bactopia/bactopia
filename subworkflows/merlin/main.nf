/**
 * MinER assisted species-specific bactopia tool seLectIoN.
 *
 * This subworkflow performs intelligent species identification and selects appropriate
 * species-specific typing tools based on the detected organism. It first identifies
 * potential species using MinHash distance estimation, then runs species-specific
 * subworkflows for detailed characterization including serotyping, MLST, virulence
 * factor detection, and antimicrobial resistance profiling.
 *
 * @status stable
 * @keywords species, identification, typing, serotype, virulence
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic,components
 * @citation mash
 *
 * @subworkflows merlindist, clermontyping, ectyper, emmtyper, genotyphi, hicap, hpsuissero, kleborate,
 *              legsta, lissero, ngmaster, pasty, pbptyper, seqsero2, seroba, shigapass,
 *              shigatyper, shigeifinder, sistr, ssuissero, staphtyper, stecfinder, tbprofiler
 *
 * @input record(meta, fna, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembly file for species identification and typing
 * - `r1?`: Illumina R1 reads (paired-end) or null
 * - `r2?`: Illumina R2 reads (paired-end) or null
 * - `se?`: Single-end Illumina reads or null
 * - `lr?`: Long reads (ONT/PacBio) or null
 *
 * @input mash_db
 * Mash sketch database for rapid species identification
 *
 * @input emmtyper_blastdb
 * EMMTyper BLAST database for Streptococcus pyogenes emm typing (optional)
 *
 * @input hicap_database_dir
 * HiCAP database directory for Haemophilus influenzae serotyping (optional)
 *
 * @input hicap_model_fp
 * HiCAP HMM model file for improved detection (optional)
 *
 * @input staphtyper_repeats
 * Staphylococcus aureus repeat sequences for spa typing (optional)
 *
 * @input staphtyper_repeat_order
 * Staphylococcus aureus repeat order file for spa typing (optional)
 *
 * @output sample_outputs
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { MERLINDIST    } from '../merlindist/main';
include { CLERMONTYPING } from '../clermontyping/main';
include { ECTYPER       } from '../ectyper/main';
include { EMMTYPER      } from '../emmtyper/main';
include { HICAP         } from '../hicap/main';
include { HPSUISSERO    } from '../hpsuissero/main';
include { GENOTYPHI     } from '../genotyphi/main';
include { KLEBORATE     } from '../kleborate/main';
include { LEGSTA        } from '../legsta/main';
include { LISSERO       } from '../lissero/main';
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

workflow MERLIN {
    take:
    assembly: Channel<Record>
    mash_db: Value<Path>
    emmtyper_blastdb: Value<Path?>
    hicap_database_dir: Value<Path?>
    hicap_model_fp: Value<Path?>
    staphtyper_repeats: Value<Path?>
    staphtyper_repeat_order: Value<Path?>

    main:
    // ID potential species
    MERLINDIST(assembly, mash_db)

    // Helper closures to build records from MERLINDIST output
    def forAssembly = { r -> record(meta: r.meta, fna: r.fna) }
    def forReads = { r -> record(meta: r.meta, r1: r.r1, r2: r.r2, se: r.se, lr: r.lr) }
    def forSeqs = { r -> record(meta: r.meta, fna: r.fna, r1: r.r1, r2: r.r2, se: r.se, lr: r.lr) }

    // Escherichia/Shigella
    ch_escherichia = MERLINDIST.out.sample_outputs.filter { r -> r.escherichia != null }
    CLERMONTYPING(ch_escherichia.map(forAssembly))
    ECTYPER(ch_escherichia.map(forAssembly))
    SHIGAPASS(ch_escherichia.map(forAssembly))
    SHIGATYPER(ch_escherichia.map(forReads))
    SHIGEIFINDER(ch_escherichia.map(forAssembly))
    STECFINDER(ch_escherichia.map(forSeqs))

    // Haemophilus
    ch_haemophilus = MERLINDIST.out.sample_outputs.filter { r -> r.haemophilus != null }
    HICAP(ch_haemophilus.map(forAssembly), hicap_database_dir, hicap_model_fp)
    HPSUISSERO(ch_haemophilus.map(forAssembly))

    // Klebsiella
    ch_klebsiella = MERLINDIST.out.sample_outputs.filter { r -> r.klebsiella != null }
    KLEBORATE(ch_klebsiella.map(forAssembly))

    // Legionella
    ch_legionella = MERLINDIST.out.sample_outputs.filter { r -> r.legionella != null }
    LEGSTA(ch_legionella.map(forAssembly))

    // Listeria
    ch_listeria = MERLINDIST.out.sample_outputs.filter { r -> r.listeria != null }
    LISSERO(ch_listeria.map(forAssembly))

    // Mycobacterium
    ch_mycobacterium = MERLINDIST.out.sample_outputs.filter { r -> r.mycobacterium != null }
    TBPROFILER(ch_mycobacterium.map(forReads))

    // Neisseria
    ch_neisseria = MERLINDIST.out.sample_outputs.filter { r -> r.neisseria != null }
    NGMASTER(ch_neisseria.map(forAssembly))

    // Pseudomonas
    ch_pseudomonas = MERLINDIST.out.sample_outputs.filter { r -> r.pseudomonas != null }
    PASTY(ch_pseudomonas.map(forAssembly))

    // Salmonella
    ch_salmonella = MERLINDIST.out.sample_outputs.filter { r -> r.salmonella != null }
    GENOTYPHI(ch_salmonella.map(forReads))
    SEQSERO2(ch_salmonella.map(forAssembly))
    SISTR(ch_salmonella.map(forAssembly))

    // Staphylococcus
    ch_staphylococcus = MERLINDIST.out.sample_outputs.filter { r -> r.staphylococcus != null }
    STAPHTYPER(ch_staphylococcus.map(forAssembly), staphtyper_repeats, staphtyper_repeat_order)

    // Streptococcus
    ch_streptococcus = MERLINDIST.out.sample_outputs.filter { r -> r.streptococcus != null }
    EMMTYPER(ch_streptococcus.map(forAssembly), emmtyper_blastdb)
    PBPTYPER(ch_streptococcus.map(forAssembly))
    SEROBA(ch_streptococcus.map(forReads))
    SSUISSERO(ch_streptococcus.map(forAssembly))

    emit:
    // Published outputs
    sample_outputs = MERLINDIST.out.sample_outputs.mix(
        CLERMONTYPING.out.sample_outputs,
        ECTYPER.out.sample_outputs,
        EMMTYPER.out.sample_outputs,
        GENOTYPHI.out.sample_outputs,
        HICAP.out.sample_outputs,
        HPSUISSERO.out.sample_outputs,
        KLEBORATE.out.sample_outputs,
        LEGSTA.out.sample_outputs,
        LISSERO.out.sample_outputs,
        NGMASTER.out.sample_outputs,
        PASTY.out.sample_outputs,
        PBPTYPER.out.sample_outputs,
        SEQSERO2.out.sample_outputs,
        SEROBA.out.sample_outputs,
        SHIGAPASS.out.sample_outputs,
        SHIGATYPER.out.sample_outputs,
        SHIGEIFINDER.out.sample_outputs,
        STECFINDER.out.sample_outputs,
        SISTR.out.sample_outputs,
        SSUISSERO.out.sample_outputs,
        STAPHTYPER.out.sample_outputs,
        TBPROFILER.out.sample_outputs
    )
    run_outputs = MERLINDIST.out.run_outputs.mix(
        CLERMONTYPING.out.run_outputs,
        ECTYPER.out.run_outputs,
        EMMTYPER.out.run_outputs,
        GENOTYPHI.out.run_outputs,
        HICAP.out.run_outputs,
        HPSUISSERO.out.run_outputs,
        KLEBORATE.out.run_outputs,
        LEGSTA.out.run_outputs,
        LISSERO.out.run_outputs,
        NGMASTER.out.run_outputs,
        PASTY.out.run_outputs,
        PBPTYPER.out.run_outputs,
        SEQSERO2.out.run_outputs,
        SEROBA.out.run_outputs,
        SHIGAPASS.out.run_outputs,
        SHIGATYPER.out.run_outputs,
        SHIGEIFINDER.out.run_outputs,
        STECFINDER.out.run_outputs,
        SISTR.out.run_outputs,
        SSUISSERO.out.run_outputs,
        STAPHTYPER.out.run_outputs,
        TBPROFILER.out.run_outputs
    )
}
