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
 * - `meta`: Groovy Record containing sample information
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
    ch_merlindist = MERLINDIST(assembly, mash_db)

    // Helper closures to build records from MERLINDIST output
    def forAssembly = { r -> record(meta: r.meta, fna: r.fna) }
    def forReads = { r -> record(meta: r.meta, r1: r.r1, r2: r.r2, se: r.se, lr: r.lr) }
    def forSeqs = { r -> record(meta: r.meta, fna: r.fna, r1: r.r1, r2: r.r2, se: r.se, lr: r.lr) }

    // Escherichia/Shigella
    ch_escherichia = ch_merlindist.sample_outputs.filter { r -> r.escherichia != null }
    ch_clermontyping = CLERMONTYPING(ch_escherichia.map(forAssembly))
    ch_ectyper = ECTYPER(ch_escherichia.map(forAssembly))
    ch_shigapass = SHIGAPASS(ch_escherichia.map(forAssembly))
    ch_shigatyper = SHIGATYPER(ch_escherichia.map(forReads))
    ch_shigeifinder = SHIGEIFINDER(ch_escherichia.map(forAssembly))
    ch_stecfinder = STECFINDER(ch_escherichia.map(forSeqs))

    // Haemophilus
    ch_haemophilus = ch_merlindist.sample_outputs.filter { r -> r.haemophilus != null }
    ch_hicap = HICAP(ch_haemophilus.map(forAssembly), hicap_database_dir, hicap_model_fp)
    ch_hpsuissero = HPSUISSERO(ch_haemophilus.map(forAssembly))

    // Klebsiella
    ch_klebsiella = ch_merlindist.sample_outputs.filter { r -> r.klebsiella != null }
    ch_kleborate = KLEBORATE(ch_klebsiella.map(forAssembly))

    // Legionella
    ch_legionella = ch_merlindist.sample_outputs.filter { r -> r.legionella != null }
    ch_legsta = LEGSTA(ch_legionella.map(forAssembly))

    // Listeria
    ch_listeria = ch_merlindist.sample_outputs.filter { r -> r.listeria != null }
    ch_lissero = LISSERO(ch_listeria.map(forAssembly))

    // Mycobacterium
    ch_mycobacterium = ch_merlindist.sample_outputs.filter { r -> r.mycobacterium != null }
    ch_tbprofiler = TBPROFILER(ch_mycobacterium.map(forReads))

    // Neisseria
    ch_neisseria = ch_merlindist.sample_outputs.filter { r -> r.neisseria != null }
    ch_ngmaster = NGMASTER(ch_neisseria.map(forAssembly))

    // Pseudomonas
    ch_pseudomonas = ch_merlindist.sample_outputs.filter { r -> r.pseudomonas != null }
    ch_pasty = PASTY(ch_pseudomonas.map(forAssembly))

    // Salmonella
    ch_salmonella = ch_merlindist.sample_outputs.filter { r -> r.salmonella != null }
    ch_genotyphi = GENOTYPHI(ch_salmonella.map(forReads))
    ch_seqsero2 = SEQSERO2(ch_salmonella.map(forAssembly))
    ch_sistr = SISTR(ch_salmonella.map(forAssembly))

    // Staphylococcus
    ch_staphylococcus = ch_merlindist.sample_outputs.filter { r -> r.staphylococcus != null }
    ch_staphtyper = STAPHTYPER(ch_staphylococcus.map(forAssembly), staphtyper_repeats, staphtyper_repeat_order)

    // Streptococcus
    ch_streptococcus = ch_merlindist.sample_outputs.filter { r -> r.streptococcus != null }
    ch_emmtyper = EMMTYPER(ch_streptococcus.map(forAssembly), emmtyper_blastdb)
    ch_pbptyper = PBPTYPER(ch_streptococcus.map(forAssembly))
    ch_seroba = SEROBA(ch_streptococcus.map(forReads))
    ch_ssuissero = SSUISSERO(ch_streptococcus.map(forAssembly))

    emit:
    // Published outputs
    sample_outputs = ch_merlindist.sample_outputs.mix(
        ch_clermontyping.sample_outputs,
        ch_ectyper.sample_outputs,
        ch_emmtyper.sample_outputs,
        ch_genotyphi.sample_outputs,
        ch_hicap.sample_outputs,
        ch_hpsuissero.sample_outputs,
        ch_kleborate.sample_outputs,
        ch_legsta.sample_outputs,
        ch_lissero.sample_outputs,
        ch_ngmaster.sample_outputs,
        ch_pasty.sample_outputs,
        ch_pbptyper.sample_outputs,
        ch_seqsero2.sample_outputs,
        ch_seroba.sample_outputs,
        ch_shigapass.sample_outputs,
        ch_shigatyper.sample_outputs,
        ch_shigeifinder.sample_outputs,
        ch_stecfinder.sample_outputs,
        ch_sistr.sample_outputs,
        ch_ssuissero.sample_outputs,
        ch_staphtyper.sample_outputs,
        ch_tbprofiler.sample_outputs
    )
    run_outputs = ch_merlindist.run_outputs.mix(
        ch_clermontyping.run_outputs,
        ch_ectyper.run_outputs,
        ch_emmtyper.run_outputs,
        ch_genotyphi.run_outputs,
        ch_hicap.run_outputs,
        ch_hpsuissero.run_outputs,
        ch_kleborate.run_outputs,
        ch_legsta.run_outputs,
        ch_lissero.run_outputs,
        ch_ngmaster.run_outputs,
        ch_pasty.run_outputs,
        ch_pbptyper.run_outputs,
        ch_seqsero2.run_outputs,
        ch_seroba.run_outputs,
        ch_shigapass.run_outputs,
        ch_shigatyper.run_outputs,
        ch_shigeifinder.run_outputs,
        ch_stecfinder.run_outputs,
        ch_sistr.run_outputs,
        ch_ssuissero.run_outputs,
        ch_staphtyper.run_outputs,
        ch_tbprofiler.run_outputs
    )
}
