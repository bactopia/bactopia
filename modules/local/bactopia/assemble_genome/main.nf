nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'assemble_genome')
options.ignore = ["-genome-size.txt", ".fastq.gz"]

process ASSEMBLE_GENOME {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${meta.id}"
    label "max_cpu_75"
    label "assemble_genome"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode,  overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(fq), path(extra), path(genome_size)

    output:
    tuple val(meta), path(genome_size), path("results/${meta.id}.{fna,fna.gz}"), file("total_contigs_*"), emit: fna, optional: true
    tuple val(meta), path("results/${meta.id}.{fna,fna.gz}"), emit: fna_only, optional: true
    tuple val(meta), path("results/${meta.id}.{fna,fna.gz}"), path(fq), emit: fna_fastq, optional: true
    tuple val(meta), path("blastdb/*"), emit: blastdb, optional: true
    path "results/*"
    path "${meta.id}-assembly-error.txt", optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    // Unicycler
    no_miniasm = params.no_miniasm ? "--no_miniasm" : ""
    no_rotate = params.no_rotate ? "--no_rotate" : ""
    keep = params.keep_all_files ? "--keep 3" : "--keep 1"
    is_hybrid = meta.runtype == "hybrid" ? "-l ${extra}" : ""
    unicycler_opts = "${no_miniasm} ${no_rotate} ${is_hybrid} ${keep}"

    // Shovill
    contig_namefmt = params.contig_namefmt ? params.contig_namefmt : "${meta.id}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0].toInteger()-1
    sopts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.no_stitch ? "--nostitch" : ""
    nocorr = params.no_corr ? "--nocorr" : ""
    shovill_mode = meta.single_end == false ? "shovill --R1 ${fq[0]} --R2 ${fq[1]} ${nostitch}" : "shovill-se --SE ${fq[0]}"
    shovill_opts = "--assembler ${params.shovill_assembler} --depth 0 --noreadcorr ${sopts} ${kmers} ${nocorr}"
    
    // Dragonflye
    dopts = params.dragonflye_opts ? "${params.dragonflye_opts}" : ""
    nopolish = params.no_polish ? "--nopolish" : ""
    medaka_model = params.medaka_model ? "--model ${params.medaka_model}" : ""
    pilon_rounds = params.pilon_rounds ? "--pilon ${params.pilon_rounds}" : ""
    dragonflye_fastq = meta.runtype == "short_polish" ? "--reads ${extra} --R1 ${fq[0]} --R2 ${fq[1]} --polypolish ${params.polypolish_rounds} ${pilon_rounds}" : "--reads ${fq[0]}"
    dragonflye_opts = "--assembler ${params.dragonflye_assembler} --depth 0 --minreadlen 0 --minquality 0 --racon ${params.racon_steps} --medaka ${params.medaka_steps} ${medaka_model} ${nopolish} ${dopts}"

    // Merge Shovill/Dragonflye opts
    assembler_wf = meta.runtype == "ont" || meta.runtype == "short_polish" ? "dragonflye" : "shovill"
    assembler_mode = meta.runtype == "ont" || meta.runtype == "short_polish" ? "dragonflye ${dragonflye_fastq}" : shovill_mode
    assembler_opts = meta.runtype == "ont" || meta.runtype == "short_polish" ? dragonflye_opts : shovill_opts

    // Assembly inputs
    use_original_assembly = null
    if (meta.runtype.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    '''
    OUTDIR=results
    GENOME_SIZE=`head -n 1 !{genome_size}`
    if [ "!{use_original_assembly}" == "true" ]; then
        mkdir ${OUTDIR}
        gzip -cd !{extra} > ${OUTDIR}/!{meta.id}.fna
    elif [[ "!{meta.runtype}" == "hybrid"  || "!{params.use_unicycler}" == "true" ]]; then
        unicycler -1 !{fq[0]} -2 !{fq[1]} -o ${OUTDIR}/ \
            --min_fasta_length !{params.min_contig_len} \
            --threads !{task.cpus} \
            --mode !{params.unicycler_mode} \
            --min_component_size !{params.min_component_size} \
            --min_dead_end_size !{params.min_dead_end_size} !{unicycler_opts}
        sed -r 's/^>([0-9]+)(.*)/>!{meta.id}_\\1\\2/' ${OUTDIR}/assembly.fasta > ${OUTDIR}/!{meta.id}.fna
        mv ${OUTDIR}/assembly.gfa ${OUTDIR}/unicycler.gfa
    else
        # Shovill or Dragonflye
        if ! !{assembler_mode} --gsize ${GENOME_SIZE} \
            --outdir ${OUTDIR} \
            --force \
            --minlen !{params.min_contig_len} \
            --mincov !{params.min_contig_cov} \
            --namefmt "!{contig_namefmt}" \
            --keepfiles \
            --cpus !{task.cpus} \
            --ram !{shovill_ram} !{assembler_opts}; then

            # Check if error is due to no contigs
            if grep "has zero contigs" results/!{assembler_wf}.log; then
                touch ${OUTDIR}/contigs.fa
            else
                exit 1
            fi
        fi

        mv ${OUTDIR}/contigs.fa ${OUTDIR}/!{meta.id}.fna

        # Rename Graphs
        if [ "!{assembler_wf}" == "shovill" ]; then 
            if [ -f "${OUTDIR}/contigs.gfa" ]; then
                mv ${OUTDIR}/contigs.gfa ${OUTDIR}/!{params.shovill_assembler}-unpolished.gfa
            elif [ -f "${OUTDIR}/contigs.fastg" ]; then
                mv ${OUTDIR}/contigs.fastg ${OUTDIR}/!{params.shovill_assembler}-unpolished.gfa
            elif [ -f "${OUTDIR}/contigs.LastGraph" ]; then
                mv ${OUTDIR}/contigs.LastGraph ${OUTDIR}/!{params.shovill_assembler}-unpolished.gfa
            fi
        fi

        if [ -f "${OUTDIR}/flye-info.txt" ]; then
            mv ${OUTDIR}/flye-info.txt ${OUTDIR}/flye.log
        fi
    fi


    # Check quality of assembly
    TOTAL_CONTIGS=`grep -c "^>" ${OUTDIR}/!{meta.id}.fna || true`
    touch "total_contigs_${TOTAL_CONTIGS}"
    if [ "${TOTAL_CONTIGS}" -gt 0 ]; then
        assembly-scan ${OUTDIR}/!{meta.id}.fna > ${OUTDIR}/!{meta.id}.json
        TOTAL_CONTIG_SIZE=$(grep "total_contig_length" ${OUTDIR}/!{meta.id}.json | sed -r 's/.*: ([0-9]+),/\\1/')
        if [ "${TOTAL_CONTIG_SIZE}" -lt !{params.min_genome_size} ]; then
            mv ${OUTDIR}/!{meta.id}.fna ${OUTDIR}/!{meta.id}-error.fna
            mv ${OUTDIR}/!{meta.id}.json ${OUTDIR}/!{meta.id}-error.json
            echo "!{meta.id} assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                    size (!{params.min_genome_size} bp). If this is unexpected, please investigate !{meta.id} to
                    determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of !{meta.id} will be discontinued." | \
            sed 's/^\\s*//' > !{meta.id}-assembly-error.txt
        else
            # Make BLASTDB
            mkdir blastdb
            cat ${OUTDIR}/!{meta.id}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}
        fi
    else
        mv ${OUTDIR}/!{meta.id}.fna ${OUTDIR}/!{meta.id}-error.fna
        echo "!{meta.id} assembled successfully, but 0 contigs were formed. Please investigate
                !{meta.id} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
                outcome. Further assembly-based analysis of !{meta.id} will be discontinued." | \
        sed 's/^\\s*//' > !{meta.id}-assembly-error.txt
    fi

    # Cleanup and compress
    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate files
        rm -rfv ${OUTDIR}/shovill.bam* ${OUTDIR}/shovill-se.bam* ${OUTDIR}/flash.extendedFrags*  \
                ${OUTDIR}/flash.notCombined* ${OUTDIR}/skesa.fasta* ${OUTDIR}/*.fq.gz ${OUTDIR}/00*.gfa \
                ${OUTDIR}/pilon_polish* ${OUTDIR}/flye/ ${OUTDIR}/flye.fasta* ${OUTDIR}/raven/  \
                ${OUTDIR}/raven.fasta* ${OUTDIR}/raven.cereal ${OUTDIR}/miniasm/ ${OUTDIR}/miniasm.fasta* \
                ${OUTDIR}/spades/ ${OUTDIR}/spades.fasta* ${OUTDIR}/megahit/ ${OUTDIR}/megahit.fasta* \
                ${OUTDIR}/velvet.fasta* ${OUTDIR}/velvet/
    fi

    if [[ !{params.skip_compression} == "false" ]]; then
        # Compress based on matched extensions
        find ${OUTDIR}/ -type f | \
            grep -E "\\.fna$|\\.fasta$|\\.fa$|\\.gfa$" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
    fi
    find ${OUTDIR} -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

    # Capture versions
    if [[ "$OSTYPE" == "darwin"* ]]; then
    
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        assembly-scan: $(echo $(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
        bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
        flash: $(echo $(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*$//')
        makeblastdb: $(echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
        megahit: $(echo $(megahit --version 2>&1) | sed 's/MEGAHIT v//')
        miniasm: $(echo $(miniasm -V))
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
        pilon: $(echo $(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*$//')
        racon: $(echo $(racon --version 2>&1) | sed 's/v//')
        samclip: $(echo $(samclip --version 2>&1) | sed 's/samclip //')
        samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
        shovill: $(echo $(shovill --version 2>&1) | sed 's/shovill //')
        shovill-se: $(echo $(shovill-se --version 2>&1) | sed 's/shovill-se //')
        skesa: $(echo $(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*$//')
        spades.py: $(echo $(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
        velvetg: $(echo $(velvetg 2>&1) | sed 's/^.*Version //;s/ .*$//')
        velveth: $(echo $(velveth 2>&1) | sed 's/^.*Version //;s/ .*$//')
        unicycler: $(echo $(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*$//')
    END_VERSIONS
    
    else
    
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        any2fasta:  $(echo $(any2fasta -v 2>&1) | sed 's/any2fasta //')
        assembly-scan: $(echo $(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
        bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
        flash: $(echo $(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*$//')
        flye: $(echo $(flye --version))
        makeblastdb: $(echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
        medaka: $(echo $(medaka --version 2>&1) | sed 's/medaka //')
        megahit: $(echo $(megahit --version 2>&1) | sed 's/MEGAHIT v//')
        miniasm: $(echo $(miniasm -V))
        minimap2: $(echo $(minimap2 --version))
        nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
        pilon: $(echo $(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*$//')
        racon: $(echo $(racon --version 2>&1) | sed 's/v//')
        rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
        raven: $(echo $(raven --version))
        samclip: $(echo $(samclip --version 2>&1) | sed 's/samclip //')
        samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
        shovill: $(echo $(shovill --version 2>&1) | sed 's/shovill //')
        shovill-se: $(echo $(shovill-se --version 2>&1) | sed 's/shovill-se //')
        skesa: $(echo $(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*$//')
        spades.py: $(echo $(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
        velvetg: $(echo $(velvetg 2>&1) | sed 's/^.*Version //;s/ .*$//')
        velveth: $(echo $(velveth 2>&1) | sed 's/^.*Version //;s/ .*$//')
        unicycler: $(echo $(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*$//')
    END_VERSIONS
    
    fi
    '''
}
