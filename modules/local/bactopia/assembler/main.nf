// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options        = initOptions(params.options ? params.options : [:], 'assembler')
options.ignore = [".fastq.gz"]
options.btype  = "main"
conda_tools    = "bioconda::bactopia-assembler=1.0.4"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ASSEMBLER {
    tag "${meta.id}"
    label (params.use_unicycler ? "process_medium" : "process_low")

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-assembler:1.0.4--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-assembler:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fq), path(extra)

    output:
    tuple val(meta), path("results/${prefix}.{fna,fna.gz}"), path("fastqs/${prefix}*.fastq.gz"), emit: fna_fastq, optional: true
    tuple val(meta), path("results/${prefix}.{fna,fna.gz}"), emit: fna, optional: true
    tuple val(meta), path("results/${prefix}.tsv")         , emit: tsv, optional: true
    path "results/*"                                       , emit: results
    path "${prefix}-assembly-error.txt"                    , emit: error, optional: true
    path "*.{log,err}"                                     , emit: logs, optional: true
    path ".command.*"                                      , emit: nf_logs
    path "versions.yml"                                    , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"

    // Determine input reads
    r1 = null
    r2 = null
    se = null
    if (fq[2]) {
        r1 = fq[0].getName().endsWith('_R1.fastq.gz') ? fq[0] : fq[1].getName().endsWith('_R1.fastq.gz') ? fq[1] : fq[2]
        r2 = fq[0].getName().endsWith('_R2.fastq.gz') ? fq[0] : fq[1].getName().endsWith('_R2.fastq.gz') ? fq[1] : fq[2]
        se = !fq[0].getName().matches('.*_R[12].fastq.gz') ? fq[0] : !fq[1].getName().matches('.*_R[12].fastq.gz') ? fq[1] : fq[2]
    } else if (fq[1]) {
        r1 = fq[0]
        r2 = fq[1]
    } else {
        se = fq[0]
    }

    // Unicycler
    is_hybrid = meta.runtype == "hybrid" ? "-l ${se}" : ""

    // Shovill
    contig_namefmt = params.contig_namefmt ? params.contig_namefmt : "${prefix}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0].toInteger()-1
    shovill_mode = meta.single_end == false ? "shovill --R1 ${r1} --R2 ${r2}" : "shovill-se --SE ${se}"

    // Dragonflye
    dragonflye_fastq = meta.runtype == "short_polish" ? "--reads ${se} --R1 ${r1} --R2 ${r2}" : "--reads ${se}"

    // Assembly inputs
    use_original_assembly = null
    if (meta.runtype.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    """
    echo "R1 ${r1}"
    echo "R2 ${r2}"
    echo "SE ${se}"

    if [ "${use_original_assembly}" == "true" ]; then
        mkdir results
        gzip -cd ${extra} > results/${prefix}.fna
    elif [[ "${meta.runtype}" == "hybrid" || "${params.use_unicycler}" == "true" ]]; then
        # Unicycler
        unicycler \\
            -1 ${r1} -2 ${r2} \\
            ${options.args3} \\
            ${is_hybrid} \\
            -o results/ \\
            --threads ${task.cpus}
        sed -r 's/^>([0-9]+)(.*)/>${prefix}_\\1\\2/' results/assembly.fasta > results/${prefix}.fna
        mv results/assembly.fasta results/unicycler-unpolished.fasta
        mv results/assembly.gfa results/unicycler-unpolished.gfa
    elif [[ "${meta.runtype}" == "ont" || "${meta.runtype}" == "short_polish" ]]; then
        # Dragonflye
        if ! dragonflye \\
            ${dragonflye_fastq} \\
            --gsize ${meta.genome_size} \\
            --outdir results \\
            ${options.args2} \\
            --namefmt "${contig_namefmt}" \\
            --cpus ${task.cpus} \\
            --ram ${shovill_ram} \\
            --noreorient; then

            # Check if error is due to no contigs
            if grep "has zero contigs" results/dragonflye.log; then
                touch results/contigs.fa
            else
                exit 1
            fi
        fi
        mv results/contigs.fa results/${prefix}.fna
    else
        # Shovill
        if ! ${shovill_mode} \\
            --gsize ${meta.genome_size} \\
            --outdir results \\
            ${options.args} \\
            --namefmt "${contig_namefmt}" \\
            --cpus ${task.cpus} \\
            --ram ${shovill_ram}; then

            # Check if error is due to no contigs
            if grep "has zero contigs" results/shovill.log; then
                touch results/contigs.fa
            else
                exit 1
            fi
        fi
        mv results/contigs.fa results/${prefix}.fna

        # Rename Graphs
        if [ -f "results/contigs.gfa" ]; then
            mv results/contigs.gfa results/${params.shovill_assembler}-unpolished.gfa
        elif [ -f "results/contigs.fastg" ]; then
            mv results/contigs.fastg results/${params.shovill_assembler}-unpolished.gfa
        elif [ -f "results/contigs.LastGraph" ]; then
            mv results/contigs.LastGraph results/${params.shovill_assembler}-unpolished.gfa
        fi

        if [ -f "results/flye-info.txt" ]; then
            mv results/flye-info.txt results/flye.log
        fi
    fi

    # Check quality of assembly
    TOTAL_CONTIGS=`grep -c "^>" results/${prefix}.fna || true`
    if [ "\${TOTAL_CONTIGS}" -gt 0 ]; then
        assembly-scan results/${prefix}.fna --prefix ${prefix} > results/${prefix}.tsv
        TOTAL_CONTIG_SIZE=\$(cut -f 3 results/${prefix}.tsv | tail -n 1)
        if [ "\${TOTAL_CONTIG_SIZE}" -lt ${params.min_genome_size} ]; then
            mv results/${prefix}.fna results/${prefix}-error.fna
            mv results/${prefix}.tsv results/${prefix}-error.tsv
            echo "${prefix} assembled size (\${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                    size (${params.min_genome_size} bp). If this is unexpected, please investigate ${prefix} to
                    determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of ${prefix} will be discontinued." | \
            sed 's/^\\s*//' > ${prefix}-assembly-error.txt
        fi
    else
        mv results/${prefix}.fna results/${prefix}-error.fna
        echo "${prefix} assembled successfully, but 0 contigs were formed. Please investigate
                ${prefix} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
                outcome. Further assembly-based analysis of ${prefix} will be discontinued." | \
        sed 's/^\\s*//' > ${prefix}-assembly-error.txt
    fi

    # Cleanup and compress
    if [ "${params.keep_all_files}" == "false" ]; then
        # Remove intermediate files
        rm -rfv results/shovill.bam* \\
                results/shovill-se.bam* \\
                results/flash.extendedFrags* \\
                results/flash.notCombined* \\
                results/skesa.fasta* \\
                results/*.fq.gz \\
                results/00*.gfa \\
                results/pilon_polish* \\
                results/flye/ \\
                results/flye.fasta* \\
                results/raven/ \\
                results/raven.fasta* \\
                results/raven.cereal \\
                results/miniasm/ \\
                results/miniasm.fasta* \\
                results/spades/ \\
                results/spades.fasta* \\
                results/megahit/ \\
                results/megahit.fasta* \\
                results/velvet.fasta* \\
                results/velvet/
    fi

    if [[ ${params.skip_compression} == "false" ]]; then
        # Compress based on matched extensions
        find results/ -type f | \
            grep -E "\\.fna\$|\\.fasta\$|\\.fa\$|\\.gfa\$" | \
            xargs -I {} pigz -n --best -p ${task.cpus} {}
    fi
    find results -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

    # Move primary fastqs
    mkdir fastqs/
    if [ "${meta.runtype}" == "hybrid" ]; then
        # Used Unicycler, so Illumina reads are expected to be best
        cp ${r1} fastqs/
        cp ${r2} fastqs/
    elif [[ "${meta.runtype}" == "ont" || "${meta.runtype}" == "short_polish" || "${meta.single_end}" == "true" ]]; then
        # Used Dragonflye or Illumina single-end reads, so use these reads going forward
        cp ${se} fastqs/
    else
        # Illumina Pair-end reads
        cp ${r1} fastqs/
        cp ${r2} fastqs/
    fi

    # Capture versions
    if [[ "\$OSTYPE" == "darwin"* ]]; then
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assembly-scan: \$(echo \$(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //;s/ .*\$//')
        flash: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*\$//')
        megahit: \$(echo \$(megahit --version 2>&1) | sed 's/MEGAHIT v//')
        miniasm: \$(echo \$(miniasm -V))
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        pilon: \$(echo \$(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*\$//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/v//')
        samclip: \$(echo \$(samclip --version 2>&1) | sed 's/^.*samclip //')
        samtools: \$(echo \$(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*\$//')
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
        shovill-se: \$(echo \$(shovill-se --version 2>&1) | sed 's/^.*shovill-se //')
        skesa: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*\$//')
        spades.py: \$(echo \$(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
        velvetg: \$(echo \$(velvetg 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        velveth: \$(echo \$(velveth 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*\$//')
    END_VERSIONS
    
    else
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        any2fasta: \$(echo \$(any2fasta -v 2>&1) | sed 's/^.*any2fasta //')
        assembly-scan: \$(echo \$(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //;s/ .*\$//')
        dragonflye: \$(echo \$(dragonflye --version 2>&1) | sed 's/^.*dragonflye //' )
        flash: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*\$//')
        flye: \$(echo \$(flye --version))
        medaka: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
        megahit: \$(echo \$(megahit --version 2>&1) | sed 's/MEGAHIT v//')
        miniasm: \$(echo \$(miniasm -V))
        minimap2: \$(echo \$(minimap2 --version))
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        pilon: \$(echo \$(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*\$//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/v//')
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
        raven: \$(echo \$(raven --version))
        samclip: \$(echo \$(samclip --version 2>&1) | sed 's/^.*samclip //')
        samtools: \$(echo \$(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*\$//')
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
        shovill-se: \$(echo \$(shovill-se --version 2>&1) | sed 's/^.*shovill-se //')
        skesa: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*\$//')
        spades.py: \$(echo \$(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
        velvetg: \$(echo \$(velvetg 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        velveth: \$(echo \$(velveth 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*\$//')
    END_VERSIONS
    
    fi
    """
}
