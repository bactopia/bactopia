process ASSEMBLER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fq), path(extra)

    output:
    tuple val(meta), path("${prefix}.{fna,fna.gz}")      , emit: fna, optional: true
    tuple val(meta), path("${prefix}.tsv")               , emit: tsv, optional: true
    tuple val(meta), path("supplemental/*")              , emit: supplemental
    tuple val(meta), path("${prefix}-assembly-error.txt"), emit: error, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.runtype = _meta.runtype
    meta.genome_size = _meta.genome_size
    meta.species = _meta.species
    meta.single_end = _meta.single_end

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
    contig_namefmt = task.ext.contig_namefmt ? task.ext.contig_namefmt : "${prefix}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0].toInteger()-1
    shovill_mode = meta.single_end == false ? "shovill --R1 ${r1} --R2 ${r2}" : "shovill-se --SE ${se}"

    // Dragonflye
    dragonflye_fastq = meta.runtype == "short_polish" ? "--reads ${se} --R1 ${r1} --R2 ${r2}" : "--reads ${se}"

    // Assembly inputs
    use_original_assembly = null
    if (meta.runtype.startsWith('assembly')) {
        use_original_assembly = task.ext.reassemble ? false : true
    }
    """
    echo "R1 ${r1}"
    echo "R2 ${r2}"
    echo "SE ${se}"

    if [ "${use_original_assembly}" == "true" ]; then
        mkdir supplemental
        gzip -cd ${extra} > supplemental/${prefix}.fna
    elif [[ "${meta.runtype}" == "hybrid" || "${task.ext.use_unicycler}" == "true" ]]; then
        # Unicycler
        unicycler \\
            -1 ${r1} -2 ${r2} \\
            ${task.ext.args3} \\
            ${is_hybrid} \\
            -o supplemental/ \\
            --threads ${task.cpus}
        sed -r 's/^>([0-9]+)(.*)/>${prefix}_\\1\\2/' supplemental/assembly.fasta > supplemental/${prefix}.fna
        mv supplemental/assembly.fasta supplemental/unicycler-unpolished.fasta
        mv supplemental/assembly.gfa supplemental/unicycler-unpolished.gfa
    elif [[ "${meta.runtype}" == "ont" || "${meta.runtype}" == "short_polish" ]]; then
        # Dragonflye
        if ! dragonflye \\
            ${dragonflye_fastq} \\
            --gsize ${meta.genome_size} \\
            --outdir supplemental \\
            ${task.ext.args2} \\
            --namefmt "${contig_namefmt}" \\
            --cpus ${task.cpus} \\
            --ram ${shovill_ram} \\
            --noreorient; then

            # Check if error is due to no contigs
            if grep "has zero contigs" supplemental/dragonflye.log; then
                touch supplemental/contigs.fa
            else
                exit 1
            fi
        fi
        mv supplemental/contigs.fa supplemental/${prefix}.fna
    else
        # Shovill
        if ! ${shovill_mode} \\
            --gsize ${meta.genome_size} \\
            --outdir supplemental \\
            ${task.ext.args} \\
            --namefmt "${contig_namefmt}" \\
            --cpus ${task.cpus} \\
            --ram ${shovill_ram}; then

            # Check if error is due to no contigs
            if grep "has zero contigs" supplemental/shovill.log; then
                touch supplemental/contigs.fa
            else
                exit 1
            fi
        fi
        mv supplemental/contigs.fa supplemental/${prefix}.fna

        # Rename Graphs
        if [ -f "supplemental/contigs.gfa" ]; then
            mv supplemental/contigs.gfa supplemental/${task.ext.shovill_assembler}-unpolished.gfa
        elif [ -f "supplemental/contigs.fastg" ]; then
            mv supplemental/contigs.fastg supplemental/${task.ext.shovill_assembler}-unpolished.gfa
        elif [ -f "supplemental/contigs.LastGraph" ]; then
            mv supplemental/contigs.LastGraph supplemental/${task.ext.shovill_assembler}-unpolished.gfa
        fi

        if [ -f "supplemental/flye-info.txt" ]; then
            mv supplemental/flye-info.txt supplemental/flye.log
        fi
    fi

    # Check quality of assembly
    TOTAL_CONTIGS=`grep -c "^>" supplemental/${prefix}.fna || true`
    if [ "\${TOTAL_CONTIGS}" -gt 0 ]; then
        assembly-scan supplemental/${prefix}.fna --prefix ${prefix} > supplemental/${prefix}.tsv
        TOTAL_CONTIG_SIZE=\$(cut -f 3 supplemental/${prefix}.tsv | tail -n 1)
        if [ "\${TOTAL_CONTIG_SIZE}" -lt ${task.ext.min_genome_size} ]; then
            mv supplemental/${prefix}.fna supplemental/${prefix}-error.fna
            mv supplemental/${prefix}.tsv supplemental/${prefix}-error.tsv
            echo "${prefix} assembled size [\${TOTAL_CONTIG_SIZE} bp] is less than the minimum allowed genome
                    size [${task.ext.min_genome_size} bp]. If this is unexpected, please investigate ${prefix} to
                    determine a cause [e.g. metagenomic, contaminants, etc...] for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of ${prefix} will be discontinued." | \
            sed 's/^\\s*//' > ${prefix}-assembly-error.txt
        fi
    else
        mv supplemental/${prefix}.fna supplemental/${prefix}-error.fna
        echo "${prefix} assembled successfully, but 0 contigs were formed. Please investigate
                ${prefix} to determine a cause [e.g. metagenomic, contaminants, etc...] for this
                outcome. Further assembly-based analysis of ${prefix} will be discontinued." | \
        sed 's/^\\s*//' > ${prefix}-assembly-error.txt
    fi

    # Cleanup and compress
    if [ "${task.ext.keep_all_files}" == "false" ]; then
        # Remove intermediate files
        rm -rfv supplemental/shovill.bam* \\
                supplemental/shovill-se.bam* \\
                supplemental/flash.extendedFrags* \\
                supplemental/flash.notCombined* \\
                supplemental/skesa.fasta* \\
                supplemental/*.fq.gz \\
                supplemental/00*.gfa \\
                supplemental/pilon_polish* \\
                supplemental/flye/ \\
                supplemental/flye.fasta* \\
                supplemental/raven/ \\
                supplemental/raven.fasta* \\
                supplemental/raven.cereal \\
                supplemental/miniasm/ \\
                supplemental/miniasm.fasta* \\
                supplemental/spades/ \\
                supplemental/spades.fasta* \\
                supplemental/megahit/ \\
                supplemental/megahit.fasta* \\
                supplemental/velvet.fasta* \\
                supplemental/velvet/
    fi

    if [[ "${task.ext.skip_compression}" == "false" ]]; then
        # Compress based on matched extensions
        find supplemental/ -type f | \
            grep -E "\\.fna\$|\\.fasta\$|\\.fa\$|\\.gfa\$" | \
            xargs -I {} pigz -n --best -p ${task.cpus} {}
    fi
    find supplemental/ -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

    if [ -f "supplemental/${prefix}.tsv" ]; then
        mv supplemental/${prefix}.tsv ./
    fi

    if [ -f "supplemental/${prefix}.fna" ]; then
        mv supplemental/${prefix}.fna ./
    fi

    if [ -f "supplemental/${prefix}.fna.gz" ]; then
        mv supplemental/${prefix}.fna.gz ./
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
