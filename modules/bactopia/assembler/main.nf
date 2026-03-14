/**
 * Assemble bacterial genomes using short read, long read, or hybrid strategies.
 *
 * Automatically selects the appropriate assembler based on input read types:
 * - **Short Paired-End Reads:** Uses [Shovill](https://github.com/tseemann/shovill) (SKESA/SPAdes wrapper).
 * - **Short Single-End Reads:** Uses [Shovill](https://github.com/rpetit3/shovill) (SKESA/SPAdes wrapper).
 * - **Long Reads:** Uses [Dragonflye](https://github.com/rpetit3/dragonflye) (Flye/Miniasm wrapper).
 * - **Hybrid:** Uses [Unicycler](https://github.com/rrwick/Unicycler) or Dragonflye (with polishing).
 *
 * Summary statistics for each assembly are generated using [assembly-scan](https://github.com/rpetit3/assembly-scan).
 *
 * Uses named record input with explicit read slots (r1, r2, se, lr, assembly) as Path?.
 *
 * @status stable
 * @keywords bacteria, assembly, hybrid, shovill, dragonflye, unicycler, illumina, nanopore
 * @tags complexity:complex input-type:multiple output-type:multiple features:conditional-logic,alternative-execution
 * @citation any2fasta, assembly_scan, bwa, dragonflye, flash, flye, medaka, megahit, miniasm, minimap2, nanoq, pigz, pilon, racon, rasusa, raven, samclip, samtools, shovill, shovill_se, skesa, spades, unicycler, velvet
 *
 * @note When runtype is 'assembly' or 'assembly_accession' and --reassemble is not set,
 * the original assembly is used without re-assembly.
 *
 * @input record(meta, r1, r2, se, lr, assembly)
 * - `meta`    : Groovy Map containing sample information
 * - `r1`      : Illumina R1 reads (paired-end)
 * - `r2`      : Illumina R2 reads (paired-end)
 * - `se`      : Single-end Illumina reads
 * - `lr`      : Long reads (ONT/PacBio) for long-read or hybrid assembly
 * - `assembly`: Assembly file (FASTA) for assembly-based runtypes
 *
 * @output record(meta, assembly, r1, r2, se, lr, tsv, supplemental, error, results, logs, nf_logs, versions)
 * - `assembly`: Assembled contigs in FASTA format
 * - `r1`: Passthrough Illumina R1 reads
 * - `r2`: Passthrough Illumina R2 reads
 * - `se`: Passthrough single-end reads
 * - `lr`: Passthrough long reads
 * - `tsv`: Tab-delimited report of assembly statistics (N50, length, coverage)
 * - `supplemental`: Supplemental files including assembly graphs (*.gfa) and tool-specific logs
 * - `error`: Captured error messages if assembly fails
 */
nextflow.preview.types = true

process ASSEMBLER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?, assembly: Path?): Record

    stage:
    stageAs 'input-assembly/*', assembly

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        assembly: file("${prefix}.{fna,fna.gz}", optional: true),
        r1: file(r1),
        r2: file(r2),
        se: file(se),
        lr: file(lr),
        tsv: file("${prefix}.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.{fna,fna.gz}", optional: true),
            files("${prefix}.tsv", optional: true),
            files("supplemental/*"),
            files("${prefix}-*-error.*", optional: true),
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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

    // Unicycler (hybrid: PE + long reads)
    def String is_hybrid = meta.runtype == "hybrid" ? "-l ${lr}" : ""

    // Shovill (short reads)
    def String contig_namefmt = task.ext.contig_namefmt ? task.ext.contig_namefmt : "${prefix}_%05d"
    def Integer shovill_ram   = task.memory.toGiga() - 1
    def String shovill_mode   = meta.single_end == false ? "shovill --R1 ${r1} --R2 ${r2}" : "shovill-se --SE ${se}"

    // Dragonflye (long reads, with optional short read polishing)
    def String dragonflye_fastq = meta.runtype == "short_polish" ? "--reads ${lr} --R1 ${r1} --R2 ${r2}" : "--reads ${lr}"

    // Assembly inputs (for assembly runtype, the original assembly is in lr slot)
    def Boolean use_original_assembly = (meta.runtype.startsWith('assembly') && !task.ext.reassemble) ? true : false
    """
    #==========================================================================================
    # Assemble based on runtype
    #==========================================================================================
    mkdir supplemental
    if [ "${use_original_assembly}" == "true" ]; then
        # Skip assembly and use provided assembly
        echo "Using provided assembly for ${prefix} without re-assembly." > supplemental/assembly-info.txt
        gzip -cd ${assembly} > ${prefix}.fna
    elif [[ "${meta.runtype}" == "hybrid" || "${task.ext.use_unicycler}" == "true" ]]; then
        #======================================================================================
        # Unicycler Assembler
        #
        # Unicycler is used for hybrid assemblies (short + long reads) or when explicitly
        # requested via the --use_unicycler parameter.
        #======================================================================================
        unicycler \\
            -1 ${r1} -2 ${r2} \\
            ${task.ext.args3} \\
            ${is_hybrid} \\
            -o supplemental/ \\
            --threads ${task.cpus}
        sed -r 's/^>([0-9]+)(.*)/>${prefix}_\\1\\2/' supplemental/assembly.fasta > ${prefix}.fna
        mv supplemental/assembly.fasta supplemental/unicycler-unpolished.fasta
        mv supplemental/assembly.gfa supplemental/unicycler-unpolished.gfa
    elif [[ "${meta.runtype}" == "ont" || "${meta.runtype}" == "short_polish" ]]; then
        #======================================================================================
        # Dragonflye Assembler
        #
        # Dragonflye is used for long read assemblies (ONT) with optional short read polishing.
        #======================================================================================
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
        mv supplemental/contigs.fa ${prefix}.fna
    else
        #======================================================================================
        # Shovill Assembler
        #
        # Shovill is used for short read assemblies (Illumina PE or SE).
        #======================================================================================
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
        mv supplemental/contigs.fa ${prefix}.fna

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

    #==========================================================================================
    # Assembly Quality Check
    #==========================================================================================
    TOTAL_CONTIGS=`grep -c "^>" ${prefix}.fna || true`
    if [ "\${TOTAL_CONTIGS}" -gt 0 ]; then
        assembly-scan ${prefix}.fna --prefix ${prefix} > ${prefix}.tsv
        TOTAL_CONTIG_SIZE=\$(cut -f 3 ${prefix}.tsv | tail -n 1)
        if [ "\${TOTAL_CONTIG_SIZE}" -lt ${task.ext.min_genome_size} ]; then
            mv ${prefix}.fna ${prefix}-error.fna
            mv ${prefix}.tsv ${prefix}-error.tsv
            echo "${prefix} assembled size [\${TOTAL_CONTIG_SIZE} bp] is less than the minimum allowed genome
                    size [${task.ext.min_genome_size} bp]. If this is unexpected, please investigate ${prefix} to
                    determine a cause [e.g. metagenomic, contaminants, etc...] for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of ${prefix} will be discontinued." | \
            sed 's/^\\s*//' > ${prefix}-assembly-error.txt
        fi
    else
        mv ${prefix}.fna ${prefix}-error.fna
        echo "${prefix} assembled successfully, but 0 contigs were formed. Please investigate
                ${prefix} to determine a cause [e.g. metagenomic, contaminants, etc...] for this
                outcome. Further assembly-based analysis of ${prefix} will be discontinued." | \
        sed 's/^\\s*//' > ${prefix}-assembly-error.txt
    fi

    #==========================================================================================
    # Cleanup intermediate files
    #==========================================================================================
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
                supplemental/velvet/ \\
                supplemental/assembly.fasta
    fi

    #==========================================================================================
    # Compress and move final outputs
    #==========================================================================================
    if [[ "${task.ext.skip_compression}" == "false" ]]; then
        # Compress based on matched extensions
        pigz -n --best -p ${task.cpus} ${prefix}*.fna
        find supplemental/ -type f | \
            grep -E "\\.fna\$|\\.fasta\$|\\.fa\$|\\.gfa\$" | \
            xargs -I {} pigz -n --best -p ${task.cpus} {}
    fi
    find supplemental/ -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

    # Capture versions (common tools available on all platforms)
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

    # Linux-only tools (long-read assemblers not available on macOS)
    if [[ "\$OSTYPE" != "darwin"* ]]; then
    cat <<-END_LINUX >> versions.yml
        any2fasta: \$(echo \$(any2fasta -v 2>&1) | sed 's/^.*any2fasta //')
        dragonflye: \$(echo \$(dragonflye --version 2>&1) | sed 's/^.*dragonflye //' )
        flye: \$(echo \$(flye --version))
        medaka: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
        minimap2: \$(echo \$(minimap2 --version))
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
        raven: \$(echo \$(raven --version))
    END_LINUX
    fi
    """
}
