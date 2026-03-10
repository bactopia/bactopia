/**
 * Automated quality control, error correction, and read subsampling.
 *
 * A comprehensive QC pipeline that adapts to the input read type:
 * - **Illumina:** Adapter/PhiX removal ([Fastp](https://github.com/OpenGene/fastp) or
 *   [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)), Error Correction
 *   ([Lighter](https://github.com/mourisl/Lighter)), and Subsampling ([Rasusa](https://github.com/mbhall88/rasusa))
 * - **Nanopore:** Adapter removal ([Porechop](https://github.com/rrwick/Porechop)), Quality filtering
 *   ([Nanoq](https://github.com/esteinig/nanoq)), and Subsampling ([Rasusa](https://github.com/mbhall88/rasusa))
 * - **Hybrid:** Processes both short and long reads through their respective pipelines
 * - **Assembly:** Passes through simulated reads from assemblies
 *
 * Generates quality metrics using [fastq-scan](https://github.com/rpetit3/fastq-scan) and optional
 * quality reports using [FastQC](https://github.com/s-andrews/FastQC) (Illumina) and
 * [NanoPlot](https://github.com/wdecoster/NanoPlot) (ONT).
 *
 * @status stable
 * @keywords fastq, qc, adapter removal, error correction, subsampling, fastp, bbduk, lighter, porechop, nanoq, fastqc, nanoplot
 * @tags complexity:complex input-type:multiple output-type:multiple features:conditional-logic,compression,path-workarounds
 * @citation bbtools, fastp, fastqc, fastq_scan, lighter, nanoplot, nanoq, porechop, rasusa
 *
 * @note Uses EMPTY_* placeholder files for optional parameters (adapters, phix)
 *
 * @input record(meta, r1, r2, se, lr, assembly)
 * - `meta`: Groovy Map containing sample information (must include `runtype`, `genome_size`, `species`)
 * - `r1`: Illumina R1 reads (paired-end forward)
 * - `r2`: Illumina R2 reads (paired-end reverse)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT)
 * - `assembly`: Assembly file (FASTA) for assembly-based simulations
 *
 * @input adapters
 * Optional filepath for custom adapter sequences (FASTA)
 *
 * @input phix
 * Optional filepath for custom PhiX sequences (FASTA)
 *
 * @output record(meta, r1, r2, se, lr, assembly, reads_grouped, supplemental, error, results, logs, nf_logs, versions)
 * - `r1`: QC'd Illumina R1 reads (paired-end forward)
 * - `r2`: QC'd Illumina R2 reads (paired-end reverse)
 * - `se`: QC'd single-end Illumina reads
 * - `lr`: QC'd long reads (ONT)
 * - `assembly`: Assembly file (FASTA)
 * - `reads_grouped`: All output FASTQs for publishing
 * - `supplemental`: QC reports (FastQC/NanoPlot), JSON metrics, and error FASTQs if QC failed
 * - `error`: Captured error messages if QC failed (e.g., reads empty after trimming)
 */
nextflow.preview.types = true

process QC {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?, assembly: Path?): Record
    adapters                                                                 : Path?
    phix                                                                     : Path?

    stage:
    stageAs 'input-r1/*', r1
    stageAs 'input-r2/*', r2
    stageAs 'input-se/*', se
    stageAs 'input-lr/*', lr
    stageAs 'input-assembly/*', assembly

    output:
    record(
        meta: meta,
        r1: file("${prefix}_R1.fastq.gz", optional: true),
        r2: file("${prefix}_R2.fastq.gz", optional: true),
        se: file("${prefix}_SE.fastq.gz", optional: true),
        lr: file("${prefix}_ONT.fastq.gz", optional: true),
        assembly: file("assembly/${prefix}.fna.gz", optional: true),
        reads_grouped: files("${prefix}*.fastq.gz", optional: true),
        supplemental: files("supplemental/*", optional: true),
        error: files("*-error.txt", optional: true),
        results: [files("${prefix}*.fastq.gz", optional: true)],
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
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/logs/"
    meta.process_name = task.ext.process_name
    meta.genome_size = _meta.genome_size ?: 0
    meta.species = _meta.species ?: null
    meta.runtype = _meta.runtype

    // Map explicit slots to legacy variable names for minimal shell script changes
    // PE: use r1 (fq1) and r2 (fq2), SE: use se (fq1), ONT: use lr (ont_fq)
    def Path fq1 = r1 ? r1 : se
    def Path fq2 = r2
    def Path ont_fq = lr

    // WF opts
    def Boolean is_assembly = meta.runtype.startsWith('assembly') ? true : false
    def Boolean single_end = r1 && r2 ? false : true
    def String qin = meta.runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    def String adapter_file = adapters.getName() == 'EMPTY_ADAPTERS' ? 'adapters' : adapters.getName()
    def String phix_file = phix.getName() == 'EMPTY_PHIX' ? 'phix' : phix.getName()
    def String adapter_opts = single_end ? "" : "in2=repair-r2.fq out2=adapter-r2.fq"
    def String phix_opts = single_end ? "" : "in2=adapter-r2.fq out2=phix-r2.fq"
    def String lighter_opts = single_end ? "" : "-r phix-r2.fq"
    def String rasusa_opts = single_end ? "-o subsample-r1.fq filt-r1.fq" : "-o subsample-r1.fq -o subsample-r2.fq filt-r1.fq filt-r2.fq"
    def String fastp_fqs = single_end ? "" : "--in2 ${fq2} --out2 filt-r2.fq --detect_adapter_for_pe"

    // For short_polish, inform downstream processes that ONT is the primary (single-end) input
    meta.single_end = meta.runtype == 'short_polish' ? true : single_end

    // set Xmx to 95% of what was allocated, to avoid going over
    xmx = Math.round(task.memory.toBytes() * 0.95)
    """
    #==========================================================================================
    # Helper Functions
    #==========================================================================================
    is_assembly()        { [[ "${is_assembly}" == "true" ]]; }
    is_hybrid()          { [[ "${meta.runtype}" =~ ^(hybrid|short_polish)\$ ]]; }
    is_ont()             { [[ "${meta.runtype}" == "ont" ]]; }
    is_paired()          { [[ "${single_end}" == "false" ]]; }
    is_illumina_se()     { [[ "${single_end}" == "true" ]] && ! is_ont; }
    has_ont_reads()      { is_ont || is_hybrid; }
    has_illumina_reads() { is_paired || is_illumina_se; }
    fastq_scan()         { gzip -cd "\$1" | fastq-scan -g ${meta.genome_size} > "\$2"; }
    compress()           { pigz -p ${task.cpus} -c -n "\$1" > "\$2"; }

    run_nanoplot() {
        local input="\$1"
        local suffix="\$2"
        local outdir="supplemental/${prefix}-\${suffix}"
        
        mkdir -p "\${outdir}"
        NanoPlot ${task.ext.nanoplot_opts} \
            --threads ${task.cpus} \
            --fastq "\${input}" \
            --outdir "\${outdir}/" \
            --prefix "${prefix}-\${suffix}_"
        
        cp "\${outdir}/${prefix}-\${suffix}_NanoPlot-report.html" "supplemental/${prefix}-\${suffix}_NanoPlot-report.html"
        tar -cvf - "\${outdir}/" | pigz --best -p ${task.cpus} > "supplemental/${prefix}-\${suffix}_NanoPlot.tar.gz"
        rm -rf "\${outdir}"
    }

    report_empty_error() {
        echo "After \${1}, ${prefix} FASTQs are empty. Please check the input FASTQs.
            Further analysis is discontinued." | sed 's/^\s*//' >> "${prefix}-\${2}-qc-error.txt"
        ERROR=1
    }

    #==========================================================================================
    # QC Initialization
    #==========================================================================================
    mkdir assembly supplemental
    ERROR=0

    if is_assembly; then
        # Copy the assembly over to assembly
        cp ${assembly} assembly/${prefix}.fna.gz
    else
        touch assembly/EMPTY_ASSEMBLY
    fi

    #==========================================================================================
    # Begin QC
    #==========================================================================================
    if [[ "${task.ext.skip_qc}" == "true" ]]; then
        #======================================================================================
        # User skipped QC, just copy input reads to output
        #======================================================================================
        echo "Sequence QC was skipped for ${prefix}" > ${prefix}-qc-skipped.txt
        if is_paired; then
            cp ${fq1} ${prefix}_R1.fastq.gz
            cp ${fq2} ${prefix}_R2.fastq.gz
        fi

        if has_ont_reads; then
            cp ${ont_fq} ${prefix}_ONT.fastq.gz
        fi

        if is_illumina_se; then
            cp ${fq1} ${prefix}_SE.fastq.gz
        fi
    else
        #======================================================================================
        # Nanopore Read QC
        #======================================================================================
        if has_ont_reads; then
            NANOQ_FQ="${ont_fq}"
            if [[ "${task.ext.use_porechop}" == "true" ]]; then
                # Remove Adapters
                NANOQ_FQ="adapter-ont.fq"
                porechop \
                    --input ${ont_fq} ${task.ext.porechop_opts} \
                    --format fastq \
                    --threads ${task.cpus} > adapter-ont.fq
            fi

            # Quality filter
            nanoq \
                --min-len ${task.ext.ont_minlength} \
                --min-qual ${task.ext.ont_minqual} \
                --input \${NANOQ_FQ} 1> filt-ont.fq
        fi

        #======================================================================================
        # Illumina Read QC
        #
        # Defaults to using Fastp unless BBMap is specified (--use_bbmap)
        #
        # ***** Note: Any errors in this block will set ERROR=1 and skip remaining steps *****
        #======================================================================================
        if has_illumina_reads; then
            if [[ "${task.ext.use_bbmap}" == "true" ]]; then
                # Use BBMap for cleaning reads
                if is_paired; then
                    # Make sure paired-end reads have matching IDs
                    repair.sh \
                        in=${fq1} \
                        in2=${fq2} \
                        out=repair-r1.fq \
                        out2=repair-r2.fq \
                        outs=repair-singles.fq \
                        ain=${task.ext.ain}

                    if [ ! -s repair-r1.fq ]; then
                        report_empty_error "validating read pairs" "paired-match"
                    fi
                else
                    # No need to validate single-end read IDs
                    gunzip -c ${fq1} > repair-r1.fq 
                fi

                if [ "\${ERROR}" -eq "0" ]; then
                    # Remove Adapters
                    bbduk.sh -Xmx${xmx} \
                        in=repair-r1.fq out=adapter-r1.fq ${adapter_opts} \
                        ref=${adapter_file} \
                        k=${task.ext.adapter_k} \
                        ktrim=${task.ext.ktrim} \
                        mink=${task.ext.mink} \
                        hdist=${task.ext.hdist} \
                        tpe=${task.ext.tpe} \
                        tbo=${task.ext.tbo} \
                        threads=${task.cpus} \
                        ftm=${task.ext.ftm} \
                        ${qin} ordered=t ${task.ext.bbduk_opts}

                    if [ ! -s adapter-r1.fq ]; then
                        report_empty_error "adapter removal" "adapter"
                    fi
                fi

                if [ "\${ERROR}" -eq "0" ]; then
                    # Remove PhiX
                    bbduk.sh -Xmx${xmx} \
                        in=adapter-r1.fq out=phix-r1.fq ${phix_opts} \
                        ref=${phix_file} \
                        k=${task.ext.phix_k} \
                        hdist=${task.ext.hdist} \
                        tpe=${task.ext.tpe} \
                        tbo=${task.ext.tbo} \
                        qtrim=${task.ext.qtrim} \
                        trimq=${task.ext.trimq} \
                        minlength=${task.ext.minlength} \
                        minavgquality=${task.ext.maq} \
                        ${qin} qout=${task.ext.qout} \
                        tossjunk=${task.ext.tossjunk} \
                        threads=${task.cpus} \
                        ordered=t ${task.ext.bbduk_opts}

                    if [ ! -s phix-r1.fq ]; then
                        report_empty_error "PhiX removal" "phix"
                    fi
                fi
            else
                # QC with fastp
                fastp \
                    --in1 ${fq1} --out1 filt-r1.fq ${fastp_fqs} \
                    --thread ${task.cpus} \
                    --json supplemental/${prefix}.fastp.json \
                    --html supplemental/${prefix}.fastp.html ${task.ext.fastp_opts} 2> ${prefix}-fastp.log
            fi
        fi

        # Error Correction
        if [ "\${ERROR}" -eq "0" ]; then
            if [[ "${task.ext.use_bbmap}" == "true" ]]; then
                if [ "${task.ext.skip_error_correction}" == "false" ] && [ "${meta.genome_size}" -gt "0" ]; then
                    lighter -od . -r phix-r1.fq ${lighter_opts} -K 31 ${meta.genome_size} -maxcor 1 -zlib 0 -t ${task.cpus}
                    mv phix-r1.cor.fq filt-r1.fq
                    if is_paired; then
                        mv phix-r2.cor.fq filt-r2.fq
                    fi
                else
                    echo "Skipping error correction"
                    ln -s phix-r1.fq filt-r1.fq
                    if is_paired; then
                        ln -s phix-r2.fq filt-r2.fq
                    fi
                fi
            elif [[ "${meta.runtype}" == "ont" ]]; then
                echo "Skipping error correction. Have a recommended ONT error corrector? Let me know!"
            else
                echo "Skipping error correction since fastp was used"
            fi
        fi

        # Reduce Coverage
        if [ "\${ERROR}" -eq "0" ]; then
            TOTAL_BP=\$(( ${task.ext.coverage}*${meta.genome_size} ))
            if (( \${TOTAL_BP} > 0 )); then
                if has_ont_reads; then
                    rasusa reads \
                        -c ${task.ext.coverage} \
                        -g ${meta.genome_size} \
                        -s ${task.ext.sampleseed} \
                        filt-ont.fq 1> subsample-ont.fq
                fi

                if has_illumina_reads; then
                    if [ -f filt-r1.fq ]; then
                        rasusa reads \
                            -c ${task.ext.coverage} \
                            -g ${meta.genome_size} \
                            -s ${task.ext.sampleseed} \
                            ${rasusa_opts}
                    fi
                fi
            else
                echo "Skipping coverage reduction"
                if has_ont_reads; then
                    ln -s filt-ont.fq subsample-ont.fq
                fi

                if has_illumina_reads; then
                    ln -s filt-r1.fq subsample-r1.fq
                    if is_paired; then
                        ln -s filt-r2.fq subsample-r2.fq
                    fi
                fi
            fi
        fi

        # Compress
        if [ "\${ERROR}" -eq "0" ]; then
            if is_paired; then
                compress subsample-r1.fq ${prefix}_R1.fastq.gz
                compress subsample-r2.fq ${prefix}_R2.fastq.gz
            fi 

            if has_ont_reads; then
                compress subsample-ont.fq ${prefix}_ONT.fastq.gz
            fi

            if is_illumina_se; then
                compress subsample-r1.fq ${prefix}_SE.fastq.gz
            fi

            if [ "${task.ext.keep_all_files}" == "false" ]; then
                # Remove remaining intermediate FASTQ files
                rm *.fq
            fi
        fi
    fi

    #==========================================================================================
    # Calculate quality stats before and after QC, optionally include plots (Default: true)
    #==========================================================================================
    if [ "\${ERROR}" -eq "0" ]; then
        # fastq-scan
        if is_paired; then
            fastq_scan ${fq1} supplemental/${prefix}_R1-original.json
            fastq_scan ${fq2} supplemental/${prefix}_R2-original.json
            fastq_scan ${prefix}_R1.fastq.gz supplemental/${prefix}_R1-final.json
            fastq_scan ${prefix}_R2.fastq.gz supplemental/${prefix}_R2-final.json
        fi

        if has_ont_reads; then
            fastq_scan ${ont_fq} supplemental/${prefix}_ONT-original.json
            fastq_scan ${prefix}_ONT.fastq.gz supplemental/${prefix}_ONT-final.json
        fi

        if is_illumina_se; then
            fastq_scan ${fq1} supplemental/${prefix}_SE-original.json
            fastq_scan ${prefix}_SE.fastq.gz supplemental/${prefix}_SE-final.json
        fi

        # FastQC and NanoPlot
        if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
            if has_ont_reads; then
                run_nanoplot "${ont_fq}" "original"
                run_nanoplot "${prefix}_ONT.fastq.gz" "final"
            fi

            if has_illumina_reads; then
                # Fix issue on HPC when /tmp not writable
                mkdir -p ./tmp
                fastqs=()
                if is_paired; then
                    ln -s ${fq1} ${prefix}_R1-original.fastq.gz
                    ln -s ${fq2} ${prefix}_R2-original.fastq.gz
                    ln -s ${prefix}_R1.fastq.gz ${prefix}_R1-final.fastq.gz
                    ln -s ${prefix}_R2.fastq.gz ${prefix}_R2-final.fastq.gz
                    fastqs=(${prefix}_R1-original.fastq.gz ${prefix}_R2-original.fastq.gz ${prefix}_R1-final.fastq.gz ${prefix}_R2-final.fastq.gz)
                else
                    ln -s ${fq1} ${prefix}_SE-original.fastq.gz
                    ln -s ${prefix}_SE.fastq.gz ${prefix}_SE-final.fastq.gz
                    fastqs=(${prefix}_SE-original.fastq.gz ${prefix}_SE-final.fastq.gz)
                fi

                # Run FastQC
                fastqc \
                    --noextract \
                    --dir ./tmp \
                    -f fastq \
                    -t ${task.cpus} "\${fastqs[@]}"
                mv *_fastqc.html *_fastqc.zip supplemental/
                rm -rf tmp/ *-original.fastq.gz *-final.fastq.gz
            fi
        fi
    fi

    #==========================================================================================
    # Final QC check
    #
    # This step ensures that the final reads meet minimum requirements. Users can skip this
    # step with --skip_fastq_check
    #==========================================================================================
    if [ "${task.ext.skip_fastq_check}" == "false" ]; then
        # Only check for errors if we haven't already found them
        if [ "\${ERROR}" -eq "0" ]; then
            # Check minimum coverage is met
            if [ "${meta.runtype}" == "hybrid" ]; then
                # Hybrid assembly, base checks on Illumina Reads
                gzip -cd ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            elif [ "${meta.runtype}" == "short_polish" ]; then
                # Short read polishing, base checks on Nanopore Reads
                fastq_scan ${prefix}_ONT.fastq.gz temp.json
            else
                # base checks on which ever reads are available
                gzip -cd *.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            fi
            FINAL_BP=\$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\\1/')
            MIN_COVERAGE=\$(( ${task.ext.min_coverage}*${meta.genome_size} ))
            if [ \${FINAL_BP} -lt \${MIN_COVERAGE} ]; then
                echo "After QC, ${prefix} FASTQ(s) contain \${FINAL_BP} total basepairs. This does
                        not exceed the required minimum \${MIN_COVERAGE} bp [${task.ext.min_coverage}x coverage]. Further analysis
                        is discontinued." | \\
                sed 's/^\\s*//' > ${prefix}-low-sequence-coverage-error.txt
                ERROR=2
            fi
            rm temp.json

            # Final QC checks using check-fastqs.py
            OPTS="--sample ${prefix} --min_basepairs ${task.ext.min_basepairs} --min_reads ${task.ext.min_reads} --min_proportion ${task.ext.min_proportion} --runtype ${meta.runtype}"
            if is_paired; then
                # Check paired-end reads have same read counts
                if ! reformat.sh in1=${prefix}_R1.fastq.gz in2=${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                    echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                        Further analysis is discontinued." | \\
                    sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
                    ERROR=2
                else
                    rm -f ${prefix}-paired-end-error.txt
                fi
                check-fastqs.py --fq1 supplemental/${prefix}_R1-final.json  --fq2 supplemental/${prefix}_R2-final.json \${OPTS}
            fi

            if has_ont_reads; then
                check-fastqs.py --fq1 supplemental/${prefix}_ONT-final.json \${OPTS}
            fi

            if is_illumina_se; then
                check-fastqs.py --fq1 supplemental/${prefix}_SE-final.json \${OPTS}
            fi

            if ls *-error.txt 1>/dev/null 2>&1; then
                ERROR=2
            fi

        fi
    fi

    if is_assembly; then
        touch supplemental/reads-simulated-from-assembly.txt
    fi

    #==========================================================================================
    # Handle Errors
    # 
    # ERROR=0 : QC process successful
    # ERROR=1 : Original reads could not be processed (e.g. empty after trimming)
    # ERROR=2 : Final QC checks failed (e.g. insufficient reads after QC)
    #==========================================================================================
    if [ "\${ERROR}" -eq "1" ]; then
        # Use original reads since we couldn't process them
        if is_paired; then
            cp ${fq1} supplemental/${prefix}_R1.error-fastq.gz
            cp ${fq2} supplemental/${prefix}_R2.error-fastq.gz
            
            if [ -s repair-singles.fq ]; then
                compress repair-singles.fq supplemental/${prefix}.singles.fastq.gz
            fi
        fi

        if has_ont_reads; then
            cp ${ont_fq} supplemental/${prefix}_ONT.error-fastq.gz
        fi

        if is_illumina_se; then
            cp ${fq1} supplemental/${prefix}_SE.error-fastq.gz
        fi
    elif [ "\${ERROR}" -eq "2" ]; then
        # Use the final processed reads which failed QC checks
        if is_paired; then
            mv ${prefix}_R1.fastq.gz supplemental/${prefix}_R1.error-fastq.gz
            mv ${prefix}_R2.fastq.gz supplemental/${prefix}_R2.error-fastq.gz

            if [ -s repair-singles.fq ]; then
                compress repair-singles.fq supplemental/${prefix}.singles-error.fastq.gz
            fi
        fi

        if has_ont_reads; then
            mv ${prefix}_ONT.fastq.gz supplemental/${prefix}_ONT.error-fastq.gz
        fi

        if is_illumina_se; then
            mv ${prefix}_SE.fastq.gz supplemental/${prefix}_SE.error-fastq.gz
        fi
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBTools version //;s/ .*\$//')
        fastp: \$(echo \$(fastp --version 2>&1) | sed -e "s/fastp //g")
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        lighter: \$(echo \$(lighter -v 2>&1) | sed 's/Lighter v//')
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        porechop: \$(echo \$(porechop --version 2>&1))
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
    END_VERSIONS

    if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
    cat <<-END_EXTRA >> versions.yml
    "${task.process}":
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
        nanoplot: \$(echo \$(NanoPlot -v 2>&1) | sed 's/.*NanoPlot //')
    END_EXTRA
    fi
    """
}
