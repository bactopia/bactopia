process QC {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fq, stageAs: 'inputs/*'), path(extra)
    path adapters
    path phix

    output:
    tuple val(meta), path("${prefix}*.fastq.gz"), path("extra/*"), emit: fastq, optional: true
    tuple val(meta), path("${prefix}*.fastq.gz")                 , emit: fastq_only, optional: true
    tuple val(meta), path("${prefix}*-fastq.gz")                 , emit: error_fastq, optional: true
    tuple val(meta), path("${prefix}.txt") , emit: txt, optional: true
    tuple val(meta), path("supplemental/*"), emit: supplemental, optional: true
    tuple val(meta), path("*-error.txt")   , emit: error, optional: true
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
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.single_end = fq[1] == null ? true : false
    meta.genome_size = _meta.genome_size ?: 0
    meta.runtype = _meta.runtype

    // WF specific parameters
    meta.single_end = fq[1] == null ? true : false
    is_assembly = meta.runtype.startsWith('assembly') ? true : false
    qin = meta.runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapter_file = adapters.getName() == 'EMPTY_ADAPTERS' ? 'adapters' : adapters
    phix_file = phix.getName() == 'EMPTY_PHIX' ? 'phix' : phix
    adapter_opts = meta.single_end ? "" : "in2=repair-r2.fq out2=adapter-r2.fq"
    phix_opts = meta.single_end ? "" : "in2=adapter-r2.fq out2=phix-r2.fq"
    lighter_opts = meta.single_end ? "" : "-r phix-r2.fq"
    reformat_opts = meta.single_end ? "" : "in2=filt-r2.fq out2=subsample-r2.fq"
    fastp_fqs = meta.single_end ? "" : "--in2 ${fq[1]} --out2 filt-r2.fq --detect_adapter_for_pe"
    ont_fq = meta.runtype == 'ont' ? fq[0] : extra[0]
    meta.single_end = meta.runtype == 'short_polish' ? true : meta.single_end

    // set Xmx to 95% of what was allocated, to avoid going over
    xmx = Math.round(task.memory.toBytes()*0.95)
    """
    mkdir -p supplemental
    touch supplemental/.${meta.runtype}
    ERROR=0
    MIN_COVERAGE=\$(( ${task.ext.min_coverage}*${meta.genome_size} ))
    TOTAL_BP=\$(( ${task.ext.coverage}*${meta.genome_size} ))
    ENABLE_ONT=0
    ENABLE_ILLUMINA=0
    IS_HYBRID=0

    if [ "${meta.runtype}" == "hybrid" ] || [ "${meta.runtype}" == "short_polish" ]; then
        ENABLE_ONT=1
        ENABLE_ILLUMINA=1
        IS_HYBRID=1
    elif [ "${meta.runtype}" == "ont" ]; then
        ENABLE_ONT=1
    else
        ENABLE_ILLUMINA=1
    fi

    mkdir extra/
    if [[ "${is_assembly}" == "true" ]]; then
        # Copy the assembly over to extra
        cp ${extra[0]} extra/
    else
        touch extra/EMPTY_EXTRA
    fi

    if [[ "${task.ext.skip_qc}" == "true" ]]; then
        echo "Sequence QC was skipped for ${prefix}" > ${prefix}-qc-skipped.txt
        if [ "\${IS_HYBRID}" -eq "1" ]; then
            # Illumina Paired-End Reads and Nanopore Reads
            cp ${fq[0]} ${prefix}_R1.fastq.gz
            cp ${fq[1]} ${prefix}_R2.fastq.gz
            cp ${extra[0]} ${prefix}.fastq.gz
        elif [ "${meta.single_end}" == "false" ]; then
            # Paired-End Reads
            cp ${fq[0]} ${prefix}_R1.fastq.gz
            cp ${fq[1]} ${prefix}_R2.fastq.gz
        else
            # Single-End Reads
            cp ${fq[0]} ${prefix}.fastq.gz
        fi
    else
        if [ "\${ENABLE_ONT}" -eq "1" ]; then
            # QC the Nanopore reads
            if [[ "${task.ext.use_porechop}" == "true" ]]; then
                # Remove Adapters
                porechop --input ${ont_fq} ${task.ext.porechop_opts} \
                    --format fastq \
                    --threads ${task.cpus} > adapter-ont.fq

                # Quality filter
                nanoq --min-len ${task.ext.ont_minlength} \
                    --min-qual ${task.ext.ont_minqual} \
                    --input adapter-r1.fq 1> filt-ont.fq
            else 
                # Quality filter
                nanoq --min-len ${task.ext.ont_minlength} \
                    --min-qual ${task.ext.ont_minqual} \
                    --input  ${ont_fq} 1> filt-ont.fq
            fi
        fi

        if [ "\${ENABLE_ILLUMINA}" -eq "1" ]; then
            if [[ "${task.ext.use_bbmap}" == "true" ]]; then
                # Use BBMap for cleaning reads
                # Illumina Reads
                # Validate paired-end reads if necessary
                if [[ "${meta.single_end}" == "false" || "\${IS_HYBRID}" -eq "1" ]]; then
                    # Make sure paired-end reads have matching IDs
                    repair.sh \
                        in=${fq[0]} \
                        in2=${fq[1]} \
                        out=repair-r1.fq \
                        out2=repair-r2.fq \
                        outs=repair-singles.fq \
                        ain=${task.ext.ain}

                    if [ ! -s repair-r1.fq ]; then
                        ERROR=1
                        echo "After validating read pairs, ${prefix} FASTQs are empty. Please check the input FASTQs.
                            Further analysis is discontinued." | \
                        sed 's/^\\s*//' >> ${prefix}-paired-match-error.txt
                    fi
                else
                    gunzip -c ${fq[0]} > repair-r1.fq 
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
                        ERROR=1
                        echo "After adapter removal, ${prefix} FASTQs are empty. Please check the input FASTQs.
                            Further analysis is discontinued." | \
                        sed 's/^\\s*//' >> ${prefix}-adapter-qc-error.txt
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
                        trimq=${task.trimq} \
                        minlength=${task.ext.minlength} \
                        minavgquality=${task.ext.maq} \
                        ${qin} qout=${task.ext.qout} \
                        tossjunk=${task.ext.tossjunk} \
                        threads=${task.cpus} \
                        ordered=t ${task.ext.bbduk_opts}

                    if [ ! -s phix-r1.fq ]; then
                        ERROR=1
                        echo "After PhiX removal, ${prefix} FASTQs are empty. Please check the input FASTQs.
                            Further analysis is discontinued." | \
                        sed 's/^\\s*//' >> ${prefix}-phix-qc-error.txt
                    fi
                fi
            else
                # QC with fastp
                mkdir -p supplemental/summary/
                fastp \
                    --in1 ${fq[0]} --out1 filt-r1.fq ${fastp_fqs} \
                    --thread ${task.cpus} \
                    --json supplemental/summary/${prefix}.fastp.json \
                    --html supplemental/summary/${prefix}.fastp.html ${task.ext.fastp_opts} 2> ${prefix}-fastp.log
            fi
        fi

        # Error Correction
        if [ "\${ERROR}" -eq "0" ]; then
            if [[ "${task.ext.use_bbmap}" == "true" ]]; then
                if [ "${task.ext.skip_error_correction}" == "false" ] && [ "${meta.genome_size}" -gt "0" ]; then
                    lighter -od . -r phix-r1.fq ${lighter_opts} -K 31 ${meta.genome_size} -maxcor 1 -zlib 0 -t ${task.cpus}
                    mv phix-r1.cor.fq filt-r1.fq
                    if [[ "${meta.single_end}" == "false" || "\${IS_HYBRID}" -eq "1" ]]; then
                        mv phix-r2.cor.fq filt-r2.fq
                    fi
                else
                    echo "Skipping error correction"
                    ln -s phix-r1.fq filt-r1.fq
                    if [[ "${meta.single_end}" == "false" || "\${IS_HYBRID}" -eq "1" ]]; then
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
            if (( \${TOTAL_BP} > 0 )); then
                if [ "\${ENABLE_ONT}" -eq "1" ]; then
                    rasusa reads \
                        -c ${task.ext.coverage} \
                        -g ${meta.genome_size} \
                        -s ${task.ext.sampleseed} \
                        filt-ont.fq 1> subsample-ont.fq
                fi

                if [ "\${ENABLE_ILLUMINA}" -eq "1" ]; then
                    if [ -f filt-r1.fq ]; then
                        reformat.sh -Xmx${xmx} \
                            in=filt-r1.fq out=subsample-r1.fq ${reformat_opts} \
                            samplebasestarget=\${TOTAL_BP} \
                            sampleseed=${task.ext.sampleseed} \
                            overwrite=t
                    fi
                fi
            else
                echo "Skipping coverage reduction"
                if [ "\${ENABLE_ONT}" -eq "1" ]; then
                    ln -s filt-ont.fq subsample-ont.fq
                fi

                if [ "\${ENABLE_ILLUMINA}" -eq "1" ]; then
                    ln -s filt-r1.fq subsample-r1.fq
                    if [[ "${meta.single_end}" == "false" || "\${IS_HYBRID}" -eq "1" ]]; then
                        ln -s filt-r2.fq subsample-r2.fq
                    fi
                fi
            fi
        fi

        # Compress
        if [ "\${ERROR}" -eq "0" ]; then
            if [ "\${IS_HYBRID}" -eq "1" ]; then
                # Illumina Paired-End Reads and Nanopore Reads
                pigz -p ${task.cpus} -c -n subsample-r1.fq > ${prefix}_R1.fastq.gz
                pigz -p ${task.cpus} -c -n subsample-r2.fq > ${prefix}_R2.fastq.gz
                pigz -p ${task.cpus} -c -n subsample-ont.fq > ${prefix}.fastq.gz
            elif [ "\${ENABLE_ONT}" -eq "1" ]; then
                # Nanopore Reads
                pigz -p ${task.cpus} -c -n subsample-ont.fq > ${prefix}.fastq.gz
            else
                # Illumina Reads
                if [ "${meta.single_end}" == "false" ]; then
                    pigz -p ${task.cpus} -c -n subsample-r1.fq > ${prefix}_R1.fastq.gz
                    pigz -p ${task.cpus} -c -n subsample-r2.fq > ${prefix}_R2.fastq.gz
                else
                    pigz -p ${task.cpus} -c -n subsample-r1.fq > ${prefix}.fastq.gz
                fi
            fi

            if [ "${task.ext.keep_all_files}" == "false" ]; then
                # Remove remaining intermediate FASTQ files
                rm *.fq
            fi
        fi
    fi

    # Quality stats before and after QC
    if [ "\${ERROR}" -eq "0" ]; then
        mkdir -p supplemental/summary/
        # fastq-scan
        if [ "\${IS_HYBRID}" -eq "1" ]; then
            # Illumina Paired-End Reads
            gzip -cd ${fq[0]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R1-original.json
            gzip -cd ${fq[1]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R2-original.json
            gzip -cd ${prefix}_R1.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R1-final.json
            gzip -cd ${prefix}_R2.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R2-final.json

            # Nanopore Reads
            gzip -cd ${extra[0]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}-original.json
            gzip -cd ${prefix}.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}-final.json
        elif [[ "${meta.single_end}" == "false" ]]; then
            # Paired-End Reads
            gzip -cd ${fq[0]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R1-original.json
            gzip -cd ${fq[1]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R2-original.json
            gzip -cd ${prefix}_R1.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R1-final.json
            gzip -cd ${prefix}_R2.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}_R2-final.json
        else
            # Single-End Reads
            gzip -cd ${fq[0]} | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}-original.json
            gzip -cd ${prefix}.fastq.gz | fastq-scan -g ${meta.genome_size} > supplemental/summary/${prefix}-final.json
        fi

        # FastQC and NanoPlot
        if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
            if [ "\${ENABLE_ONT}" -eq "1" ]; then
                # Nanopore Plots
                mkdir supplemental/summary/${prefix}-original supplemental/summary/${prefix}-final
                NanoPlot ${task.ext.nanoplot_opts} \
                    --threads ${task.cpus} \
                    --fastq ${ont_fq} \
                    --outdir supplemental/summary/${prefix}-original/ \
                    --prefix ${prefix}-original_
                cp supplemental/summary/${prefix}-original/${prefix}-original_NanoPlot-report.html supplemental/summary/${prefix}-original_NanoPlot-report.html
                tar -cvf - supplemental/summary/${prefix}-original/ | pigz --best -p ${task.cpus} > supplemental/summary/${prefix}-original_NanoPlot.tar.gz

                NanoPlot ${task.ext.nanoplot_opts} \
                    --threads ${task.cpus} \
                    --fastq ${prefix}.fastq.gz \
                    --outdir supplemental/summary/${prefix}-final/ \
                    --prefix ${prefix}-final_
                cp supplemental/summary/${prefix}-final/${prefix}-final_NanoPlot-report.html supplemental/summary/${prefix}-final_NanoPlot-report.html
                tar -cvf - supplemental/summary/${prefix}-final/ | pigz --best -p ${task.cpus} > supplemental/summary/${prefix}-final_NanoPlot.tar.gz
                rm -rf supplemental/summary/${prefix}-original/ supplemental/summary/${prefix}-final/
            fi

            if [ "\${ENABLE_ILLUMINA}" -eq "1" ]; then
                # Fix issue on HPC when /tmp not writable
                mkdir -p ./tmp
                if [[ "${meta.single_end}" == "false" || "\${IS_HYBRID}" -eq "1" ]]; then
                    # Paired-End Reads
                    ln -s ${fq[0]} ${prefix}_R1-original.fastq.gz
                    ln -s ${fq[1]} ${prefix}_R2-original.fastq.gz
                    ln -s ${prefix}_R1.fastq.gz ${prefix}_R1-final.fastq.gz
                    ln -s ${prefix}_R2.fastq.gz ${prefix}_R2-final.fastq.gz
                    fastqc \
                        --noextract \
                        --dir ./tmp \
                        -f fastq \
                        -t ${task.cpus} \
                        ${prefix}_R1-original.fastq.gz ${prefix}_R2-original.fastq.gz ${prefix}_R1-final.fastq.gz ${prefix}_R2-final.fastq.gz
                else
                    # Single-End Reads
                    ln -s ${fq[0]} ${prefix}-original.fastq.gz
                    ln -s ${prefix}.fastq.gz ${prefix}-final.fastq.gz
                    fastqc \
                        --noextract \
                        --dir ./tmp \
                        -f fastq \
                        -t ${task.cpus} \
                        ${prefix}-original.fastq.gz ${prefix}-final.fastq.gz
                fi
                mv *_fastqc.html *_fastqc.zip supplemental/summary/
                rm -rf tmp/ *-original.fastq.gz *-final.fastq.gz
            fi
        fi
    fi

    # Final QC check
    if [ "${task.ext.skip_fastq_check}" == "false" ]; then
        # Only check for errors if we haven't already found them
        if [ "\${ERROR}" -eq "0" ]; then
            if [ "${meta.run_type}" == "hybrid" ]; then
                # Hybrid assembly, base checks on Illumina Reads
                gzip -cd ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            elif [ "${meta.run_type}" == "short_polish" ]; then
                # Short read polishing, base checks on Nanopore Reads
                gzip -cd ${prefix}.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            else
                # base checks on which ever reads are available
                gzip -cd *.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            fi
            FINAL_BP=\$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\\1/')
            rm temp.json

            if [ \${FINAL_BP} -lt \${MIN_COVERAGE} ]; then
                ERROR=1
                echo "After QC, ${prefix} FASTQ(s) contain \${FINAL_BP} total basepairs. This does
                        not exceed the required minimum \${MIN_COVERAGE} bp [${task.ext.min_coverage}x coverage]. Further analysis 
                        is discontinued." | \
                sed 's/^\\s*//' > ${prefix}-low-sequence-coverage-error.txt
                ERROR=2
            fi

            # Check paired-end reads have same read counts
            OPTS="--sample ${prefix} --min_basepairs ${task.ext.min_basepairs} --min_reads ${task.ext.min_reads} --min_proportion ${task.ext.min_proportion} --runtype ${meta.runtype}"
            if [ "\${IS_HYBRID}" -eq "1" ]; then
                # Illumina Paired-End Reads and Nanopore Reads for hybrid assembly
                gzip -cd ${prefix}_R1.fastq.gz | fastq-scan > r1.json
                gzip -cd ${prefix}_R2.fastq.gz | fastq-scan > r2.json
                gzip -cd ${prefix}.fastq.gz | fastq-scan > ont.json
                if ! reformat.sh in1=${prefix}_R1.fastq.gz in2=${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                    ERROR=2
                    echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
                else
                    rm -f ${prefix}-paired-end-error.txt
                fi

                if [ "${meta.run_type}" == "hybrid" ]; then
                    # Base check on Illumina Reads
                    if ! check-fastqs.py --fq1 r1.json --fq2 r2.json \${OPTS}; then
                        ERROR=2
                    fi
                else
                    # Base check on Oxford Nanopore Reads
                    if ! check-fastqs.py --fq1 r1.json \${OPTS}; then
                        ERROR=2
                    fi
                fi

                rm -f r1.json r2.json ont.json
            elif [ -f  "${prefix}_R2.fastq.gz" ]; then
                # Paired-end
                gzip -cd ${prefix}_R1.fastq.gz | fastq-scan > r1.json
                gzip -cd ${prefix}_R2.fastq.gz | fastq-scan > r2.json
                if ! reformat.sh in1=${prefix}_R1.fastq.gz in2=${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                    ERROR=2
                    echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
                else
                    rm -f ${prefix}-paired-end-error.txt
                fi
                if ! check-fastqs.py --fq1 r1.json --fq2 r2.json \${OPTS}; then
                    ERROR=2
                fi
                rm r1.json r2.json
            else
                # Single-end
                gzip -cd ${prefix}.fastq.gz | fastq-scan > r1.json
                if ! check-fastqs.py --fq1 r1.json \${OPTS}; then
                    ERROR=2
                fi
                rm r1.json
            fi
        fi
    fi

    if [ "${is_assembly}" == "true" ]; then
        touch supplemental/reads-simulated-from-assembly.txt
    fi

    if [ "\${ERROR}" -eq "1" ]; then
        if [ "\${IS_HYBRID}" -eq "1" ]; then
            cp ${fq[0]} supplemental/${prefix}_R1.error-fastq.gz
            cp ${fq[1]} supplemental/${prefix}_R2.error-fastq.gz
            cp ${extra[0]} supplemental/${prefix}.error-fastq.gz
            if [ ! -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > supplemental/${prefix}.error-fastq.gz
            fi
        elif [ "${meta.single_end}" == "false" ]; then
            cp ${fq[0]} supplemental/${prefix}_R1.error-fastq.gz
            cp ${fq[1]} supplemental/${prefix}_R2.error-fastq.gz
            if [ ! -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > supplemental/${prefix}.singles-fastq.gz
            fi
        else
            cp ${fq[0]} supplemental/${prefix}.error-fastq.gz
        fi
    elif [ "\${ERROR}" -eq "2" ]; then
        if [ -f ${prefix}_R1.fastq.gz ]; then
            mv ${prefix}_R1.fastq.gz ${prefix}_R1.error-fastq.gz
            mv ${prefix}_R2.fastq.gz ${prefix}_R2.error-fastq.gz

            if [ -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > ${prefix}.singles-fastq.gz
            fi
        fi

        if [ -f ${prefix}.fastq.gz ]; then
            mv ${prefix}.fastq.gz ${prefix}.error-fastq.gz
        fi
    fi

    # Capture versions
    if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBTools version //;s/ .*\$//')
        fastp: \$(echo \$(fastp --version 2>&1) | sed -e "s/fastp //g")
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        lighter: \$(echo \$(lighter -v 2>&1) | sed 's/Lighter v//')
        nanoplot: \$(echo \$(NanoPlot -v 2>&1) | sed 's/.*NanoPlot //')
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        porechop: \$(echo \$(porechop --version 2>&1))
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
    END_VERSIONS
    else
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
    fi
    """
}
