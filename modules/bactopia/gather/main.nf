/**
 * Search, validate, gather, or simulate input samples.
 *
 * This process is the entry point for data ingestion. It handles:
 * - **Validation:** Verifies FASTQ formatting and gzip integrity.
 * - **Merging:** Combines multiple runs (lanes) into a single sample.
 * - **Downloading:** Fetches reads (SRA/ENA) or assemblies (NCBI) from accessions.
 * - **Simulation:** Generates synthetic reads from assemblies using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) to enable read-based analysis.
 *
 * Uses explicit named slots for input and output reads:
 * - Input accepts Set<Path> for each slot (pre-merge, supports multiple files)
 * - Output emits Path? for each slot (post-merge, single consolidated file or null)
 *
 * @status stable
 * @keywords fastq, validation, sra, ena, download, merging, simulation, art, ncbi
 * @tags complexity:complex input-type:multiple output-type:multiple features:internet-access,resource-download,conditional-logic
 * @citation bactopia, art, fastq_dl, fastq_scan, ncbigenomedownload, pigz
 *
 * @input record(meta, r1_files, r2_files, se_files, lr_files, fna_files)
 * - `meta`: Groovy Map containing sample information
 * - `r1_files`: Illumina R1 read files (Set, elements may be null)
 * - `r2_files`: Illumina R2 read files (Set, elements may be null)
 * - `se_files`: Single-end read files (Set, elements may be null)
 * - `lr_files`: Long read files (ONT) or assembly for simulation (Set, elements may be null)
 * - `fna_files`: Input or downloaded assembly file (Set, elements may be null)
 *
 * @output record(meta, r1?, r2?, se?, lr?, fna?, tsv, results, logs, nf_logs, versions)
 * - `r1?`: Merged Illumina R1 read file
 * - `r2?`: Merged Illumina R2 read file
 * - `se?`: Merged single-end read file
 * - `lr?`: Merged long read file (ONT)
 * - `fna?`: Assembly file
 * - `tsv`: A tab-delimited metadata file describing the valid samples
 *
 * @results additional
 * - `*-error.txt`: Error files from failed downloads, validation, or accession checks
 * - `*-merged.txt`: Log of merged read sets when multiple runs are combined
 */
nextflow.preview.types = true

// bactopia-lint: ignore M026
process GATHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        r1_files: Set<Path?>,
        r2_files: Set<Path?>,
        se_files: Set<Path?>,
        lr_files: Set<Path?>,
        fna_files: Set<Path?>
    )

    stage:
    stageAs r1_files, "staging/r1/*???-r1"
    stageAs r2_files, "staging/r2/*???-r2"
    stageAs se_files, "staging/se/*???-se"
    stageAs lr_files, "staging/lr/*???-lr"
    stageAs fna_files, "staging/fna/*???-assembly"

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        r1: file("fastqs/${prefix}_R1.fastq.gz", optional: true),
        r2: file("fastqs/${prefix}_R2.fastq.gz", optional: true),
        se: file("fastqs/${prefix}_SE.fastq.gz", optional: true),
        lr: file("fastqs/${prefix}_ONT.fastq.gz", optional: true),
        fna: file("assembly/${prefix}.fna.gz", optional: true),
        tsv: file("${prefix}-meta.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}-meta.tsv"),
            files("*-{error,merged}.txt", optional: true),
            files("failed-tests-fastqs/*", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml", optional: true)
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Handler variables
    def String wfPath = task.ext.wf == "teton" ? "teton/main" : "main"
    def String runtype = _meta.runtype
    def Map runtypes = [
        'hybrid-merge-pe': 'hybrid',
        'short_polish-merge-pe': 'short_polish',
        'merge-pe': 'paired-end',
        'merge-se': 'single-end',
        'sra_accession_ont': 'ont'
    ]

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/${wfPath}/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/${wfPath}/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.original_runtype = _meta.runtype
    meta.genome_size = _meta.genome_size
    meta.species = _meta.species
    meta.is_compressed = task.ext.skip_compression ? false : true
    meta.runtype = runtypes[runtype] ?: runtype

    // WF specific parameters
    def String no_cache = task.ext.no_cache ? '-N' : ''
    def String archive = task.ext.use_ena ? (task.attempt >= 4 ? "SRA" : "ENA") : "SRA"
    def String section = runtype == 'assembly_accession' ? (prefix.startsWith('GCF') ? 'refseq' : 'genbank') : ''
    def Integer fcov = task.ext.coverage.toInteger() == 0 ? 150 : Math.round(task.ext.coverage.toInteger() * 1.5)

    // Determine what reads we have based on the explicit slots
    // For paired-end, get the first file from each set (single file case) or handle merge
    def Path r1_first = r1_files.size() > 0 ? r1_files.toList()[0] : null
    def Path r2_first = r2_files.size() > 0 ? r2_files.toList()[0] : null
    def Path se_first = se_files.size() > 0 ? se_files.toList()[0] : null
    def Path lr_first = lr_files.size() > 0 ? lr_files.toList()[0] : null
    def Path fna_first = fna_files.size() > 0 ? fna_files.toList()[0] : null
    def String qin = runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    """
    #==========================================================================================
    # Helper functions
    #===========================================================================================
    check_gzip_or_empty() {
        local fq="\$1"
        if ! gzip -t "\${fq}"; then
            # Failed Gzip test
            echo "${prefix} FASTQs failed Gzip tests. Please check the input FASTQs. Further analysis is discontinued." > "${prefix}-gzip-error.txt"
            return 1
        fi
        # Empty FASTQ file
        echo "${prefix} FASTQs are empty. Please check the input FASTQs. Further analysis is discontinued." > "${prefix}-empty-error.txt"
        return 1
    }

    MERGED="multiple-read-sets-merged.txt"
    merge_runs() {
        local label="\$1"
        local pattern="\$2"
        local output="\$3"

        echo "\${label}:" >> \${MERGED}
        find -name "\${pattern}" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "\${pattern}" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > "\${output}"
        echo "Merged \${label}:" >> \${MERGED}
        ls -l "\${output}" | awk '{print \$5"\t"\$9}' >> \${MERGED}
    }

    write_download_error() {
        local source="\$1"  # "SRA and ENA" or "NCBI Assembly"
        local error_file="\$2"
        echo "Unable to download ${prefix} from \${source} ${task.ext.max_retry} times. This may or may not be a temporary connection issue. Rather than stop the whole Bactopia run, further analysis of ${prefix} will be discontinued." > "\${error_file}"
    }

    validate_single() {
        local fq="\$1"
        gzip -cd "\${fq}" | fastq-scan > r1.json

        if [[ -s r1.json ]]; then
            bactopia-check-fastqs --fq1 r1.json \${OPTS} || { rm -f r1.json; return 1; }
        else
            check_gzip_or_empty "\${fq}"
            rm -f r1.json
            return 1
        fi

        rm -f r1.json
        return 0
    }

    #==========================================================================================
    # Gather or simulate reads based on runtype
    #==========================================================================================
    mkdir -p fastqs assembly
    case "${runtype}" in
        paired-end)
            # Paired-End Reads
            cp -L ${r1_first} fastqs/${prefix}_R1.fastq.gz
            cp -L ${r2_first} fastqs/${prefix}_R2.fastq.gz
            ;;
        single-end)
            # Single-End Reads
            cp -L ${se_first} fastqs/${prefix}_SE.fastq.gz
            ;;
        ont)
            # Long reads (Nanopore)
            cp -L ${lr_first} fastqs/${prefix}_ONT.fastq.gz
            ;;
        hybrid|short_polish)
            # Paired-End Reads + Long Reads
            cp -L ${r1_first} fastqs/${prefix}_R1.fastq.gz
            cp -L ${r2_first} fastqs/${prefix}_R2.fastq.gz
            cp -L ${lr_first} fastqs/${prefix}_ONT.fastq.gz
            ;;
        merge-pe|hybrid-merge-pe|short_polish-merge-pe)
            # Merged Paired-End Reads
            echo "This sample had reads merged." > \${MERGED}
            merge_runs "R1" "*-r1" "fastqs/${prefix}_R1.fastq.gz"
            merge_runs "R2" "*-r2" "fastqs/${prefix}_R2.fastq.gz"
            if [ "${runtype}" == "hybrid-merge-pe" ] || [ "${runtype}" == "short_polish-merge-pe" ]; then
                cp -L ${lr_first} fastqs/${prefix}_ONT.fastq.gz
            fi
            ;;
        merge-se)
            # Merged Single-End Reads
            echo "This sample had reads merged." > \${MERGED}
            merge_runs "SE" "*-se" "fastqs/${prefix}_SE.fastq.gz"
            ;;
        sra_accession|sra_accession_ont)
            if [ "${task.attempt}" == "${task.ext.max_retry}" ]; then
                write_download_error "SRA and ENA" "${prefix}-fastq-download-error.txt"
                exit
            else
                # Download accession from ENA/SRA
                fastq-dl \\
                    --accession ${prefix} \\
                    --provider ${archive} \\
                    --cpus ${task.cpus} \\
                    --outdir fastqs/ \\
                    --prefix ${prefix} \\
                    --group-by-experiment

                # Rename single FASTQ files based on runtype
                if [ -f "fastqs/${prefix}.fastq.gz" ]; then
                    if [ "${runtype}" == "sra_accession_ont" ]; then
                        # ONT reads
                        mv fastqs/${prefix}.fastq.gz fastqs/${prefix}_ONT.fastq.gz
                    else
                        # Single-End
                        mv fastqs/${prefix}.fastq.gz fastqs/${prefix}_SE.fastq.gz
                    fi
                fi
            fi
            ;;
        assembly_accession|assembly)
            if [ "${runtype}" == "assembly_accession" ]; then
                if [ "${task.attempt}" == "${task.ext.max_retry}" ]; then
                    write_download_error "NCBI Assembly" "${prefix}-assembly-download-error.txt"
                    exit
                fi

                # Verify Assembly accession
                bactopia-check-assembly-accession ${prefix} > accession.txt 2> check-assembly-accession.txt

                if [ ! -s "accession.txt" ]; then
                    mv check-assembly-accession.txt ${prefix}-assembly-accession-error.txt
                    exit
                else
                    rm check-assembly-accession.txt
                fi

                # Download from NCBI assembly and simulate reads
                ncbi-genome-download bacteria -o ./ -F fasta -p ${task.cpus} \\
                                            -u "https://ftp.ncbi.nlm.nih.gov/genomes" \\
                                            -s ${section} -A accession.txt -r 50 ${no_cache}

                # Nested directories are not easy to predict, but there should only be a
                # single assembly file. The assembly version (e.g., GCF_000005845.2 --> .2)
                # is removed for consistency.
                find . -name "*${prefix}*.fna.gz" | xargs -I {} mv {} assembly/
                rename 's/(GC[AF]_\\d+).*/\$1.fna.gz/' assembly/*
                gzip -cd assembly/${prefix}.fna.gz > ${prefix}-art.fna
            elif [ "${runtype}" == "assembly" ]; then
                if [ "${meta.is_compressed}" == "true" ]; then
                    gzip -cd ${fna_first} > ${prefix}-art.fna
                else
                    cat ${fna_first} > ${prefix}-art.fna
                fi
                # Keep the original assembly
                cp -L ${fna_first} assembly/${prefix}.fna.gz
            fi

            # Simulate reads from assembly, reads are 250bp without errors
            art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov ${fcov} -ir 0 -ir2 0 -dr 0 -dr2 0 -rs ${task.ext.sampleseed}\\
                         -na -qL 33 -qU 40 -o ${prefix}_R --id ${prefix} -i ${prefix}-art.fna

            mv ${prefix}_R1.fq fastqs/${prefix}_R1.fastq
            mv ${prefix}_R2.fq fastqs/${prefix}_R2.fastq
            pigz -p ${task.cpus} --fast fastqs/*.fastq
            rm ${prefix}-art.fna
            ;;
        *)
            # This should not happen, in case it does fail the whole pipeline
            echo "Runtype ${runtype} is not recognized. Further analysis of ${prefix} will be discontinued." > ${prefix}-runtype-error.txt
            exit 1
            ;;
    esac

    #==========================================================================================
    # Validate input FASTQs
    #==========================================================================================
    if [ "${task.ext.skip_fastq_check}" == "false" ]; then
        ERROR=0
        # Check paired-end reads have same read counts
        OPTS="--sample ${prefix} --min_basepairs ${task.ext.min_basepairs} --min_reads ${task.ext.min_reads} --min_proportion ${task.ext.min_proportion} --runtype ${runtype}"
        if [ -f  "fastqs/${prefix}_R2.fastq.gz" ]; then
            # Paired-end
            gzip -cd fastqs/${prefix}_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd fastqs/${prefix}_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=fastqs/${prefix}_R1.fastq.gz in2=fastqs/${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                echo "${prefix} FASTQs contains an error. Please check the input FASTQs. Further analysis is discontinued." >> ${prefix}-paired-end-error.txt
                ERROR=1
            else
                rm -f ${prefix}-paired-end-error.txt
            fi

            if [[ -s r1.json ]] && [[ -s r2.json ]]; then
                bactopia-check-fastqs --fq1 r1.json --fq2 r2.json \${OPTS} || ERROR=1
            else
                check_gzip_or_empty "fastqs/${prefix}_R1.fastq.gz" || ERROR=1
                check_gzip_or_empty "fastqs/${prefix}_R2.fastq.gz" || ERROR=1
            fi
            rm -f r1.json r2.json
        elif [ -f "fastqs/${prefix}_SE.fastq.gz" ]; then
            # Single-end
            validate_single "fastqs/${prefix}_SE.fastq.gz" || ERROR=1
        elif [ -f "fastqs/${prefix}_ONT.fastq.gz" ]; then
            # Long reads
            validate_single "fastqs/${prefix}_ONT.fastq.gz" || ERROR=1
        fi

        # Failed validations so, let's keep them from continuing
        if [ "\${ERROR}" -eq "1" ]; then
            mv fastqs/ failed-tests-fastqs/
        fi
    fi

    #==========================================================================================
    # Determine the paired/single-end status
    #==========================================================================================
    IS_PAIRED="false"
    if [ -f "fastqs/${prefix}_R2.fastq.gz" ] && [ "${runtype}" != "short_polish" ]; then
        IS_PAIRED="true"
    fi

    #==========================================================================================
    # Create metadata file
    #==========================================================================================
    printf "sample\\truntype\\toriginal_runtype\\tis_paired\\tis_compressed\\tspecies\\tgenome_size\\n" > ${prefix}-meta.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${prefix}" "${meta.runtype}" "${meta.original_runtype}" "\$IS_PAIRED" "${meta.is_compressed}" "${meta.species}" "${meta.genome_size}" >> ${prefix}-meta.tsv

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(echo \$(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        fastq-dl: \$(echo \$(fastq-dl --version 2>&1) | sed 's/fastq-dl, version //')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        ncbi-genome-download: \$(echo \$(ncbi-genome-download --version 2>&1))
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    """
}
