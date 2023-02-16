/*
========================================================================================
    Nextflow Functions specific to Bactopia Tools
========================================================================================
/*

/*
========================================================================================
    Public: collect_samples(bactopia_dir, extension, exclude_list, include_list)

    bactopia_dir -> the path to bactopia outputs
    extension -> the file extension to collect (e.g. 'fastq', 'fna', gff', ...)
    exclude_list -> a text file with a subset of samples to exclude from analysis
    include_list -> a text file with a subset of samples to only include from analysis
========================================================================================
*/
def collect_samples(bactopia_dir, extension, include_list, exclude_list) {
    include_all = true
    inclusions = []
    exclusions = []
    IGNORE_LIST = ['.nextflow', 'bactopia-info', 'bactopia-tools', 'work',]
    if (include_list) {
        new File(include_list).eachLine { line -> 
            inclusions << line.trim().split('\t')[0]
        }
        include_all = false
        log.info "Including ${inclusions.size} samples for analysis"
    }
    else if (exclude_list) {
        new File(exclude_list).eachLine { line -> 
            exclusions << line.trim().split('\t')[0]
        }
        log.info "Excluding ${exclusions.size} samples from the analysis"
    }
    
    sample_list = []
    log.info("${bactopia_dir}/bactopia-samples")
    file("${bactopia_dir}/bactopia-samples").eachFile { item ->
        if( item.isDirectory() ) {
            sample = item.getName()
            if (!IGNORE_LIST.contains(sample)) {
                if (inclusions.contains(sample) || include_all) {
                    if (!exclusions.contains(sample)) {
                        if (_is_sample_dir(sample, bactopia_dir)) {
                            sample_list << _collect_inputs(sample, bactopia_dir, extension)
                        } else {
                            log.info "${sample} does not appear to be a Bactopia sample, skipping..."
                        }
                    }
                }
            }
        }
    }

    log.info "Found ${sample_list.size} samples to process"
    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for 5 seconds..."
    log.info "--------------------------------------------------------------------"
    sleep(5000)
    return sample_list
}

/*
========================================================================================
    Public: collect_samples(file_path)

    key -> parameter that represents the filr
    local_file -> path to the files to collect
    local_pattern -> Pattern to match incase `local_file` is a directory
========================================================================================
*/
def collect_local_files(local_file, local_pattern) {
    local_list = []
    if (local_file) {
        if (file(local_file).exists()) {
            if (file(local_file).isDirectory()) {
                file("${local_file}/${local_pattern}" ).each{ 
                    local_list << tuple([[id: it.getSimpleName()], it])                    
                }
            } else {
                local_list << tuple([id: file(local_file).getSimpleName()], file(local_file))
            }
        }
    }
    return local_list
}

/*
========================================================================================
    Private: _is_sample_dir(sample, dir)

    This function does a quick check to see if the *-genome-size.txt file exists, which
    is a required output and all samples processed by Bactopia should have it.

    sample -> The name of the sample
    dir -> The base directory of sample results
========================================================================================
*/
def _is_sample_dir(sample, dir) {
    return file("${dir}/bactopia-samples/${sample}").exists()
}

/*
========================================================================================
    Private: _collect_inputs(sample, dir, extension)

    This function collects outputs based on file extension. As more outputs are needed 
    the PATHS map will need to be updated.

    sample -> The name of the sample
    dir -> The base directory of sample results
    extension -> The file extension to match.
========================================================================================
*/
def _collect_inputs(sample, dir, extension) {
    PATHS = [:]
    PATHS.fastq = "qc"
    PATHS.fna = "assembler"
    PATHS.faa = "annotator"
    PATHS.gff = "annotator"

    base_dir = "${dir}/bactopia-samples/${sample}/bactopia-main/"
    se = "${base_dir}/${PATHS['fastq']}/${sample}.fastq.gz"
    ont = "${base_dir}/${PATHS['fastq']}/.ont"
    pe1 = "${base_dir}/${PATHS['fastq']}/${sample}_R1.fastq.gz"
    pe2 = "${base_dir}/${PATHS['fastq']}/${sample}_R2.fastq.gz"
    fna = "${base_dir}/${PATHS['fna']}/${sample}.fna"

    if (extension == 'fastq') {
        if (file(se).exists()) {
            if (file(ont).exists()) {
                return tuple([id:sample, single_end:true, runtype:'ont'], [file(se)])
            } else {
                return tuple([id:sample, single_end:true, runtype:'illumina'], [file(se)])
            }
        } else if (file(pe1).exists() && file(pe2).exists()) {
            return tuple([id:sample, single_end:false, runtype:'illumina'], [file(pe1), file(pe2)])
        } else {
            log.error("Could not locate FASTQs for ${sample}, please verify existence. Unable to continue.")
            exit 1
        }
    } else if (extension == 'fna_fastq') {
        if (file(se).exists()) {
            runtype = "illumina"
            if (file(ont).exists()) {
                runtype = "ont"
            }

            if (file("${fna}.gz").exists()) {
                return tuple([id:sample, single_end:true, is_compressed:true, runtype:runtype], [file("${fna}.gz")], [file(se)])
            } else {
                return tuple([id:sample, single_end:true, is_compressed:false, runtype:runtype], [file("${fna}")], [file(se)])
            }
        } else if (file(pe1).exists() && file(pe2).exists()) {
            if (file("${fna}.gz").exists()) {
                return tuple([id:sample, single_end:false, is_compressed:true, runtype:'illumina'], [file("${fna}.gz")], [file(pe1), file(pe2)])
            } else {
                return tuple([id:sample, single_end:false, is_compressed:false, runtype:'illumina'], [file("${fna}")], [file(pe1), file(pe2)])
            }
        } else {
            log.error("Could not locate FASTQs for ${sample}, please verify existence. Unable to continue.")
            exit 1
        }
    } else if (extension == 'fna_faa') {
        // Default to Bakta faa
        faa = "${base_dir}/${PATHS['faa']}/bakta/${sample}.faa"
        if (!file("${faa}").exists() && !file("${faa}.gz").exists()) {
            // Fall back on Prokka
            faa = "${base_dir}/${PATHS['faa']}/prokka/${sample}.faa"
        }

        if (file("${fna}.gz").exists() && file("${faa}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], [file("${fna}.gz")], [file("${faa}.gz")])
        } else {
            return tuple([id:sample, is_compressed:false], [file("${fna}")], [file("${faa}")])
        }
    } else {
        input = "${base_dir}/${PATHS[extension]}/${sample}.${extension}"
        if (extension == "faa" || extension == "gff") {
            // Default to Bakta faa
            input = "${base_dir}/${PATHS[extension]}/bakta/${sample}.${extension}"
            if (!file("${input}").exists() && !file("${input}.gz").exists()) {
                // Fall back on Prokka
                input = "${base_dir}/${PATHS[extension]}/prokka/${sample}.${extension}"
            }
        }

        if (file("${input}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], [file("${input}.gz")])
        } else if (file(input).exists()) {
            return tuple([id:sample, is_compressed:false], [file("${input}")])
        } else {
            log.error("Could not locate ${input} for ${sample}, please verify existence. Unable to continue.")
            exit 1
        }
    }
}
