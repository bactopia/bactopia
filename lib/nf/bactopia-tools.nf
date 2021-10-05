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
def collect_samples(bactopia_dir, extension, exclude_list, include_list) {
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
    file(bactopia_dir).eachFile { item ->
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
    sleep(5000)
    return sample_list
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
    return file("${dir}/${sample}/${sample}-genome-size.txt").exists()
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
    PATHS.fastq = "quality-control"
    PATHS.fna = "assembly"
    PATHS.gff = "annotation"

    if (extension == 'fastq') {
        se = "${dir}/${sample}/quality-control/${sample}.fastq.gz"
        pe1 = "${dir}/${sample}/quality-control/${sample}_R1.fastq.gz"
        pe2 = "${dir}/${sample}/quality-control/${sample}_R2.fastq.gz"

        if (path(se).exists()) {
            return tuple([id:sample, single_end:true], [path(se)])
        } else if (path(pe1).exists() && path(pe2).exists()) {
            return tuple([id:sample, single_end:false], [path(pe1), path(pe2)])
        } else {
            log.error("Could not locate FASTQs for ${sample}, please verify existence. Unable to continue.")
            exit 1
        }
    } else {
        input = "${dir}/${sample}/${PATHS[extension]}/${sample}.${extension}"
        if (file("${input}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], ("${input}.gz")
        } else if (file(assembly).exists()) {
            return tuple([id:sample, is_compressed:false], ("${input}.gz")
        } else {
            log.error("Could not locate ${input} for ${sample}, please verify existence. Unable to continue.")
            exit 1
        }
    }
}
