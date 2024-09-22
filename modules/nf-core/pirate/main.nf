// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'pirate')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::pirate=1.0.5 bioconda::perl-bioperl=1.7.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PIRATE {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate%3A1.0.5--hdfd78af_0' :
        'quay.io/biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("core-genome.aln.gz")               , emit: aln, optional: true
    tuple val(meta), path("results/gene_presence_absence.csv"), emit: csv, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs -r gunzip

    # PIRATE only supports .gff extension, will need to adjust for gff3 (Bakta) files
    # https://github.com/SionBayliss/PIRATE/blob/master/scripts/run_PIRATE.pl#L153
    # note for later: "xargs -r" will not run if no files are found
    find gff/ -name "*.gff3.gz" | xargs -r gunzip
    find gff/ -name "*.gff3" -print0 | while read -d \$'\0' file; do mv "\$file" "\${file%.gff3}.gff"; done

    PIRATE \\
        $options.args \\
        --align \\
        --threads $task.cpus \\
        --input ./gff/ \\
        --output results/
    PIRATE_to_roary.pl -i results/PIRATE.*.tsv -o results/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}

    # Only copy files if they exist
    if [[ -f "results/core_alignment.fasta.gz" ]]; then
        cp results/core_alignment.fasta.gz ./core-genome.aln.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
