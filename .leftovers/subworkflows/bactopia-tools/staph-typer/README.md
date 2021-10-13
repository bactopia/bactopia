# Bactopia Tools - *staph-typer*
The `staph-typer` tool includes multiple tools that are specific for typing certain features of *Staphylococcus aureus*. Currently `staph-typer` includes

1. [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) - *agr* locus type and *agr* operon variants.
2. [spaTyper](https://github.com/HCGB-IGTP/spaTyper) - *spa* type
3. [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec) - SCCmec type

This tool will evolve with *S. aureus* genomics, so you can expect it to add more typing methods (maybe even replace current methods) in the future. If a certain typing method for *S. aureus* please feel free to suggest it be added!~


## Example
The following command will run `staph-typer` on each available sample.
```
bactopia tools staph-typer --bactopia ~/bactopia-tutorial/bactopia
```

## Output Overview
Below is the default output structure for the `staph-typer` tool. Where possible the 
file descriptions below were modified from a tools description.
```
bactopia-tools/
└── staph-typer/
    └── ${PREFIX}/
        ├── bactopia-info
        │   ├── staph-typer-report.html
        │   ├── staph-typer-timeline.html
        │   └── staph-typer-trace.txt
        └── ${SAMPLE_NAME}
        │  ├── agrvate
        │  │   ├── ${SAMPLE_NAME}-agr_gp.tab
        │  │   ├── ${SAMPLE_NAME}-agr_operon.fna
        │  │   ├── ${SAMPLE_NAME}-agr_operon_frameshifts.tab
        │  │   ├── ${SAMPLE_NAME}-blastn_log.txt
        │  │   ├── ${SAMPLE_NAME}-mummer/
        │  │   ├── ${SAMPLE_NAME}-mummer-log.txt
        │  │   ├── ${SAMPLE_NAME}-snippy/
        │  │   ├── ${SAMPLE_NAME}-snippy-log.txt
        │  │   └── ${SAMPLE_NAME}-summary.tab
        │  ├── ${SAMPLE_NAME}-spatyper.txt
        │  └── ${SAMPLE_NAME}-sccmec.txt
        ├── agrvate-results.txt
        ├── spatyper-results.txt
        └── sccmec-results.txt
```

Below is a description of `staph-typer` outputs.

| Filename | Description |
|-----------|-------------|
| agrvate-results.txt| Merged set of outputs from `AgrVATE` |
| spatyper-results.txt | Merged set of outputs from `spaTyper` |
| sccmec-results.txt | Merged set of outputs from `staphopia-sccmec` |


### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| staph-typer-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| staph-typer-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| staph-typer-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### Per Sample

| Filename | Description |
|-----------|-------------|
| ${SAMPLE_NAME}-spatyper.txt | Predicted *spa* type from `spaTyper` |
| ${SAMPLE_NAME}-sccmec.txt | Predicted SCCmec type from `staphopia-sccmec` |

##### agrvate
`AgrVATE` includes outputs from multiple programs so it gets a separate directory within the sample directory. 

| Extension | Description |
|----------|-------------|
| -agr_gp.tab | Detailed report for agr kmer matches |
| -agr_operon.fna | Agr operon extracted from in-silico PCR |
| -agr_operon_frameshifts.tab | Frameshift mutations in CDS of extracted agr operon detected by Snippy |
| log.txt | Log files from programs called by `AgrVATE` |
| -{mummer|snippy}/ | Intermediate files from mummer and snippy |
| -summary.tab | A final summary report for agr typing |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

AgrVATE Parameters:
    --typing_only           Does agr typing only. Skips agr operon extraction and frameshift
                                detection.

spaTyper Parameters:
    --do_enrich             Do PCR product enrichment

staphopia-sccmec Parameters:
    --hamming               Report the results as hamming distances.
                                Default: True (perfect match) or False (at least one mismatch)

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: staph-typer

    --outdir DIR            Directory to write results to
                                Default: ./

    --min_time INT          The minimum number of minutes a job should run before being halted.
                                Default: 60 minutes

    --max_time INT          The maximum number of minutes a job should run before being halted.
                                Default: 120 minutes

    --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                Default: 32 Gb

    --cpus INT              Number of processors made available to a single
                                process.
                                Default: 1

Nextflow Related Parameters:
    --condadir DIR          Directory to Nextflow should use for Conda environments
                                Default: Bactopia's Nextflow directory

    --registry STR          Docker registry to pull containers from.
                                Available options: dockerhub, quay, or github
                                Default: dockerhub

    --singularity_cache STR Directory where remote Singularity images are stored. If using a cluster, it must
                                be accessible from all compute nodes.
                                Default: NXF_SINGULARITY_CACHEDIR evironment variable, otherwise /local/home/rpetit/.bactopia/singularity

    --queue STR             The name of the queue(s) to be used by a job scheduler (e.g. AWS Batch or SLURM).
                                If using multiple queues, please seperate queues by a comma without spaces.
                                Default: general

    --disable_scratch       All intermediate files created on worker nodes of will be transferred to the head node.
                                Default: Only result files are transferred back

    --cleanup_workdir       After Bactopia is successfully executed, the work directory will be deleted.
                                Warning: by doing this you lose the ability to resume workflows.

    --publish_mode          Set Nextflow's method for publishing output files. Allowed methods are:
                                'copy' (default)    Copies the output files into the published directory.

                                'copyNoFollow' Copies the output files into the published directory
                                               without following symlinks ie. copies the links themselves.

                                'link'    Creates a hard link in the published directory for each
                                          process output file.

                                'rellink' Creates a relative symbolic link in the published directory
                                          for each process output file.

                                'symlink' Creates an absolute symbolic link in the published directory
                                          for each process output file.

                                Default: copy

    --force                 Nextflow will overwrite existing output files.
                                Default: true

    --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                will wait before execution.
                                Default: 5 seconds

    --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows
                                you to create profiles specific to your environment (e.g. SGE,
                                AWS, SLURM, etc...). This config file is loaded last and will
                                overwrite existing variables if set.
                                Default: Bactopia's default configs

    -resume                 Nextflow will attempt to resume a previous run. Please notice it is
                                only a single '-'

AWS Batch Related Parameters:
    --aws_region STR        AWS Region to be used by Nextflow
                                Default: us-east-1

    --aws_volumes STR       Volumes to be mounted from the EC2 instance to the Docker container
                                Default: /opt/conda:/mnt/conda

    --aws_cli_path STR       Path to the AWS CLI for Nextflow to use.
                                Default: /home/ec2-user/conda/bin/aws

    --aws_upload_storage_class STR
                            The S3 storage slass to use for storing files on S3
                                Default: STANDARD

    --aws_max_parallel_transfers INT
                            The number of parallele transfers between EC2 and S3
                                Default: 8

    --aws_delay_between_attempts INT
                            The duration of sleep (in seconds) between each transfer between EC2 and S3
                                Default: 15

    --aws_max_transfer_attempts INT
                            The maximum number of times to retry transferring a file between EC2 and S3
                                Default: 3

    --aws_max_retry INT     The maximum number of times to retry a process on AWS Batch
                                Default: 4

    --aws_ecr_registry STR  The ECR registry containing Bactopia related containers.
                                Default: Use the registry given by --registry

Useful Parameters:
    --version               Print workflow version information
    --help                  Show this message and exit
```

## Conda Environment
Below is the command used to create the Conda environment.
```
conda create -n staph-typer -c conda-forge -c bioconda agrvate spatyper staphopia-sccmec 'snippy>=4.5.0'
```

## References
* __[AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE)__  
Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.  
_Raghuram V., [AgrVATE: Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.](https://github.com/VishnuRaghuram94/AgrVATE)_  

* __[spaTyper](https://github.com/HCGB-IGTP/spaTyper)__  
Computational method for finding spa types.  
_Harmsen D., Claus H., Witte W., Rothgänger J., Claus H., Turnwald D., and Vogel U.. [Typing of methicillin-resistant Staphylococcus aureus in a university hospital setting using a novel software for spa-repeat determination and database management.](https://doi.org/10.1128/jcm.41.12.5442-5448.2003) J. Clin. Microbiol. 41:5442-5448 (2003)._  
_Sanchez-Herrero J.F., & mjsull. (2020, October 2). [spaTyper: Staphylococcal protein A (spa) characterization pipeline](http://doi.org/10.5281/zenodo.4063625). Zenodo._  

* __[staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec)__  
A standalone version of Staphopia's SCCmec typing method.  
_Petit III R.A., Read T.D., [Staphylococcus aureus viewed from the perspective of 40,000+ genomes.](http://dx.doi.org/10.7717/peerj.5261) PeerJ 6, e5261 (2018)._  
