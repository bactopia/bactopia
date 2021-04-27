# Bactopia Tools - *hicap*
The `hicap` tool uses [hicap](https://github.com/scwatts/hicap) for the in-silico typing of the *H. influenzae* cap locus.

## Example
The following command will run `hicap` on each available sample.
```
bactopia tools hicap --bactopia ~/bactopia-tutorial/bactopia
```

## Output Overview
Below is the default output structure for the `hicap` tool. Where possible the 
file descriptions below were modified from a [tools description](https://github.com/scwatts/hicap#outputs).
```
bactopia-tools/
└── hicap/
    └── ${PREFIX}/
        ├── bactopia-info
        │   ├── hicap-report.html
        │   ├── hicap-timeline.html
        │   └── hicap-trace.txt
        └── ${SAMPLE_NAME}
        │  ├── ${SAMPLE_NAME}.gbk
        │  ├── ${SAMPLE_NAME}.svg
        │  └── ${SAMPLE_NAME}.tsv 
        └── hicap-results.txt
```

Below is a description of `hicap` outputs.

| Filename | Description |
|-----------|-------------|
| hicap-results.txt| Merged set of outputs from `hicap` |


### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| hicap-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| hicap-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| hicap-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |

#### Per Sample
| Extension | Description |
|-----------|-------------|
| .gbk | GenBank file with sequence marked up with cap locus annotations |
| .svg | visual representation of the annotated cap locus |
| .tsv | detailed summary information |


## Usage
```
Required Parameters:
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

hicap Related Parameters
    --database_dir STR      Directory containing locus database.
                                Default: Use hicap's default

    --model_fp STR    Path to prodigal model.
                                Default: Use hicap's default

    --full_sequence BOOL    Write the full input sequence out to the genbank file rather than
                                just the region surrounding and including the locus

    --gene_coverage FLOAT   Minimum percentage coverage to consider a single gene complete.
                                Default: 0.8

    --gene_identity FLOAT   Minimum percentage identity to consider a single gene complete.
                                Default: 0.7

    --broken_gene_length INT
                            Minimum length to consider a broken gene.
                                Default: 60

    --broken_gene_identity FLOAT
                            Minimum percentage identity to consider a broken gene
                                Default: 0.8

    --log_fp                Record logging messages to file

    --debug                 Print debug messages

Optional Parameters:
    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --prefix DIR            Prefix to use for final output files
                                Default: hicap

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
conda create -n hicap -c conda-forge -c bioconda hicap
```

## References
* __[hicap](https://github.com/scwatts/hicap)__  
in silico typing of the *H. influenzae* cap locus
_Watts S.C. and Holt K.E. [hicap: in silico serotyping of the Haemophilus influenzae capsule locus.](https://doi.org/10.1128/JCM.00190-19) Journal of Clinical Microbiology, JCM.00190-19 (2019). _  
