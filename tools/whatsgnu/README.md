# Bactopia Tools - *whatsgnu*


## Example
The following command will run `whatsgnu` on each available sample.
```

```

## Output Overview
Below is the default output structure for the `whatsgnu` tool. Where possible the 
file descriptions below were modified from a tools description.

```

```

| Filename | Description |
|-----------|-------------|
|  |  |
|  |  |

### Directory Description
#### bactopia-info
| Filename | Description |
|----------|-------------|
| whatsgnu-report.html | The Nextflow [Execution Report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) |
| whatsgnu-timeline.html | The Nextflow [Timeline Report](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) |
| whatsgnu-trace.txt | The Nextflow [Trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) report |


## Usage
```

```

## Conda Environment
Below is the command used to create the Conda environment.
```
conda create -n bactopia-whatsgnu -c conda-forge -c bioconda whatsgnu ncbi-genome-download wget
```

## References

* __[WhatsGNU](https://github.com/ahmedmagds/WhatsGNU)__  
_A.M. Moustafa, P.J. Planet, [WhatsGNU: a tool for identifying proteomic novelty.](https://doi.org/10.1186/s13059-020-01965-w) Genome Biol. 21, 58 (2020)._  

* __[ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)__  
_Blin, K. [ncbi-genome-download: Scripts to download genomes from the NCBI FTP 
servers](https://github.com/kblin/ncbi-genome-download)_  
