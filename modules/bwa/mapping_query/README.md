# mapping_query process testing:

This process maps FASTQ reads against a given set of FASTA files using BWA.

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run mappingg_query.nf -profile test,docker -params-file test_params.yaml -entry test


if you've used `bactopia conda activate` you can also trade `docker` by `conda` to test with conda. 
