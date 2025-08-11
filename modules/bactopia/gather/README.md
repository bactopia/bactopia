# gather_fastqs process testing:

This process handles the input files into channels for other process in the workflow.

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run gather_fastqs.nf -params-file test_params.yaml -profile test,docker -entry test


if you've used `bactopia conda activate` you can also trade `docker` by conda to test with conda. 