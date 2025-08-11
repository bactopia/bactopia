# assemble_genome process testing:

This process assemble the genome using Shovill, SKESA is used by default

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run asssemble_genome.nf -params-file test_params.yaml -profile test,docker -entry test


if you've used `bactopia conda activate` you can also trade `docker` by conda to test with conda. 

