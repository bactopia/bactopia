# download_references process testing:

This process downloads the nearest RefSeq genomes (based on Mash) to have variants called against.

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run download_references.nf -params-file test_params.yaml -profile test,docker -entry test


if you've used `bactopia conda activate` you can also trade `docker` by conda to test with conda. 

