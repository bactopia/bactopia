# `annotate_genome` module

This process annotate the assembly using Prokka, use a proteins FASTA if available

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
    - the local files in `test_data` 
    - params in  `test_params.yaml`
    - `test` profile in `nextflow.config`

## Testing `annotate_genome`

### Conda

```{bash}
BACTOPIA_ARGS="--test_data_dir /path/to/repos/bactopia-tests/data --condadir /path/to/conda/envs" \
    pytest --tag annotate_genome --symlink --keep-workflow-wd
```

### Docker

```{bash}
BACTOPIA_ARGS="-profile docker --test_data_dir /path/to/repos/bactopia-tests/data --condadir /path/to/conda/envs" \
    pytest --tag annotate_genome --symlink --keep-workflow-wd
```

### Singularity

```{bash}
BACTOPIA_ARGS="-profile singularity --test_data_dir /path/to/repos/bactopia-tests/data --condadir /path/to/conda/envs" \
    pytest --tag annotate_genome --symlink --keep-workflow-wd
```
