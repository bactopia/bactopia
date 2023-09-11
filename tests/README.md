# Testing in Bactopia

## Acknowledgements

Before I get started walking everyone through testing in Bactopia, I would first like to extend a huge: *Thank you!* to the [nf-core Team](https://nf-co.re/community). I attended a nf-core Hackathon, and was introduced to [nf-core/modules](https://github.com/nf-core/modules) and it's usage of `pytest` for testing each module. Instantly I took a liking to it, and if you take a look at Bactopia's testing setup, it is very much adapted from `nf-core/modules`. So, again big thank you to the `nf-core` team! You definitely provided a great model for testing, that allowed me to over a weekend get set up!

## `bactopia-tests` Data

Taking a cue from [nf-core/test-datasets](https://github.com/nf-core/test-datasets), a repository housing basic test data for Bactopia was created at [bactopia/bactopia-tests](https://github.com/bactopia/bactopia-tests). `bactopia-tests` includes most of the data necessary to test each Bactopia module. As you might imagine, there are some datasets (*e.g. RefSeq Mash sketch (~700mb)*) that are much to large to include in a repo. To keep the `bactopia-tests` repo small, I included seqeuncing for *Candidatus Portiera aleyrodidarum* which has a completed genome of ~350kb (*yep! its super small!*). This allows modules to be tested with real data, without long run times or extensive compute resources.

As testing expands, more test data will be added to `bactopia-tests`. Please keep in mind, much like nf-core's `test-datasets` repo, you are more than welcome to use the data in `bactopia-tests` for your own testing.

## Bactopia CI Setup

*Below are mostly notes for myself (Robert), when updating and working with the CI. If you make use of this, you will need to make modifications to fit your system.*

The testing in Bactopia allows you to test the `conda`, `docker`, or `singularity` profiles using `pytest` and `pytest-workflows`. Here are notes on how things have been set up on a self-hosted GHA runner. Worth noting also, [Mambaforge](https://github.com/conda-forge/miniforge) is being used for Conda.


```{bash}
# If no screen
screen -S bactopia-ci
# else
screen -Dr bactopia-ci

# Move to working directory
cd /data/storage/bactopia-ci

# Rebuild the bactopia-dev environment
mamba create -n bactopia-dev -y -c rpetit3 -c conda-forge -c bioconda bactopia
conda activate bactopia-dev

# Clone Bactopia repo
rm -rf bactopia/
git clone --branch dev git@github.com:bactopia/bactopia.git
cd bactopia/
git pull
cd ..

# Download necessary environments (our case all of them!)
# This will take a while, but its only needed occasionally
mkdir -p envs
cd envs
bactopia-download \
    --bactopia-path /data/storage/bactopia-ci/bactopia \
    --envtype all \
    --build-all \
    --condadir /data/storage/bactopia-ci/envs/conda \
    --singularity_cache /data/storage/bactopia-ci/envs/singularity \
    --verbose

# 

# Clone Bactopia-Tests repo (optional, but saves downloading each test)
git clone git@github.com:bactopia/bactopia-tests.git
```

I think that should be all you need to test things out.

## Testing Overview

Bactopia uses `pytest` and `pytest-workflow` to run the tests and determine if the expected files are output (and if applicable the `md5sum` matches). For each module, the structure looks like this:

```{bash}
<module_name>
├── README.md
├── bin
│   └── check-staging.py
├── main.nf
├── tests.nf
└── tests.yml
```

The two files, `tests.nf` and `tests.yml` (haha as if their names didn't give them away), are used for testing. For demonstration purposes, I'm going to use tests from the `minmer_sketch` module.

### `tests.nf`

The `tests.nf` file is a Nextflow workflow that provides the necessaty inputs to run the module.

```{nextflow}
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MINMER_SKETCH } from './main.nf'

workflow test_minmer_sketch_pe {
    inputs = tuple(
        "test_minmer_sketch_pe",
        false,
        [file(params.test_data['illumina']['r1'], checkIfExists: true), file(params.test_data['illumina']['r2'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}

workflow test_minmer_sketch_se {
    inputs = tuple(
        "test_minmer_sketch_se",
        true,
        [file(params.test_data['illumina']['se'], checkIfExists: true)]
    )

    MINMER_SKETCH ( inputs )
}
```

In this example above, there are two tests `test_minmer_sketch_pe` and `test_minmer_sketch_se`. As you might have guessed the `test_minmer_sketch_pe` tests the `minmer_sketch` module with paired-end reads, and `test_minmer_sketch_se` with single-end reads. That's pretty much it, when it comes to the `tests.nf`.

### `tests.yml`

The `tests.yml` file contains the command that `pytest` will run along with the expected outputs and if reproducible, their `md5sum` hash.

```{yaml}
- name: minmer_sketch_pe
  command: bash -c 'nextflow run ./modules/bactopia/minmer_sketch/tests.nf -entry test_minmer_sketch_pe -c nextflow.config ${BACTOPIA_ARGS}'
  tags:
    - minmer_sketch_pe
  files:
    - path: test_minmer_sketch_pe/logs/minmer_sketch/minmer_sketch.err
      md5sum: b09af7ccf1a4f6766cd1b217c01a9fa7
    - path: test_minmer_sketch_pe/logs/minmer_sketch/minmer_sketch.out
      md5sum: d41d8cd98f00b204e9800998ecf8427e
    - path: test_minmer_sketch_pe/logs/minmer_sketch/minmer_sketch.sh
      md5sum: fda1ca7d2ee3a47ecbb9980eeae1b6c6
    - path: test_minmer_sketch_pe/logs/minmer_sketch/minmer_sketch.trace
      md5sum: d41d8cd98f00b204e9800998ecf8427e
    - path: test_minmer_sketch_pe/logs/minmer_sketch/minmer_sketch.versions
    - path: test_minmer_sketch_pe/minmers/test_minmer_sketch_pe-k21.msh
      md5sum: c8bcc71befab49c30c5fcdaa432a2b79
    - path: test_minmer_sketch_pe/minmers/test_minmer_sketch_pe-k31.msh
      md5sum: 233521f8fe7eb56a40beba081f0c0842
    - path: test_minmer_sketch_pe/minmers/test_minmer_sketch_pe.sig
      md5sum: d433f60ebbc0c332e085261398c6d37d

- name: minmer_sketch_se
  command: bash -c 'nextflow run ./modules/bactopia/minmer_sketch/tests.nf -entry test_minmer_sketch_se -c nextflow.config ${BACTOPIA_ARGS}'
  tags:
    - minmer_sketch_se
  files:
    - path: test_minmer_sketch_se/logs/minmer_sketch/minmer_sketch.err
      md5sum: 63a12f74cdd2a65c4cd186576edaceb0
    - path: test_minmer_sketch_se/logs/minmer_sketch/minmer_sketch.out
      md5sum: d41d8cd98f00b204e9800998ecf8427e
    - path: test_minmer_sketch_se/logs/minmer_sketch/minmer_sketch.sh
      md5sum: 475d824115605d13f01d68711dad4c5e
    - path: test_minmer_sketch_se/logs/minmer_sketch/minmer_sketch.trace
      md5sum: d41d8cd98f00b204e9800998ecf8427e
    - path: test_minmer_sketch_se/logs/minmer_sketch/minmer_sketch.versions
    - path: test_minmer_sketch_se/minmers/test_minmer_sketch_se-k21.msh
      md5sum: af17317885feeaeb47d183a465211451
    - path: test_minmer_sketch_se/minmers/test_minmer_sketch_se-k31.msh
      md5sum: 6d46feb639587bc57dd28cd4fb855eb0
    - path: test_minmer_sketch_se/minmers/test_minmer_sketch_se.sig
      md5sum: dd9e4ae787e4eb5f6fccad820fbf853b
```

I think most of this we can wrap our head around, just but looking at the contents. But what I want to go into detail about is the `command` line.

```{yaml}}
command: bash -c 'nextflow run ./modules/bactopia/minmer_sketch/tests.nf -entry test_minmer_sketch_pe -c nextflow.config ${BACTOPIA_ARGS}'
```

What this is doing is running Nextflow using the default Bactopia config, to run the `tests.nf` for `minmer_sketch` with an entry piont at `test_minmer_sketch_pe`. But there is also this `${BACTOPIA_ARGS}` but at the end. Now what this does is allows you to provide command line arguments to the (also explains why the `nextflow run` is wrapped with `bash -c`!)

### Using `${BACTOPIA_ARGS}`

`${BACTOPIA_ARGS}` follows, nf-core/modules' `${PROFILE}` parameter, which would tell Nextflow when profile to use (`conda`, `docker`, or `singularity`). I instead took it a step further, since Bactopia includes numerous adjustable parameters I wanted to be able to change them, but I think most imporantly, I wanted to be able to test Bactopia using its default configs (*e.g. params.config*) into order to better replicate the user experience.

*How does this work?*

All you need to do is prepend your `pytest` command with `BACTOPIA_ARGS="...parameters..."` and the contents of `${BACTOPIA_ARGS}` will be appended to the end of the command in the `tests.yml`.

#### Examples

```{bash}
# Use Singularity
BACTOPIA_ARGS="-profile singularity" pytest --symlink --keep-workflow-wd --tag minmer_sketch_se
# becomes ---> nextflow run ./modules/bactopia/minmer_sketch/tests.nf -entry test_minmer_sketch_se -profile singularity

# Point to local bactopia-tests
BACTOPIA_ARGS="--test_data_dir /home/robert_petit/repos/bactopia-tests/data" pytest --symlink --keep-workflow-wd --tag minmer_sketch_se
# becomes ---> nextflow run ./modules/bactopia/minmer_sketch/tests.nf -entry test_minmer_sketch_se --test_data_dir /home/robert_petit/repos/bactopia-tests/data
```

As you might imagine, this allows you to completely control the parameters that are available in Bactopia.
