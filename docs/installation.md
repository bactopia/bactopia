# Installation
## Requirements
The main requirement for using Bactopia is to have Conda installed and setup.
If you aren't familiar with Conda, it [is an open source package management
system and environment management system that runs on Windows, macOS and
Linux](https://conda.io/en/latest/).

If you need help setting up Conda, here is a link to the official Conda
documentation for
*[Installing Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)*.
Bactopia has been tested using the
[Miniconda installer](https://conda.io/en/latest/miniconda.html),
but the [Anaconda installer](https://www.anaconda.com/distribution/) should
work the same.


## Local Install
Once you have Conda all set up, you are ready to create an environment for
Bactopia.

### Setting Up A Bactopia Conda Environment
The Conda envrionment configuration file
[bactopia.yml](https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml)
has been created to help ensure the minimum requirements to run Bactopia are
installed.

#### Download *bactopia.yml*
First you will want to download the `bactopia.yml` file.

``` bash
wget https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml

- or -

curl -o bactopia.yml https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml
```

Or, feel free to "Right-click Save As" [**bactopia.yml**](https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml)

#### Create A Bactopia Environment

``` bash
conda env create -n bactopia -f bactopia.yml
```

!!! warning "Possible warning message, you can ignore it!"

    Warning: you have pip-installed dependencies in your environment file, but
    you do not list pip itself as one of your conda dependencies. Conda may
    not use the correct pip to install your packages, and they may end up in
    the wrong place. Please add an explicit pip dependency. I'm adding one for
    you, but still nagging you.

Pip is in the `bactopia.yml`, but for some reason its not recoginzed! Again, it
is safe to disregard this warning, as the additional Python packages should still
be installed correctly.

#### Clone Bactopia
Now that a Bactopia environment has been created. You can activate it and clone the Bactopia repo.

``` bash
conda activate bactopia  # or, source activate bactopia
git clone git@github.com:staphopia/bactopia-ap.git
cd bactopia-ap

./bactopia.nf --version
N E X T F L O W  ~  version 19.01.0
Launching `./bactopia.nf` [cheesy_carson] - revision: 97a14660a2
bactopia 0.0.1
```

Everything *should* be all set! (Is it ever?) When you execute Bactopia, seperate conda environments will be created for each process automatically by Nextflow.

The nextstep is to [Build Datasets](datasets.md) for usage by Bactopia.

## Coming Soon
As the analysis pipeline for Bactopia is finalized there are plans to also
include [Nextflow Core](https://nf-co.re/) installations and containers via
[Docker](https://www.docker.com/) and [Singularity](https://www.sylabs.io/singularity/)
