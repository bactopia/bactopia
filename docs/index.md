# Get Started
## About Bactopia
We initially released [Staphopia](https://staphopia.emory.edu/) in the summer
of 2018, for the analysis of *Staphylococcus aureus* genomes. At conferences,
or via email, the most common question about Staphopia that we received was:

*"How can I use Staphopia for my bacteria of interest?"*.

We too, have always been interested in this, because we also work with other
bacteria and would like a standard pipeline for incoming sequencing. Given the
self-interest and community interest, we developed Bactopia, *a generic pipeline
for the analysis of short-read bacterial sequences*! Side-note... we
kept the *"-opia"* in the name to pay homage to Staphopia!

While the philosophy behind Staphopia and Bactopia is mostly the same. The
development of Bactopia from scratch has really allowed us to take what we
learned from Staphopia and use it to our advantage. Bactopia has been developed
with usability, portability, and speed in mind.

Bactopia uses [Nextflow](https://www.nextflow.io/) to manage the workflow. We
have also prioritized software packages available from
[Bioconda](https://bioconda.github.io/) (or other
[Anaconda channels](https://anaconda.org/)) to make installation
as simple as possible for users. Finally, we have updated our workflow to
include the latest methods.

## Workflow Overview


## Installation
### Requirements
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


### Local Install
Once you have Conda all set up, you are ready to create an environment for
Bactopia.

#### Setting Up A Bactopia Conda Environment
The Conda envrionment configuration file
[bactopia.yml](https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml)
has been created to help ensure the minimum requirements to run Bactopia are
installed.

##### Download *bactopia.yml*
First you will want to download the `bactopia.yml` file.

```bash
wget https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml

- or -

curl -o bactopia.yml https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml
```

Or, feel free to "Right-click Save As" [**bactopia.yml**](https://raw.githubusercontent.com/staphopia/bactopia-ap/master/conda/bactopia.yml)

##### Create A Bactopia Environment

```bash
conda env create -n bactopia -f bactopia.yml
```

You may get the falling warning message:
```bash
Warning: you have pip-installed dependencies in your environment file, but you do not list pip itself as one of your conda dependencies.  Conda may not use the correct pip to install your packages, and they may end up in the wrong place.  Please add an explicit pip dependency.  I'm adding one for you, but still nagging you.
```
Please disregard that message. Pip is in the `bactopia.yml`, and the additional
Python packages are still installed correctly.


```bash
git clone git@github.com:staphopia/bactopia-ap.git
```



### Coming Soon
As the analysis pipeline for Bactopia is finalized there are plans to also
include [Nextflow Core](https://nf-co.re/) installations and containers via
[Docker](https://www.docker.com/) and [Singularity](https://www.sylabs.io/singularity/)
