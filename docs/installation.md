# Installation
Bactopia has a **a lot** of tools built into its workflow. As you can imagine, all these tools lead to numerous dependencies, and navigating dependencies can often turn into a very frustrating process. With this in mind, from the onset Bactopia was developed to only include programs that are installable using [Conda](https://conda.io/en/latest/).

Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. In other words, it makes it super easy to get the tools you need installed! The [official Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is a good starting point for getting started with Conda. Bactopia has been tested using the
[Miniconda installer](https://conda.io/en/latest/miniconda.html),
but the [Anaconda installer](https://www.anaconda.com/distribution/) should
work the same.

A Docker container based on the Bioconda install is also available.

## Bioconda
Once you have Conda all set up, you are ready to create an environment for
Bactopia. To do so, you can use the following command:

```
conda create -n bactopia -c conda-forge -c bioconda bactopia
```

After a few minutes you will have a new conda environment suitably named *bactopia*. To activate this environment, you will can use the following command:

```
conda activate bactopia
```

And voilà, you are all set to get started processing your data!


But first, it is highly recommended that you take the time to [Build Datasets](datasets.md) that Bactopia can take advantage of.

## Docker
A Docker container, that is also Singularity compatible, has been created that is based off the Conda install.

```
docker pull bactopia/bactopia
```
