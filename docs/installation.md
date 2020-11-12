# Installation
Bactopia has a **a lot** of tools built into its workflow. As you can imagine, all these tools lead to numerous dependencies, and navigating dependencies can often turn into a very frustrating process. With this in mind, from the onset Bactopia was developed to only include programs that are installable using [Conda](https://conda.io/en/latest/).

Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. In other words, it makes it super easy to get the tools you need installed! The [official Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is a good starting point for getting started with Conda. Bactopia has been tested using the
[Miniconda installer](https://conda.io/en/latest/miniconda.html),
but the [Anaconda installer](https://www.anaconda.com/distribution/) should
work the same.

Containers (Docker and Singularity) are also available.

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

!!! error "OSX has limited support"
    I have developed Bactopia primarily for Linux, but I recognize it is useable on Mac OSX. Currently the support for OSX will be limited due to not having significant resources available for testing OSX extensively. My current setup, a mid-2013 MacBook, only allows me to maintain the Conda YAMLS for OSX. Please keep this in mind when using Bactopia on OSX. I will still try to help out if you run into any issues!

!!! error "Windows is not supported, please use Windows Subsystem for Linux"
    Bactopia will never support Windows natively due to dependencies. To use Bactopia on a Windows 10 machine, you will need to set up Windows Subsystem for Linux (WSL). This would allow you to run Bactopia inside the Linux subsystem. I have not tested Bactopia on WSL, but assume it should work fine. I have limited resources to test Bactopia in WSL, but if you give it a go and run into any issues please reach out!

## Container
A Docker and Singularity container has been created that is based off the Conda install.

```
# Docker 
docker pull bactopia/bactopia

# Singularity
singularity pull library://rpetit3/bactopia/bactopia
```

!!! error "These might not be available"
    Recent changes to DockerHub and Singularity policies have made the availability of these containers questionable. While I will try my best to make sure they are available, please understand if they are not. 
    
    If they are not available, I maintain the `Dockerfile` and `Singularity` files which can be used to build the containers locally.
