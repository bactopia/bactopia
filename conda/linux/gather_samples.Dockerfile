FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - gather_samples"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/gather_samples.yml"
LABEL conda.md5="e2ca67e0aae408581093c074882fc05e"

# Bactopia environment
COPY conda/linux/gather_samples.yml /
RUN mamba env create -q -f gather_samples.yml && \
    mamba clean -a -y
RUN chmod -R 777 /root && \
    mkdir -p /root/.ncbi/ && \
    mkdir /.ncbi && \
    printf '/LIBS/GUID = "%s"\n' $(uuidgen) > /root/.ncbi/user-settings.mkfg && \
    cp /root/.ncbi/user-settings.mkfg /.ncbi/user-settings.mkfg

# Add bactopia env to path
ENV PATH /opt/conda/envs/bactopia-gather_samples/bin:$PATH
