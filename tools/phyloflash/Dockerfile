# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh phyloflash 1
RUN conda create -n python3 -c conda-forge python=3 && \
    conda-pack -n python3 -o /tmp/python3.tar && \
    mkdir /python3 && \
    cd /python3 && tar xf /tmp/python3.tar && \
    rm /tmp/python3.tar && \
    /python3/bin/conda-unpack

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool phyloflash"
COPY --from=build /conda /opt/conda/envs/bactopia-phyloflash
COPY --from=build /python3 /opt/conda/envs/python3
ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:/opt/conda/envs/python3:$PATH
