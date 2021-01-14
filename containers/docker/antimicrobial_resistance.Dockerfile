# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh antimicrobial_resistance && \
    export CONDA_PREFIX=/conda && \
    /conda/bin/amrfinder -u && \
    rm /conda/share/amrfinderplus/data/latest && \
    ls /conda/share/amrfinderplus/data/ | xargs -I {} mv /conda/share/amrfinderplus/data/{} /conda/share/amrfinderplus/data/docker

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia antimicrobial_resistance process"
COPY --from=build /conda /opt/conda/envs/bactopia-antimicrobial_resistance
ENV PATH /opt/conda/envs/bactopia-antimicrobial_resistance/bin:$PATH
ENV CONDA_PREFIX /opt/conda/envs/bactopia-antimicrobial_resistance
ENV AMRFINDER_DB /opt/conda/envs/bactopia-antimicrobial_resistance/share/amrfinderplus/data/docker
