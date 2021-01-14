# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh assembly_qc

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia assembly_qc process"
COPY --from=build /conda /opt/conda/envs/bactopia-assembly_qc
ENV PATH /opt/conda/envs/bactopia-assembly_qc/bin:$PATH
