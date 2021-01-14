# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh bactopia

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.6"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia"
COPY --from=build /conda /opt/conda/envs/bactopia
ENV PATH /opt/conda/envs/bactopia/bin:$PATH
