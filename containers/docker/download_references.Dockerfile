# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh download_references

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia download_references process"
COPY --from=build /conda /opt/conda/envs/bactopia-download_references
ENV PATH /opt/conda/envs/bactopia-download_references/bin:$PATH
