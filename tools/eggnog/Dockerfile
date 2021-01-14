# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh eggnog 1

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool eggnog"
COPY --from=build /conda /opt/conda/envs/bactopia-eggnog
ENV PATH /opt/conda/envs/bactopia-eggnog/bin:$PATH
