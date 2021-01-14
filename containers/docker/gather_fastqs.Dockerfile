# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh gather_fastqs

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia gather_fastqs process"
COPY --from=build /conda /opt/conda/envs/bactopia-gather_fastqs
ENV PATH /opt/conda/envs/bactopia-gather_fastqs/bin:$PATH
