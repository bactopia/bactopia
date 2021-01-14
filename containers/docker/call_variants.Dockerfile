# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh call_variants

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia call_variants process"
COPY --from=build /conda /opt/conda/envs/bactopia-call_variants
RUN apt-get update && apt-get install -y locales && apt-get clean -y && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.UTF-8
ENV PATH /opt/conda/envs/bactopia-call_variants/bin:$PATH
