# Build then use conda-pack
FROM continuumio/miniconda3:4.9.2 AS build
COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh TOOL_NAME 1

# Use the conda-pack version
FROM debian:buster AS runtime
LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool TOOL_NAME"
COPY --from=build /conda /opt/conda/envs/bactopia-TOOL_NAME
ENV PATH /opt/conda/envs/bactopia-TOOL_NAME/bin:$PATH
