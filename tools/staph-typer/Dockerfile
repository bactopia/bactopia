FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - staph-typer"
LABEL software.version="1.7.1"
LABEL description="A workflow with set of tools specific to typing Staphylococcus aureus (e.g. agr, SCCmec, spa)."
LABEL website="https://bactopia.github.io/bactopia-tools/staph-typer/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/staph-typer/environment-linux.yml"
LABEL conda.md5="378cab315aa371a2a456d9336e228058"

COPY tools/staph-typer/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && \
    conda clean -y -a && \
    printf ">none\nAAAA\n" > ./test.fna && \
    /opt/conda/envs/bactopia-staph-typer/bin/spaTyper -f ./test.fna && \
    /opt/conda/envs/bactopia-staph-typer/bin/spaTyper -f ./test.fna && \
    rm ./test.fna
ENV PATH /opt/conda/envs/bactopia-staph-typer/bin:$PATH
