FROM amazonlinux:2

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia count_31mers process"

RUN yum install -y gzip python3 && yum clean all && rm -rf /var/cache/yum

ENV PATH /mnt/conda/envs/bactopia-count_31mers/bin:/mnt/conda/envs/utils/bin:$PATH
