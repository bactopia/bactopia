FROM amazonlinux:2

LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia to be used with AWS Batch"

ENV PATH /opt/conda/envs/bactopia/bin:$PATH
