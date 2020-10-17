FROM amazonlinux:2

LABEL version="1.4.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia (AWS Batch)"

RUN yum install -y gzip python3 && yum clean all && rm -rf /var/cache/yum

ENV PATH /home/ec2-user/miniconda3/envs/bactopia/bin:/home/ec2-user/miniconda3/envs/utils/bin:/home/ec2-user/miniconda3/bin:$PATH
