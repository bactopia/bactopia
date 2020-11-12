FROM amazonlinux:2

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia (AWS Batch)"

RUN yum install -y git gzip python3 && \
    amazon-linux-extras install docker && \
    chkconfig docker on && \
    yum clean all && \
    rm -rf /var/cache/yum

ENV PATH /home/ec2-user/miniconda3/envs/bactopia/bin:/home/ec2-user/miniconda3/envs/utils/bin:/home/ec2-user/miniconda3/bin:$PATH
