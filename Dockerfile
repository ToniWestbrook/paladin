FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y libcurl4-openssl-dev git make gcc \
					   zlib1g-dev wget unzip python2.7 \
					   python-dev python-pip bats awscli curl && \
	pip install boto3==1.4.7

# Set the default langage to C
ENV LC_ALL C

# Import the repo to /usr/bin/paladin
RUN mkdir /usr/bin/paladin
ADD . /usr/bin/paladin/
RUN cd /usr/bin/paladin && \
    make && \
    ln -s /usr/bin/paladin/paladin /usr/local/bin/ && \
    ln -s /usr/bin/paladin/run.py /usr/local/bin/

# Test the installation
RUN cd /usr/bin/paladin/sample_data && \
	bash make_test.sh && \
	cd /usr/bin/paladin && \
	rm -r sample_data
