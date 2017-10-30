FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y libcurl4-openssl-dev git make gcc \
					   zlib1g-dev wget unzip python2.7 \
					   python-dev python-pip bats awscli curl

# Import the repo to /usr/bin/paladin
RUN mkdir /usr/bin/paladin
ADD . /usr/bin/paladin/
RUN cd /usr/bin/paladin && \
    make && \
    export PATH=$PATH:/usr/bin/paladin

# Test the installation
RUN cd /usr/bin/paladin/sample_data && \
	bash make_test.sh && \
	cd /usr/bin/paladin && \
	rm -r sample_data
