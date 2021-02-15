FROM ubuntu:20.04

MAINTAINER Aman Chokshi <achokshi@student.unimelb.edu.au>

# Set timezone to utc
ENV TZ 'UTC'
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update and install packages
RUN apt-get -y update && \
    apt-get -y install git \
                          wget \
                          python3 \
                          python3-pip \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get autoremove \
    && apt-get clean

# Make "python" available.
RUN ln -s /usr/bin/python3 /usr/bin/python

# Upgrade pip
RUN pip3 install \
        --upgrade \
        --no-cache-dir \
        pip

# RMextract requirements
RUN pip3 install \
        --no-cache-dir \
        numpy \
        scipy \
        pandas \
        pyephem \
        astropy \
        Pillow

# Install RMextract
RUN git clone https://github.com/lofar-astron/RMextract.git && \
        cd RMextract && \
        python setup.py build && \
        python setup.py install && \
        cd .. && rm -rf RMextract

ENTRYPOINT /bin/bash
