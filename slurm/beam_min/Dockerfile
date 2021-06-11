FROM ubuntu:20.04

MAINTAINER Aman Chokshi <achokshi@student.unimelb.edu.au>

# Set timezone to utc
ENV TZ 'UTC'
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update and install packages
RUN apt-get -y update && \
    apt-get -y install git \
    gcc \
    wget \
    python3 \
    gfortran \
    binutils \
    python3-pip \
    openmpi-bin \
    openmpi-doc \
    libopenmpi-dev \
    openmpi-common \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get autoremove \
    && apt-get clean

ADD http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 /mwa_full_embedded_element_pattern.h5
ENV MWA_BEAM_FILE="/mwa_full_embedded_element_pattern.h5"

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
        emcee \
        healpy \
        mpi4py \
        schwimmbad \
        mwa_hyperbeam

ENTRYPOINT /bin/bash
