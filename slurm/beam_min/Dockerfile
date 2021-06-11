FROM ubuntu:20.04

MAINTAINER Aman Chokshi <achokshi@student.unimelb.edu.au>

RUN apt-get update -y && \
        DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata

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

# Enable mpi as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1 
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

ADD http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 /mwa_full_embedded_element_pattern.h5
ENV MWA_BEAM_FILE="/mwa_full_embedded_element_pattern.h5"

# Make "python" available.
RUN ln -s /usr/bin/python3 /usr/bin/python

# Upgrade pip
RUN pip3 install \
        --upgrade \
        --no-cache-dir \
        pip

# Python requirements
RUN pip3 install \
        --no-cache-dir \
        h5py \
        tqdm \
        numpy \
        emcee \
        healpy \
        mpi4py \
        schwimmbad \
        mwa_hyperbeam

ENTRYPOINT /bin/bash
