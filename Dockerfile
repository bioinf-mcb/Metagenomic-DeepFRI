# syntax=docker/dockerfile:1
FROM ubuntu:latest
RUN apt update &&  \
    apt upgrade -y &&  \
    apt-get update &&  \
    apt-get upgrade -y && \
    apt install python3-pip -y && \
    apt install mmseqs2 -y && \
    apt-get install libboost-numpy1.71 libboost-python1.71 -y && \
    apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

WORKDIR /metagenomic-deepfri

COPY requirements.txt requirements.txt
COPY setup.py setup.py
RUN pip3 install --no-cache-dir .

COPY . .
