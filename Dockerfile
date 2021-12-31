# syntax=docker/dockerfile:1
FROM tensorflow/tensorflow:2.7.0
RUN apt update
RUN apt upgrade -y

RUN apt install mmseqs2 -y
RUN apt-get install libboost-numpy1.71 libboost-python1.71 -y

WORKDIR /app
COPY . .

RUN pip install .
