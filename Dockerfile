# syntax=docker/dockerfile:1
FROM ubuntu
RUN apt update
RUN apt upgrade -y
RUN apt install python3 -y
RUN apt install python3-pip -y
RUN apt install mmseqs2 -y

WORKDIR /app
COPY . .

RUN pip install .
RUN apt-get install libboost-numpy1.71-dev libboost-python1.71-dev -y
RUN apt install wget -y