# syntax=docker/dockerfile:1
FROM tensorflow/tensorflow:2.7.0
RUN apt update
RUN apt upgrade -y

RUN apt install mmseqs2 -y
RUN apt-get install libboost-numpy1.71 libboost-python1.71 -y

WORKDIR /metagenomic-deepfri

COPY setup.py setup.py
RUN pip install .

COPY post_setup.py post_setup.py
COPY CONFIG/FOLDER_STRUCTURE.py CONFIG/FOLDER_STRUCTURE.py
COPY utils/utils.py utils/utils.py
# download weights into docker and remove unpacked files to save on docker image size
RUN python post_setup.py && rm -rf /data

COPY . .
