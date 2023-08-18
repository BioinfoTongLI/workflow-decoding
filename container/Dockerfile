FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04

RUN apt-get -y --allow-releaseinfo-change update && apt-get -y upgrade && apt-get install -y procps python3-pip python3 python-is-python3

RUN pip install --upgrade pip
COPY ./requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
