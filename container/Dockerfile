#FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04
FROM rapidsai/miniforge-cuda:cuda12.2.2-base-ubuntu22.04-py3.11

RUN apt-get -y --allow-releaseinfo-change update && apt-get -y upgrade && apt-get install -y procps

RUN pip install --upgrade pip
COPY ./requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
