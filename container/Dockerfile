FROM nvidia/cuda:11.0.3-cudnn8-runtime-ubuntu20.04
#FROM nvidia/cuda:11.0.0-runtime-ubuntu20.04

#RUN conda install -c pytorch pytorch>=1.7.0 torchvision torchaudio #cudatoolkit=10.2

RUN apt-get -y --allow-releaseinfo-change update && apt-get -y upgrade && apt-get install -y procps python3-pip python3 python-is-python3

RUN pip install --upgrade pip
COPY ./requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
