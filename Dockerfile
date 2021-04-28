FROM ubuntu:16.04
MAINTAINER xufeng <crazyhsu9527@gmail.com>

# install basic dependencies
RUN sed -i s@/archive.ubuntu.com/@/mirrors.tuna.tsinghua.edu.cn/@g /etc/apt/sources.list
RUN apt-get clean && apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    apt-transport-https \
    libz-dev \
    ca-certificates \
    vim \
    less \
    sudo \
    apt-utils

RUN wget --no-check-certificate https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /root/miniconda3.sh && \
    /bin/bash /root/miniconda3.sh -b -p /root/conda && \
    rm /root/miniconda3.sh

ENV PATH=/root/conda/bin:$PATH

RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/

RUN conda update --all -y
RUN conda create -n modas python=3.7.8 -y

ENV CONDA_DEFAULT_ENV=modas
ENV CONDA_PREFIX=/root/conda/envs/$CONDA_DEFAULT_ENV
ENV CONDA_AUTO_UPDATE_CONDA=false

RUN echo ". /root/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate modas" >> ~/.bashrc

RUN git clone https://github.com/liusy-jz/MODAS.git /root/MODAS && \
    cd /root/MODAS && \
    $CONDA_PREFIX/bin/python setup.py build && $CONDA_PREFIX/bin/python setup.py install
