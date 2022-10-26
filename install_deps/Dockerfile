FROM intel/oneapi-basekit:devel-ubuntu20.04
# FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

ENV BASIC_DEPS="build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libgit2-dev curl libpq-dev python3-pip python3 wget apt-transport-https ca-certificates gnupg-agent software-properties-common"
ARG reqs_path
RUN echo reqs_path=${reqs_path}

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y ${BASIC_DEPS}  

# This variable specifices the system packages needed to be installed
# nvidia-cuda-toolkit needed for ParBatchRCM can be ignored for now
ENV BUILD_DEPS="build-essential cmake flex bison libgl1 pkg-config" \ 
    APP_DEPS="curl wget git apt-transport-https ca-certificates gnupg-agent software-properties-common sqlite3 libsqlite3-dev libgoogle-perftools-dev numactl libnuma-dev python3.10 python3.10-distutils ninja-build g++-11 python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev libpthread-stubs0-dev" 

RUN add-apt-repository ppa:deadsnakes/ppa \
    && add-apt-repository -y ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y ${BUILD_DEPS} ${APP_DEPS}     

# after install, update g++, gcc path
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 40 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 60 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 40 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 

# install intel-basekit manually 
# use wget to fetch the Intel repository public key
# add to your apt sources keyring so that archives signed with this key will be trusted.
# remove the public key
# RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
#     && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
#     && rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
#     && echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list \
#     && add-apt-repository "deb https://apt.repos.intel.com/oneapi all main" \
#     && apt install -y intel-basekit

# RUN ./opt/intel/oneapi/setvars.sh

# Install igraph from source
RUN git clone https://github.com/igraph/igraph.git \
    && cd igraph \
    && mkdir build && cd build \
    && cmake .. \
    && cmake --build . \
    && cmake --install .  

# Install boost 1.58.0 from source
ENV BOOST_VERSION_MAJOR=1
ENV BOOST_VERSION_MINOR=58
ENV BOOST_VERSION_REVISION=0

ENV BOOST_VERSION_FNAME="boost_${BOOST_VERSION_MAJOR}_${BOOST_VERSION_MINOR}_${BOOST_VERSION_REVISION}"
ENV BOOST_URL="http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_REVISION}/${BOOST_VERSION_FNAME}.tar.gz"

RUN wget ${BOOST_URL} \
    && tar -xzvf ${BOOST_VERSION_FNAME}.tar.gz \
    && cd ${BOOST_VERSION_FNAME} \
    && ./bootstrap.sh --prefix=/usr/local --with-libraries=atomic,serialization,graph,filesystem,system
RUN cd ${BOOST_VERSION_FNAME} \
    && cpuCores=`cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}'`\
    && echo "Available CPU cores: "$cpuCores \
    && ./b2 --with-program_options --with-atomic --with-serialization --with-graph --with-filesystem --with-system -j $cpuCores install 


# copy reqs txt from local to image
COPY requirements.txt requirements.txt
# install all python3.10 requirements
RUN which pip3
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10
RUN pip3.10 install -r requirements.txt



