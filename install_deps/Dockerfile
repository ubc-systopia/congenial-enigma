FROM intel/oneapi-basekit:devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

ENV BASIC_DEPS="build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libgit2-dev curl libpq-dev python3-pip python3 wget apt-transport-https ca-certificates gnupg-agent software-properties-common"

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y ${BASIC_DEPS}  

# This variable specifices the system packages needed to be installed
# nvidia-cuda-toolkit needed for ParBatchRCM can be ignored for now
ENV BUILD_DEPS="build-essential cmake flex bison libgl1 pkg-config" \ 
    APP_DEPS="curl wget git apt-transport-https ca-certificates gnupg-agent software-properties-common sqlite3 libsqlite3-dev libgoogle-perftools-dev numactl" 

RUN apt-get update \
    &&  apt-get upgrade -y \
    && apt-get install -y ${BUILD_DEPS} ${APP_DEPS} 

# Install igraph from source
RUN git clone https://github.com/igraph/igraph.git \
    && cd igraph \
    && mkdir build && cd build \
    && cmake .. \
    && cmake --build . \
    && cmake --install . 

# Install boost 1.74.0 from source
ENV BOOST_VERSION="1_74_0"
ENV BOOST_URL="http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION}/boost_${BOOST_VERSION}.tar.gz"
ENTRYPOINT echo $BOOST_URL
# install python3.10


