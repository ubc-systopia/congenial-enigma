#!/bin/bash

# like the Dockerfile, this script will install all congenial enigma dependencies from an (assumed) bare
# metal ubuntu installation
# basically just replicates the RUN instructions in Dockerfile

BASIC_DEPS="build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libgit2-dev curl libpq-dev python3-pip python3 wget apt-transport-https ca-certificates gnupg-agent software-properties-common build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev"

sudo apt-get update \
    && sudo apt-get upgrade -y \
    && sudo apt-get install -y ${BASIC_DEPS}

# Install newer (>3.16) version of cmake (needed for igraph)
BUILD_DEPS="build-essential flex bison libgl1 pkg-config" 
APP_DEPS="curl wget apt-transport-https ca-certificates gnupg-agent software-properties-common sqlite3 libsqlite3-dev libgoogle-perftools-dev numactl intel-basekit" 

wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
    # add to your apt sources keyring so that archives signed with this key will be trusted.
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
    # remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo add-apt-repository ppa:deadsnakes/ppa

sudo apt-get update \
    && sudo apt-get upgrade -y \
    && sudo apt-get install -y ${BUILD_DEPS} ${APP_DEPS} \
    && sudo /opt/intel/oneapi/setvars.sh

# install cmake from source (cmake > 3.16 needed for igraph)
cmake_fname="cmake-3.25.0-rc1"
cmake_url="https://github.com/Kitware/CMake/releases/download/v3.25.0-rc1/cmake-3.25.0-rc1.tar.gz"
cd ~ 
wget ${cmake_url}
tar -xzvf ${cmake_fname}.tar.gz
cd ${cmake_fname}
./configure
make
sudo make install

# Install igraph from source
cd ~
git clone https://github.com/igraph/igraph.git 
cd igraph 
mkdir build 
cd build 
cmake .. 
cmake --build . 
cmake --install . 

# Install boost 1.74.0 from source
BOOST_VERSION_MAJOR=1
BOOST_VERSION_MINOR=74
BOOST_VERSION_REVISION=0
BOOST_VERSION_FNAME="boost_${BOOST_VERSION_MAJOR}_${BOOST_VERSION_MINOR}_${BOOST_VERSION_REVISION}"
BOOST_URL="http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_REVISION}/${BOOST_VERSION_FNAME}.tar.gz"
cd ~
wget ${BOOST_URL}
tar -xzvf ${BOOST_VERSION_FNAME}.tar.gz
cd ${BOOST_VERSION_FNAME}
# get the no of cpucores to make faster
cpuCores=`cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}'`
echo "Available CPU cores: "$cpuCores
./bootstrap.sh  # this will generate ./b2
sudo ./b2 --with=all -j $cpuCores install
echo "Verifying boost version: "
cat /usr/local/include/boost/version.hpp | grep "BOOST_LIB_VERSION"

# install python3.10 

# create virtual environment