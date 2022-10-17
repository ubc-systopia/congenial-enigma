#!/bin/bash

# comparably to the Dockerfile, this script will install all congenial enigma dependencies from an (assumed) bare
# metal ubuntu installation
# basically just replicates the RUN instructions in Dockerfile


BASIC_DEPS="build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libgit2-dev curl libpq-dev python3-pip python3 wget apt-transport-https ca-certificates gnupg-agent software-properties-common"

sudo apt-get update \
    && sudo apt-get upgrade -y \
    && sudo apt-get install -y ${BASIC_DEPS}

BUILD_DEPS="build-essential cmake flex bison libgl1 pkg-config" 
APP_DEPS="curl wget apt-transport-https ca-certificates gnupg-agent software-properties-common sqlite3 libsqlite3-dev libgoogle-perftools-dev numactl intel-basekit" 

wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
    # add to your apt sources keyring so that archives signed with this key will be trusted.
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
    # remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB 
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"

sudo apt-get update \
    && sudo apt-get upgrade -y \
    && sudo apt-get install -y ${BUILD_DEPS} ${APP_DEPS} \
    && sudo /opt/intel/oneapi/setvars.sh