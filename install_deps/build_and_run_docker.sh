#!/bin/bash

docker build \
    -t install-deps \
    .
    # --build-arg arg=2.3 \
    # --build-arg pip_reqs=${}
    # .
docker run install-deps
# sudo apt-get -y install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev

# sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev

# docker run -it \
#    install-deps:latest /bin/bash \
#    --mount type=bind,source=/home/atrostan/Workspace/repos/congenial-enigma/,target=/app
#    /bin/bash


git clone https://github.com/ubc-systopia/congenial-enigma.git && cd congenial-enigma
git submodule update --force --recursive --init --remote
cd par_slashburn
git submodule update --force --recursive --init --remote 
cd ../


docker run -it \
   --mount type=bind,source=/home/atrostan/workspace/repos/tmp/congenial-enigma,target=/congenial-enigma \
   install-deps:latest /bin/bash 

   