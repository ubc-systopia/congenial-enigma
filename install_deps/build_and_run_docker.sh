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


# python -m konect_scraper.cluster.main --mode download --directed --graph-numbers 20 30 --data-dir /home/atrostan/scratch/data

#    python -m konect_scraper.cluster.slurm.download /home/atrostan/workspace/repos/congenial-enigma/konect_scraper/cluster/csvs/download.csv 1


singularity shell --bind /home/atrostan/scratch/data,/scratch/atrostan/repos/congenial-enigma/konect_scraper:/congenial-enigma  /home/atrostan/singularity-images/congenial_enigma.sif