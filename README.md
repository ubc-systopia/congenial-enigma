# congenial-enigma

This repo contains scripts for pre-processing graphs including re-ordering, cleaning, ... etc

## pre-process-big-files.sh

input: an edge list (input.txt)
output: three files

- input.txt.sorted.uniq: a list of uniq node ids

# Cloning

<`repo-name`> makes use of submodules to compute various vertex reordering and graph statistics. To clone <`repo-name`>:  
```
git clone https://github.com/ubc-systopia/congenial-enigma.git && cd congenial-enigma
git submodule update --force --recursive --init --remote
cd par_slashburn
git submodule update --force --recursive --init --remote 
```

# Requirements

## Docker
All the project's dependencies have been included in [`./install_deps/Dockerfile`](./install_deps/Dockerfile).  
You can create a singularity image by donwloading the existing container from Docker Hub to run the various executables on Compute Canada:
```
$ sudo singularity build congenial_enigma.sif docker://atrostan/install-deps:latest
```

This project's build relies on:

- [cmake](https://cmake.org/install/)
- ninja
    - ```sudo apt update && sudo apt install ninja-build```
- (If using an intellij ide (e.g. PyCharm, Clion), cmake, ninja should be included)

## Python

### Install Python3.10 virtualenv

> sudo apt install python3.10-venv

### Create a virtual env
Recommended to create a `virtualenv >= Python3.10`.  
`$ virtualenv --python="/usr/bin/python3.10" "./venv"`  
`$ pip install -r requirements.txt`

See [requirements.txt](./requirements.txt)

## C++

### g++11
```
add-apt-repository -y ppa:ubuntu-toolchain-r/test
apt-get update
apt-get install g++-11
```

- [Boost](https://www.boost.org/) (==1.5.8)
- [igraph](https://igraph.org/c/)
- [oneDPL](https://www.intel.com/content/www/us/en/developer/articles/guide/installation-guide-for-oneapi-toolkits.html)

### Rabbit Order

- g++ (4.9.2)
- Boost C++ library (1.58.0)
- libnuma (2.0.9)
- libtcmalloc_minimal in google-perftools (2.1)

# Build and Install
`$ <venv> python setup.py`
## Clean
`$ <venv> python setup.py --clean`

# Execution
Supply `driver.py` with a configuration file and an optional cluster mode (`--slurm`).  
See [`local_config.yaml`](./local_config.yaml) and [`cluster_config.yaml`](./konect_scraper/cluster/cluster_config.yaml) for examples of configuration files for  local and cluster executions.  
e.g.
```
$ python driver.py --config-path <path-to-config-file> --slurm
```
