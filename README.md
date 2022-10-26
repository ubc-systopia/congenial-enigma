# congenial-enigma

This repo contains scripts for pre-processing graphs including re-ordering, cleaning, ... etc

## pre-process-big-files.sh

input: an edge list (input.txt)
output: three files

- input.txt.sorted.uniq: a list of uniq node ids

# Cloning

<`repo-name`> makes use of submodules to compute various vertex reordering. To clone <`repo-name`>:  
`$ git clone --recurse-submodules git@github.com:ubc-systopia/congenial-enigma.git`

# Requirements

- [cmake](https://cmake.org/install/)
- ninja
    - ```sudo apt update && sudo apt install ninja-build```
- (If using an intellij ide (e.g. PyCharm, Clion), cmake, ninja should be included)

## Python

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

# Build

`$ cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S ./graph_preprocess -B ./graph_preprocess/<cmake build directory name>`

# Install

- `graph_preprocess`:  
  `$ cmake --build ./graph_preprocess/<cmake build directory name> --target graph_preprocess -j <num threads>`
- `slashburn`:  
  `$ cmake --build ./graph_preprocess/<cmake build directory name> --target slashshburn -j <num threads>`

# Examples

1. Download and compress `petster-hamster-household` and `opsahl-powergrid` graphs from konect, compute the `slashburn`
   order of those graphs, and plot the spy plots of the results.

```
python main.py \
	--initialize \
	--graph-names petster-hamster-household opsahl-powergrid \
	--download \
	--reorder sb \
	--plot
```

2. Compute the `slashburn` order of `petster-hamster-household` and `opsahl-powergrid` graphs and plot the spy plots of
   the results.

```
python main.py \
	--no-initialize \
	--graph-names petster-hamster-household opsahl-powergrid \
	--no-download \
	--reorder sb \
	--plot
```

# SQL

## All Possible Network Categories

- Affiliation network
- Animal network
- Authorship network
- Citation network
- Co-authorship network
- Co-citation network
- Communication network
- Computer network
- Feature network
- Human contact network
- Human social network
- Hyperlink network
- Infrastructure network
- Interaction network
- Lexical network
- Metabolic network
- Miscellaneous network
- Neural network
- Online contact network
- Online social network
- Rating network
- Software network
- Text network
- Trophic network


