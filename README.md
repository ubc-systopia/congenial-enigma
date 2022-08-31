# congenial-enigma
This repo contains scripts for pre-processing graphs including re-ordering, cleaning, ... etc

## pre-process-big-files.sh 

input: an edge list (input.txt)
output: three files 
		- input.txt.sorted.uniq: a list of uniq node ids


# Requirements
- [cmake](https://cmake.org/install/)
- ninja
	- ```sudo apt update && sudo apt install ninja-build```
- (If using an intellij ide (e.g. PyCharm, Clion), cmake, ninja should be included)
## Python
## C++
- [Boost](https://www.boost.org/)
  - `sudo apt-get install libboost-all-dev`
- [igraph](https://igraph.org/c/)
- [oneDPL](https://www.intel.com/content/www/us/en/developer/articles/guide/installation-guide-for-oneapi-toolkits.html)

# Examples

1. Download and compress `petster-hamster-household` and `opsahl-powergrid` graphs from konect, compute the `slashburn` order of those graphs, and plot the spy plots of the results.
```
python main.py \
	--initialize \
	--graph-names petster-hamster-household opsahl-powergrid \
	--download \
	--preproces sb \
	--plot
```
