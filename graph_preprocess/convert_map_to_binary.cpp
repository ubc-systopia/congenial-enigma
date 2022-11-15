#include <cstdint>
#include <csignal>
#include <string>
#include <vector>
#include "io.h"

//
// Created by atrostan on 14/11/22.
//
/**
 * All reorderings are originally saved text files
 * Given an input path to a vertex reordering, read it into a vector, and save it using a
 * boost::serialization binary archive.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint64_t num_vertices = 0;
	std::string input_path;
	std::string output_path;

	while ((opt = getopt(argc, argv, "n:i:o:")) != -1) {
		switch (opt) {
			case 'n':
				num_vertices = atol(optarg);
				break;
			case 'i':
				input_path = optarg;
				break;
			case 'o':
				output_path = optarg;
				break;
			default:
				abort();
		}
	}

	std::vector<uint32_t> orig_map(num_vertices);
	read_map(input_path, orig_map);
	write_binary_container<std::vector<uint32_t>>(output_path, orig_map);


}