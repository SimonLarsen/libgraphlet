#include <fstream>
#include <tclap/CmdLine.h>
#include <graph/GraphReader.hpp>
#include <graph/Algorithms.hpp>
#include <orca/Orca.hpp>
#include <libgraphlet/GDD.hpp>
#include "Graph.hpp"

int main(int argc, const char **argv) {
	TCLAP::CmdLine cmd(
		"gdd",
		"Compute graphlet degree distribution (GDD).",
		"0.1",
		"Simon Larsen <simonhffh@gmail.com>"
	);

	TCLAP::ValueArg<int> graphletSizeArg("s", "size", "Graphlet size. 2-5 supported. Default: 4", false, 4, "size", cmd);
	TCLAP::UnlabeledValueArg<std::string> graphArg("graph", "Path to graph file", true, "", "GRAPH", cmd);
	TCLAP::UnlabeledValueArg<std::string> outputArg("output", "Output file", true, "", "FILE", cmd);
	TCLAP::SwitchArg normalizeSwitch("n", "normalize", "Normalize distribution", cmd, false);

	cmd.parse(argc, argv);

	// Read graph
	std::cerr << "Loading graph" << std::endl;
	Graph g;
	graph::readGraph(graphArg.getValue(), g);
	graph::removeEdgeLoops(g);
	std::vector<std::pair<size_t,size_t>> edges;
	graph::get_edges(g, edges);

	// Compute GDVs
	std::cerr << "Computing graphlet degree vectors" << std::endl;
	orca::Orca orca(num_vertices(g), edges, graphletSizeArg.getValue());
	orca.compute();

	// Compute GDD
	std::cerr << "Computing graphlet degree distribution" << std::endl;
	libgraphlet::GDD gdd;
	libgraphlet::gdd(orca, gdd, normalizeSwitch.getValue());

	// Write GDD to file
	std::ofstream file(outputArg.getValue());
	for(auto &v : gdd) {
		size_t max_key = (v.size() > 0 ? v.rbegin()->first : 0);

		for(size_t i = 0; i <= max_key; ++i) {
			if(v.find(i) != v.end()) {
				file << v[i];
			} else {
				file << "0";
			}
			if(i < max_key) file << " ";
		}
		file << "\n";
	}
	file.close();

	return 0;
}
