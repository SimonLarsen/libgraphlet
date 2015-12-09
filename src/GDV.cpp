#include <fstream>
#include <tclap/CmdLine.h>
#include <graph/GraphReader.hpp>
#include <orca/Orca.hpp>
#include "Graph.hpp"

int main(int argc, const char **argv) {
	TCLAP::CmdLine cmd(
		"gdv",
		"Count graphlet degree vectors.",
		"0.1",
		"Simon Larsen <simonhffh@gmail.com>"
	);

	TCLAP::ValueArg<int> graphletSizeArg("s", "size", "Graphlet size. 2-5 supported. Default: 4", false, 4, "size", cmd);
	TCLAP::UnlabeledValueArg<std::string> graphArg("graph", "Path to graph file", true, "", "GRAPH", cmd);
	TCLAP::UnlabeledValueArg<std::string> outputArg("output", "Output file", true, "", "FILE", cmd);

	cmd.parse(argc, argv);

	// Read graph
	Graph g;
	graph::readGraph(graphArg.getValue(), g);
	remove_edge_loops(g);
	std::vector<std::pair<size_t,size_t>> edges;
	get_edges(g, edges);

	// Compute GDVs
	orca::Orca orca(num_vertices(g), edges, graphletSizeArg.getValue());
	orca.compute();

	// Write to file
	const auto &orbits = orca.getOrbits();
	std::ofstream file(outputArg.getValue());
	for(auto it1 = orbits.begin1(); it1 != orbits.end1(); ++it1) {
		for(auto it = it1.begin(); it != it1.end(); ++it) {
			file << *it << " ";
		}
		file << "\n";
	}
	file.close();

	return 0;
}
