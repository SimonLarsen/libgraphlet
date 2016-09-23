#include <fstream>
#include <numeric>
#include <tclap/CmdLine.h>
#include <graph/GraphReader.hpp>
#include <orca/Orca.hpp>
#include <libgraphlet/GDD.hpp>
#include "Graph.hpp"

int main(int argc, const char **argv) {
	TCLAP::CmdLine cmd(
		"gdd_agreement",
		"Computes GDD-agreement between two networks.",
		"0.1",
		"Simon Larsen <simonhffh@gmail.com>"
	);

	TCLAP::ValueArg<int> graphletSizeArg("s", "size", "Graphlet size. 2-5 supported. Default: 4", false, 4, "size", cmd);
	TCLAP::UnlabeledValueArg<std::string> graph1Arg("graph1", "Path to first graph file", true, "", "GRAPH", cmd);
	TCLAP::UnlabeledValueArg<std::string> graph2Arg("graph2", "Path to second graph file", true, "", "GRAPH", cmd);
	TCLAP::UnlabeledValueArg<std::string> outputArg("output", "Output file", true, "", "FILE", cmd);

	cmd.parse(argc, argv);

	// Load graphs
	std::cerr << "Loading graphs (1/2)";
	Graph g1;
	graph::readGraph(graph1Arg.getValue(), g1);
	remove_edge_loops(g1);
	std::vector<std::pair<size_t,size_t>> edges1;
	get_edges(g1, edges1);

	std::cerr << "\rLoading graphs (2/2)" << std::endl;
	Graph g2;
	graph::readGraph(graph2Arg.getValue(), g2);
	remove_edge_loops(g2);
	std::vector<std::pair<size_t,size_t>> edges2;
	get_edges(g2, edges2);

	// Compute GDVs
	std::cerr << "Computing graphlet degree vectors (1/2)";
	orca::Orca orca1(num_vertices(g1), edges1, graphletSizeArg.getValue());
	orca1.compute();

	std::cerr << "\rComputing graphlet degree vectors (2/2)" << std::endl;
	orca::Orca orca2(num_vertices(g2), edges2, graphletSizeArg.getValue());
	orca2.compute();

	// Compute GDDs
	std::cerr << "Computing graphlet degree distributions (1/2)";
	libgraphlet::GDD gdd1;
	libgraphlet::gdd(orca1, gdd1, true);

	std::cerr << "\rComputing graphlet degree distributions (2/2)" << std::endl;
	libgraphlet::GDD gdd2;
	libgraphlet::gdd(orca2, gdd2, true);

	// Compute GGD-agreement
	std::cerr << "Computing GDD agreement" << std::endl;
	std::vector<float> gdda;
	libgraphlet::gdd_agreement(gdd1, gdd2, gdda);

	// Calculate overall agreement
	float sum = std::accumulate(gdda.begin(), gdda.end(), 0.0f);
	float mean = sum / (float)gdda.size();

	std::cerr << "Mean agreement: " << mean << std::endl;

	// Output to file
	std::ofstream file(outputArg.getValue());
	for(size_t i = 0; i < gdda.size(); ++i) {
		file << gdda[i] << "\n";
	}
	file.close();

	return 0;
}
