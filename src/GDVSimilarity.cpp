#include <fstream>
#include <tclap/CmdLine.h>
#include <boost/format.hpp>
#include <graph/GraphReader.hpp>
#include <graph/Algorithms.hpp>
#include <orca/Orca.hpp>
#include <orca/Similarity.hpp>
#include "Graph.hpp"

int main(int argc, const char **argv) {
	TCLAP::CmdLine cmd(
		"gdv_similarity",
		"Compute GDV similarity matrix of two networks.",
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
	graph::removeEdgeLoops(g1);
	std::vector<std::pair<size_t,size_t>> edges1;
	graph::get_edges(g1, edges1);

	std::cerr << "\rLoading graphs (2/2)" << std::endl;
	Graph g2;
	graph::readGraph(graph2Arg.getValue(), g2);
	graph::removeEdgeLoops(g2);
	std::vector<std::pair<size_t,size_t>> edges2;
	graph::get_edges(g2, edges2);

	// Compute GDVs
	std::cerr << "Computing graphlet degree vectors (1/2)";
	orca::Orca orca1(num_vertices(g1), edges1, graphletSizeArg.getValue());
	orca1.compute();

	std::cerr << "\rComputing graphlet degree vectors (2/2)" << std::endl;
	orca::Orca orca2(num_vertices(g2), edges2, graphletSizeArg.getValue());
	orca2.compute();

	// Compute similarity matrix
	std::cerr << "Computing similarity matrix" << std::endl;
	boost::numeric::ublas::matrix<float> sim;
	orca::similarity(orca1, orca2, sim);

	// Write to file
	std::ofstream file(outputArg.getValue());
	for(size_t i = 0; i < num_vertices(g1); ++i) {
		for(size_t j = 0; j < num_vertices(g2); ++j) {
			file << boost::format("%s\t%s\t%f\n") % g1[i].label % g2[j].label % sim(i, j);
		}
	}
	file.close();

	std::cerr << "Done!" << std::endl;

	return 0;
}
