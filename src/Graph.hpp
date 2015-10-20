#ifndef ORCA_GRAPH_HPP
#define ORCA_GRAPH_HPP 

#include <graph/Graph.hpp>

typedef typename boost::adjacency_list<
	boost::setS,
	boost::vecS,
	boost::undirectedS,
	graph::LabeledVertex,
	graph::LabeledEdge,
	graph::LabeledGraph
> Graph;

#endif
