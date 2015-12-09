#ifndef ORCA_GRAPH_HPP
#define ORCA_GRAPH_HPP 

#include <graph/Graph.hpp>
#include <vector>
#include <utility>

typedef typename boost::adjacency_list<
	boost::setS,
	boost::vecS,
	boost::undirectedS,
	graph::LabeledVertex,
	graph::LabeledEdge,
	graph::LabeledGraph
> Graph;

template<typename G, typename V>
inline void get_edges(const G &g, std::vector<std::pair<V,V>> &out) {
	out.clear();
	for(auto it = boost::edges(g); it.first != it.second; ++it.first) {
		V u = source(*it.first, g);
		V v = target(*it.first, g);

		out.push_back(std::pair<V,V>(u, v));
	}
}

template<typename G>
inline void remove_edge_loops(G &g) {
	typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
	remove_edge_if([&g](
		const edge_descriptor &e) {
			return source(e, g) == target(e, g);
		},
		g
	);
}

#endif
