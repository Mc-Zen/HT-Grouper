#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include <iostream>
#include "sector_length_distribution.h"


using namespace Q;

TEST_CASE("Sector length distribution") {
	auto graph = Graph<4>::starGraph(0);
	graph.removeEdge(0,2);
	graph.removeEdge(0,3);
	graph.addEdge(1,2);
	graph.addEdge(2,3);
	auto sld = sectorLengthDistribution(graph);
	for (int i : sld) {
		std::cout << i << " ";
	}
}