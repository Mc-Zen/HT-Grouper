#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "graph.h"
#include "formatting.h"


using namespace Q;


TEST_CASE("Graph") {
	auto g = Graph<3>::fullyConnectedGraph();
	REQUIRE(g.adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,1,1,1,0,1,1,1,0 });
	REQUIRE(Graph<3>::starGraph(1).adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,1,0,1,0,1,0,1,0 });
	REQUIRE(Graph<3>::starGraph(2).adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,0,1,0,0,1,1,1,0 });
}


TEST_CASE("Efficient graph") {
	auto g = Graph<3>::starGraph(2);
	efficient::Graph eg(g);
	println("{}", g.getAdjacencyMatrix());
	println("{}", eg.getAdjacencyMatrix());
	REQUIRE(eg.adjacencyMatrix.rows[0] == 4);
	REQUIRE(eg.adjacencyMatrix.rows[1] == 4);
	REQUIRE(eg.adjacencyMatrix.rows[2] == 3);

	REQUIRE(eg.toGraph().adjacencyMatrix == g.adjacencyMatrix);

}

TEST_CASE("Compress/Decompress") {
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::starGraph())) == Graph<9>::starGraph());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::starGraph(4))) == Graph<9>::starGraph(4));
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::starGraph(6))) == Graph<9>::starGraph(6));
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::pusteblumeGraph())) == Graph<9>::pusteblumeGraph());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::fullyConnectedGraph())) == Graph<9>::fullyConnectedGraph());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::emptyGraph())) == Graph<9>::emptyGraph());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::linearGraph())) == Graph<9>::linearGraph());
}