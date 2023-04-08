#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "graph.h"
#include "formatting.h"


using namespace Q;


TEST_CASE("Graph") {
	auto g = Graph<3>::fullyConnected();
	REQUIRE(g.adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,1,1,1,0,1,1,1,0 });
	REQUIRE(Graph<3>::star(1).adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,1,0,1,0,1,0,1,0 });
	REQUIRE(Graph<3>::star(2).adjacencyMatrix == Graph<3>::AdjacencyMatrix{ 0,0,1,0,0,1,1,1,0 });
}


TEST_CASE("Efficient graph") {
	auto g = Graph<3>::star(2);
	efficient::Graph eg(g);
	println("{}", g.getAdjacencyMatrix());
	println("{}", eg.getAdjacencyMatrix());
	REQUIRE(eg.adjacencyMatrix.rows[0] == 4);
	REQUIRE(eg.adjacencyMatrix.rows[1] == 4);
	REQUIRE(eg.adjacencyMatrix.rows[2] == 3);

	REQUIRE(eg.toGraph().adjacencyMatrix == g.adjacencyMatrix);

}

TEST_CASE("Compress/Decompress") {
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::star())) == Graph<9>::star());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::star(4))) == Graph<9>::star(4));
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::star(6))) == Graph<9>::star(6));
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::pusteblume())) == Graph<9>::pusteblume());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::fullyConnected())) == Graph<9>::fullyConnected());
	REQUIRE(Graph<9>::decompress(Graph<9>::compress(Graph<9>::linear())) == Graph<9>::linear());

	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::star(9))) == Graph<>::star(9));
	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::star(9, 4))) == Graph<>::star(9, 4));
	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::star(9, 6))) == Graph<>::star(9, 6));
	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::pusteblume(9))) == Graph<>::pusteblume(9));
	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::fullyConnected(9))) == Graph<>::fullyConnected(9));
	REQUIRE(Graph<>::decompress(9, Graph<>::compress(Graph<>::linear(9))) == Graph<>::linear(9));
}