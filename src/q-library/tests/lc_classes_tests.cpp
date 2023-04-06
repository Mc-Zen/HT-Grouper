#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "lc_classes.h"
#include "graph.h"


using namespace Q;



TEST_CASE("LC class 2") {
	//Graph<2> g;
	determine_lc_class<2>(getStabilizer(Graph<2>{}));
	determine_lc_class<3>(getStabilizer(Graph<3>{}));
	determine_lc_class<4>(getStabilizer(Graph<4>{}));

}