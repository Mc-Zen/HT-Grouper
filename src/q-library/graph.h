
#pragma once
#include "efficient_binary_math.h"
#include <iostream>
#include <vector>

namespace Q {




	template<bool empty = true>
	struct GraphSize {
	};


	template<>
	struct GraphSize<false> {
		int n{};
	};



	template<size_t n = Math::dynamic>
	class Graph {
	public:
		using AdjacencyMatrix = Math::Matrix<Binary, n, n>;
		static constexpr bool is_dynamic = (n == Math::dynamic);
		GraphSize<!is_dynamic> graphSize;


		AdjacencyMatrix adjacencyMatrix;

		const AdjacencyMatrix& getAdjacencyMatrix() const { return adjacencyMatrix; }

		constexpr Graph() requires !is_dynamic = default;
		constexpr Graph(GraphSize<!is_dynamic> graphSize) : graphSize(graphSize) {
			if constexpr (is_dynamic) {
				adjacencyMatrix = AdjacencyMatrix(numVertices(), numVertices());
			}
		}
		explicit constexpr Graph(int numVertices) requires is_dynamic : graphSize{ numVertices }, adjacencyMatrix(numVertices, numVertices) {}


		constexpr int numVertices() const {
			if constexpr (is_dynamic) return graphSize.n;
			else return n;
		}


		constexpr static auto emptyGraph() {
			return Graph{};
		}

		constexpr static auto fullyConnectedGraph() requires !is_dynamic{
			Graph graph;
			graph.adjacencyMatrix.fill(1);
			for (int i = 0; i < numVertices(); ++i) graph.adjacencyMatrix(i, i) = 0;
			return graph;
		}
			constexpr static auto fullyConnectedGraph(int n) requires is_dynamic{
			Graph graph(n);
			graph.adjacencyMatrix.fill(1);
			for (int i = 0; i < graph.numVertices(); ++i) graph.adjacencyMatrix(i, i) = 0;
			return graph;
		}

			constexpr static auto starGraph(int center = 0) {
			Graph graph;
			auto& mat = graph.adjacencyMatrix;
			std::for_each(mat.col_begin(center), mat.col_end(center), [](Binary& el) { el.negate(); });
			std::for_each(mat.row_begin(center), mat.row_end(center), [](Binary& el) { el.negate(); });
			return graph;
		}

		constexpr static auto linearGraph() {
			Graph graph;
			for (int i = 0; i < n - 1; ++i) {
				graph.addEdge(i, i + 1);
			}
			return graph;
		}
		

		constexpr static auto linearGraph(int numVertices) requires is_dynamic {
			Graph graph(numVertices);
			for (int i = 0; i < numVertices - 1; ++i) {
				graph.addEdge(i, i + 1);
			}
			return graph;
		}

		constexpr static auto cycleGraph() {
			auto graph = linearGraph();
			graph.addEdge(0, n - 1);
			return graph;
		}

		constexpr static auto pusteblumeGraph() {
			static_assert(n >= 5, "The Pusteblume graph is only possible for at least 5 vertices");
			Graph graph;
			auto& mat = graph.adjacencyMatrix;
			for (size_t i = 1; i < 4; ++i) {
				graph.addEdge(0, i);
			}
			for (size_t i = 4; i < n; ++i) {
				graph.addEdge(3, i);
			}
			return graph;
		}

		constexpr bool hasEdge(int vertex1, int vertex2) const {
			return adjacencyMatrix(vertex1, vertex2) == 1;
		}

		constexpr int edgeCount() const {
			return std::accumulate(adjacencyMatrix.begin(), adjacencyMatrix.end(), 0, [](int count, Binary b) { return count + b.toInt(); }) / 2;
		}

		constexpr void addEdge(int vertex1, int vertex2) {
			if (vertex1 == vertex2) return;
			adjacencyMatrix(vertex1, vertex2) = 1;
			adjacencyMatrix(vertex2, vertex1) = 1;
		}

		constexpr void addPath(const std::initializer_list<int>& vertices) {
			if (vertices.size() < 2) return;
			int previousVertex = *vertices.begin();
			for (auto it = vertices.begin() + 1; it != vertices.end(); ++it) {
				addEdge(*it, previousVertex);
				previousVertex = *it;
			}
		}

		constexpr void removeEdge(int vertex1, int vertex2) {
			adjacencyMatrix(vertex1, vertex2) = 0;
			adjacencyMatrix(vertex2, vertex1) = 0;
		}

		constexpr void removeEdgesTo(int vertex) {
			for (int i = 0; i < numVertices(); ++i) removeEdge(i, vertex);
		}

		constexpr void toggleEdge(int vertex1, int vertex2) {
			adjacencyMatrix(vertex1, vertex2).negate();
			adjacencyMatrix(vertex2, vertex1).negate();
		}

		constexpr void localComplementation(int vertex) {
			Math::Vector<Binary, numVertices()> ithColumn = adjacencyMatrix.col(vertex);
			adjacencyMatrix += ithColumn * ithColumn.transpose();
			for (int i = 0; i < numVertices(); ++i) {
				adjacencyMatrix(i, i) = 0;
			}
		}

		constexpr void swap(int vertex1, int vertex2) {
			Math::Vector<Binary, numVertices()> col = adjacencyMatrix.col(vertex1);
			adjacencyMatrix.col(vertex1) = adjacencyMatrix.col(vertex2);
			adjacencyMatrix.col(vertex2) = col;
			Math::RowVector<Binary, numVertices()> row = adjacencyMatrix.row(vertex1);
			adjacencyMatrix.row(vertex1) = adjacencyMatrix.row(vertex2);
			adjacencyMatrix.row(vertex2) = row;
		}

		/// @brief Perform a series of local complementations
		/// @param vertices Local complementations will be executed for vertices in the given order
		constexpr void localComplementation(const std::initializer_list<int>& vertices) {
			std::ranges::for_each(vertices, [this](int vertex) { localComplementation(vertex); });
		}


		/// @brief Perform a vertex permutation from {0,1,2,3,...} to mapping. 
		/// @param mapping Each number from 0 to numVertices-1 needs to occur exactly once. 
		/// @return permuted graph
		constexpr Graph graphIsomorphism(const std::vector<int>& mapping) const {
			Graph result;
			for (int i = 0; i < numVertices() - 1; ++i) {
				for (int j = i + 1; j < numVertices(); ++j) {
					result.adjacencyMatrix(mapping[i], mapping[j]) = adjacencyMatrix(i, j);
				}
			}
			return result;
		}

		constexpr void clear() {
			adjacencyMatrix.fill(0);
		}


		template<class Predicate>
		static constexpr Graph transform(const Graph& g1, const Graph& g2, Predicate predicate) {
			Graph result;
			std::transform(g1.adjacencyMatrix.begin(), g1.adjacencyMatrix.end(), g2.adjacencyMatrix.begin(),
				result.adjacencyMatrix.begin(), predicate);
			return result;
		}

		template<class Predicate>
		constexpr void transform(const Graph& g, Predicate predicate) {
			Graph result;
			std::transform(adjacencyMatrix.begin(), adjacencyMatrix.end(), g.adjacencyMatrix.begin(), adjacencyMatrix.begin(), predicate);
		}

		static constexpr Graph add(const Graph& g1, const Graph& g2) {
			return Graph::transform(g1, g2, [](Binary b1, Binary b2) { return b1 | b2; });
		}

		static constexpr Graph intersect(const Graph& g1, const Graph& g2) {
			return Graph::transform(g1, g2, [](Binary b1, Binary b2) { return b1 & b2; });
		}

		/// @brief Subtract edges of g2 from g1
		static constexpr Graph subtract(const Graph& g1, const Graph& g2) {
			return Graph::transform(g1, g2, [](Binary b1, Binary b2) { return b1 * b2.negate(); });
		}

		/// @brief Add edges from other graph to this graph
		constexpr void add(const Graph& g) {
			transform(g, [](Binary b1, Binary b2) { return b1 | b2; });
		}

		/// @brief Form intersection of this graphs and the other graphs edges
		constexpr void intersect(const Graph& g) {
			transform(g, [](Binary b1, Binary b2) { return b1 & b2; });
		}

		/// @brief Remove all edges of this graph that occur in the other graph
		constexpr void subtract(const Graph& g) {
			transform(g, [](Binary b1, Binary b2) { return b1 * b2.negate(); });
		}

		/// @brief Get all edges in form of integer pairs
		auto getEdges() const {
			std::vector<std::pair<int, int>> edges;
			for (int i = 0; i < numVertices() - 1; ++i) {
				for (int j = i + 1; j < numVertices(); ++j) {
					if (hasEdge(i, j)) edges.emplace_back(i, j);
				}
			}
			return edges;
		}

		constexpr friend bool operator==(const Graph& g1, const Graph& g2) = default;

		/// @brief Compress the graph into a single 64-bit integer. 
		static int64_t compress(const Graph& graph) {
			static_assert(numVertices() * (numVertices() - 1) / 2 <= 64);
			int64_t code{};
			int index{};
			for (int i = 0; i < numVertices() - 1; ++i) {
				for (int j = i + 1; j < numVertices(); ++j) {
					if (graph.hasEdge(i, j)) code |= (1ULL << index);
					++index;
				}
			}
			return code;
		}

		/// @brief Restore a graph from its compressed form. 
		static Graph decompress(int64_t code) {
			static_assert(numVertices() * (numVertices() - 1) / 2 <= 64);
			Graph graph;
			int index{};
			for (int i = 0; i < numVertices() - 1; ++i) {
				for (int j = i + 1; j < numVertices(); ++j) {
					if (code & (1ULL << index)) graph.addEdge(i, j);
					++index;
				}
			}
			return graph;
		}
	};


	//     0 o--o 1
	//        \ |
	//         \|
	//          o 2
	inline void printGraph(const Graph<3>& graph) {
		using std::cout;
		const auto& g = graph.getAdjacencyMatrix();
		cout << "o" << (g(0, 1) ? "--" : "  ") << "o\n";
		cout << ' ' << (g(0, 2) ? '\\' : ' ') << ' ' << (g(1, 2) ? '|' : ' ') << '\n';
		cout << ' ' << ' ' << (g(0, 2) ? '\\' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
		cout << "    o\n";
	}


	//     0 o--o 1
	//       |\/|
	//       |/\|
	//     3 o--o 2
	inline void printGraph(const Graph<4>& graph) {
		using std::cout;
		const auto& g = graph.getAdjacencyMatrix();
		cout << "o" << (g(0, 1) ? "--" : "  ") << "o\n";
		cout << (g(0, 3) ? '|' : ' ') << (g(0, 2) ? '\\' : ' ') << (g(1, 3) ? '/' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
		cout << (g(0, 3) ? '|' : ' ') << (g(1, 3) ? '/' : ' ') << (g(0, 2) ? '\\' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
		cout << "o" << (g(2, 3) ? "--" : "  ") << "o\n";
	}

	//          0     1
	//           o---o     
	//           |   |\
	//           |\ /|.\
	//           | X |  o 2
	//           |/ \|./
	//           | . |/	
	//           o---o	
	//          4     3
	inline void printGraph(const Graph<5>& graph) {
		using std::cout;
		const auto& g = graph.getAdjacencyMatrix();
		cout << "   o" << (g(0, 1) ? "---" : "   ") << "o\n  ";
		cout << ' ' << (g(0, 4) ? '|' : g(0, 3) ? '\\' : ' ') << ' ' << (g(0, 2) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(1, 4) ? '/' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n' << ' ';
		cout << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << (g(0, 3) ? '\\' : ' ') << ' ' << (g(1, 4) ? '/' : ' ') << (g(1, 3) ? '|' : ' ') << (g(0, 2) ? '.' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n';
		cout << ' ' << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << ' ' << (g(0, 3) && g(1, 4) ? 'X' : (g(0, 3) ? '\\' : (g(1, 4) ? '/' : ' ')))
			<< ' ' << (g(1, 3) ? '|' : ' ') << ' ' << ' ' << "o\n ";
		cout << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << (g(1, 4) ? '/' : ' ') << ' ' << (g(0, 3) ? '\\' : ' ') << (g(1, 3) ? '|' : ' ') << (g(2, 4) ? '.' : ' ') << (g(2, 3) ? '/' : ' ') << '\n' << ' ' << ' ';
		cout << ' ' << (g(0, 4) ? '|' : g(1, 4) ? '/' : ' ') << ' ' << (g(2, 4) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(0, 3) ? '\\' : ' ') << (g(2, 3) ? '/' : ' ') << '\n';
		cout << "   o" << (g(3, 4) ? "---" : "   ") << "o\n";
	}
	//          0     1
	//           o---o             o---o            o---o
	//          /| . |\		      /| .		        |   /\
	//         /.|\ /|.\	     / |\   .	        |  /  \
	//      5 o--|-X-|--o 2     o  | \    o	     o  | /    o
	//         \.|/ \|./	       |  \		      \ |/    /
	//          \| . |/		       |   \	       \|    /
	//           o---o		       o   o	        o---o
	//          4     3
	inline void printGraph(const Graph<6>& graph) {
		using std::cout;
		const auto& g = graph.getAdjacencyMatrix();
		cout << "   o" << (g(0, 1) ? "---" : "   ") << "o\n  ";
		cout << (g(0, 5) ? '/' : ' ') << (g(0, 4) ? '|' : g(0, 3) ? '\\' : ' ') << ' ' << (g(0, 2) || g(1, 5) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(1, 4) ? '/' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n' << ' ';
		cout << (g(0, 5) ? '/' : ' ') << (g(1, 5) ? '.' : ' ') << (g(0, 4) ? '|' : ' ') << (g(0, 3) ? '\\' : ' ') << ' ' << (g(1, 4) ? '/' : ' ') << (g(1, 3) ? '|' : ' ') << (g(0, 2) ? '.' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n';
		cout << 'o' << (g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(0, 4) ? '|' : g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(0, 3) && g(1, 4) ? 'X' : (g(0, 3) ? '\\' : (g(1, 4) ? '/' : (g(2, 5) ? '-' : ' '))))
			<< (g(2, 5) ? '-' : ' ') << (g(1, 3) ? '|' : g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << "o\n ";
		cout << (g(4, 5) ? '\\' : ' ') << (g(3, 5) ? '.' : ' ') << (g(0, 4) ? '|' : ' ') << (g(1, 4) ? '/' : ' ') << ' ' << (g(0, 3) ? '\\' : ' ') << (g(1, 3) ? '|' : ' ') << (g(2, 4) ? '.' : ' ') << (g(2, 3) ? '/' : ' ') << '\n' << ' ' << ' ';
		cout << (g(4, 5) ? '\\' : ' ') << (g(0, 4) ? '|' : g(1, 4) ? '/' : ' ') << ' ' << (g(2, 4) || g(3, 5) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(0, 3) ? '\\' : ' ') << (g(2, 3) ? '/' : ' ') << '\n';
		cout << "   o" << (g(3, 4) ? "---" : "   ") << "o\n";
	}

	template<int n>
	auto generateSubgraphs(const Graph<n>& graph, int minEdges, int maxEdges) {
		std::vector<std::pair<size_t, size_t>> edges;
		std::vector<Graph<n>> subgraphs;
		for (int i = 0; i < graph.numVertices(); ++i) {
			for (int j = i + 1; j < graph.numVertices(); ++j) {
				if (graph.adjacencyMatrix(i, j) == 1) edges.push_back({ i,j });
			}
		}
		assert(edges.size() < 64 && "this algorithm only works with less than 64 edges");
		const auto end = 1 << edges.size();
		for (size_t i = 0; i < end; ++i) {
			if (const auto edgeCount = std::popcount(i); edgeCount < minEdges || edgeCount > maxEdges) continue;

			Graph<n> subgraph(graph.graphSize);
			for (size_t j = 0; j < edges.size(); ++j) {
				if (i & (1ULL << j)) {
					subgraph.adjacencyMatrix(edges[j].first, edges[j].second) = 1;
					subgraph.adjacencyMatrix(edges[j].second, edges[j].first) = 1;
				}
			}
			subgraphs.push_back(subgraph);
		}
		return subgraphs;
	}

	template<int n>
	auto generateSubgraphs(const Graph<n>& graph, int maxEdges = std::numeric_limits<int>::max()) {
		return generateSubgraphs(graph, 0, maxEdges);
	}


	namespace efficient {

		// Space- (and often time-) efficient representation using bitstrings
		template<int numVertices>
		struct Graph {
			static_assert(numVertices < 64);
			using AdjacencyMatrix = BinaryRowMatrix<numVertices, numVertices>;
			AdjacencyMatrix adjacencyMatrix;

			explicit constexpr Graph(const Q::Graph<numVertices>& graph) : adjacencyMatrix(graph.adjacencyMatrix) {}

			Q::Graph<numVertices> toGraph() const {
				Q::Graph<numVertices> graph;
				graph.adjacencyMatrix = adjacencyMatrix.toMatrix();
				return graph;
			}

			const AdjacencyMatrix& getAdjacencyMatrix() const { return adjacencyMatrix; }


		};
	}
}
