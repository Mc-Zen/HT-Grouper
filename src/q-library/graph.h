
#pragma once
#include "efficient_binary_math.h"
#include <iostream>
#include <vector>

namespace Q {




	template<bool empty = true>
	struct GraphSize {
		constexpr friend bool operator==(const GraphSize& a, const GraphSize& b) = default;
	};


	template<>
	struct GraphSize<false> {
		int n{};
		constexpr friend bool operator==(const GraphSize& a, const GraphSize& b) = default;
	};


	template<size_t n = Math::dynamic>
	class Graph {
	public:
		using AdjacencyMatrix = Math::Matrix<Binary, n, n>;
		static constexpr bool is_dynamic = (n == Math::dynamic);
		GraphSize<!is_dynamic> graphSize;


		AdjacencyMatrix adjacencyMatrix;


		constexpr Graph() requires !is_dynamic = default;

		explicit constexpr Graph(GraphSize<!is_dynamic> graphSize) : graphSize(graphSize) {
			if constexpr (is_dynamic) {
				adjacencyMatrix = AdjacencyMatrix(numVertices(), numVertices());
			}
		}

		explicit constexpr Graph(int numVertices) requires is_dynamic : graphSize{ numVertices }, adjacencyMatrix(numVertices, numVertices) {}


		constexpr int numVertices() const {
			if constexpr (is_dynamic) return graphSize.n;
			else return n;
		}

		constexpr static auto fullyConnected() requires (!is_dynamic) { return Graph{}.fullyConnect(); }
		constexpr static auto fullyConnected(int n) requires (is_dynamic) { return Graph{ n }.fullyConnect(); }

		constexpr static auto star(int center = 0) requires (!is_dynamic) { return Graph{}.makeStar(center); }
		constexpr static auto star(int n, int center = 0) requires (is_dynamic) { return Graph{ n }.makeStar(center); }

		constexpr static auto linear() requires (!is_dynamic) { return Graph{}.makeLinear(); }
		constexpr static auto linear(int n) requires (is_dynamic) { return Graph{ n }.makeLinear(); }

		constexpr static auto cycle() requires (!is_dynamic) { return Graph{}.makeCycle(); }
		constexpr static auto cycle(int n) requires (is_dynamic) { return Graph{ n }.makeCycle(); }

		constexpr static auto pusteblume() requires (!is_dynamic) { return Graph{}.makePusteblume(); }
		constexpr static auto pusteblume(int n) requires (is_dynamic) { return Graph{ n }.makePusteblume(); }

		constexpr const AdjacencyMatrix& getAdjacencyMatrix() const { return adjacencyMatrix; }

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
		constexpr auto getEdges() const {
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
		static constexpr int64_t compress(const Graph& graph) {
			static_assert(n * (n - 1) / 2 <= 64 || n != Math::dynamic, "Compression is not supported for graphs of this size");
			assert(graph.numVertices() * (graph.numVertices() - 1) / 2 <= 64 && "Compression is not supported for graphs of this size");
			int64_t code{};
			int index{};
			for (int i = 0; i < graph.numVertices() - 1; ++i) {
				for (int j = i + 1; j < graph.numVertices(); ++j) {
					if (graph.hasEdge(i, j)) code |= (1ULL << index);
					++index;
				}
			}
			return code;
		}

		/// @brief Restore a graph from its compressed form. 
		static constexpr Graph decompress(int64_t code) requires(!is_dynamic) {
			static_assert(n * (n - 1) / 2 <= 64, "Deompression is not supported for graphs of this size");
			Graph graph;
			decompressImpl(graph, code);
			return graph;
		}

		/// @brief Restore a graph from its compressed form. 
		static constexpr Graph decompress(int numVertices, int64_t code) requires(is_dynamic) {
			assert(numVertices * (numVertices - 1) / 2 <= 64 && "Deompression is not supported for graphs of this size");
			Graph graph(numVertices);
			decompressImpl(graph, code);
			return graph;
		}


	private:

		static constexpr void decompressImpl(Graph& graph, int64_t code) {
			int index{};
			for (int i = 0; i < graph.numVertices() - 1; ++i) {
				for (int j = i + 1; j < graph.numVertices(); ++j) {
					if (code & (1ULL << index)) graph.addEdge(i, j);
					++index;
				}
			}
		}

		constexpr Graph& fullyConnect() {
			adjacencyMatrix.fill(1);
			for (int i = 0; i < numVertices(); ++i) adjacencyMatrix(i, i) = 0;
			return *this;
		}

		constexpr Graph& makeStar(int center) {
			std::for_each(adjacencyMatrix.col_begin(center), adjacencyMatrix.col_end(center), [](Binary& el) { el.negate(); });
			std::for_each(adjacencyMatrix.row_begin(center), adjacencyMatrix.row_end(center), [](Binary& el) { el.negate(); });
			return *this;
		}

		constexpr Graph& makeLinear() {
			for (int i = 0; i < numVertices() - 1; ++i) {
				addEdge(i, i + 1);
			}
			return *this;
		}

		constexpr Graph& makeCycle() {
			makeLinear();
			addEdge(0, numVertices() - 1);
			return *this;
		}

		constexpr Graph& makePusteblume() {
			static_assert(n >= 5 || n == Math::dynamic, "The Pusteblume graph is only possible for at least 5 vertices");
			assert(numVertices() >= 5 && "The Pusteblume graph is only possible for at least 5 vertices");
			for (size_t i = 1; i < 4; ++i) {
				addEdge(0, i);
			}
			for (size_t i = 4; i < numVertices(); ++i) {
				addEdge(3, i);
			}
			return *this;
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

	inline void printGraph(const Graph<Math::dynamic>& graph) {
		using std::cout;
		const auto& g = graph.getAdjacencyMatrix();
		switch (graph.numVertices()) {
		case 3:
			cout << "o" << (g(0, 1) ? "--" : "  ") << "o\n";
			cout << ' ' << (g(0, 2) ? '\\' : ' ') << ' ' << (g(1, 2) ? '|' : ' ') << '\n';
			cout << ' ' << ' ' << (g(0, 2) ? '\\' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
			cout << "    o\n";
			break;
		case 4:
			cout << "o" << (g(0, 1) ? "--" : "  ") << "o\n";
			cout << (g(0, 3) ? '|' : ' ') << (g(0, 2) ? '\\' : ' ') << (g(1, 3) ? '/' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
			cout << (g(0, 3) ? '|' : ' ') << (g(1, 3) ? '/' : ' ') << (g(0, 2) ? '\\' : ' ') << (g(1, 2) ? '|' : ' ') << '\n';
			cout << "o" << (g(2, 3) ? "--" : "  ") << "o\n";
			break;
		case 5:
			cout << "   o" << (g(0, 1) ? "---" : "   ") << "o\n  ";
			cout << ' ' << (g(0, 4) ? '|' : g(0, 3) ? '\\' : ' ') << ' ' << (g(0, 2) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(1, 4) ? '/' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n' << ' ';
			cout << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << (g(0, 3) ? '\\' : ' ') << ' ' << (g(1, 4) ? '/' : ' ') << (g(1, 3) ? '|' : ' ') << (g(0, 2) ? '.' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n';
			cout << ' ' << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << ' ' << (g(0, 3) && g(1, 4) ? 'X' : (g(0, 3) ? '\\' : (g(1, 4) ? '/' : ' ')))
				<< ' ' << (g(1, 3) ? '|' : ' ') << ' ' << ' ' << "o\n ";
			cout << ' ' << ' ' << (g(0, 4) ? '|' : ' ') << (g(1, 4) ? '/' : ' ') << ' ' << (g(0, 3) ? '\\' : ' ') << (g(1, 3) ? '|' : ' ') << (g(2, 4) ? '.' : ' ') << (g(2, 3) ? '/' : ' ') << '\n' << ' ' << ' ';
			cout << ' ' << (g(0, 4) ? '|' : g(1, 4) ? '/' : ' ') << ' ' << (g(2, 4) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(0, 3) ? '\\' : ' ') << (g(2, 3) ? '/' : ' ') << '\n';
			cout << "   o" << (g(3, 4) ? "---" : "   ") << "o\n";
		case 6:
			cout << "   o" << (g(0, 1) ? "---" : "   ") << "o\n  ";
			cout << (g(0, 5) ? '/' : ' ') << (g(0, 4) ? '|' : g(0, 3) ? '\\' : ' ') << ' ' << (g(0, 2) || g(1, 5) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(1, 4) ? '/' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n' << ' ';
			cout << (g(0, 5) ? '/' : ' ') << (g(1, 5) ? '.' : ' ') << (g(0, 4) ? '|' : ' ') << (g(0, 3) ? '\\' : ' ') << ' ' << (g(1, 4) ? '/' : ' ') << (g(1, 3) ? '|' : ' ') << (g(0, 2) ? '.' : ' ') << (g(1, 2) ? '\\' : ' ') << '\n';
			cout << 'o' << (g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(0, 4) ? '|' : g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(0, 3) && g(1, 4) ? 'X' : (g(0, 3) ? '\\' : (g(1, 4) ? '/' : (g(2, 5) ? '-' : ' '))))
				<< (g(2, 5) ? '-' : ' ') << (g(1, 3) ? '|' : g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << (g(2, 5) ? '-' : ' ') << "o\n ";
			cout << (g(4, 5) ? '\\' : ' ') << (g(3, 5) ? '.' : ' ') << (g(0, 4) ? '|' : ' ') << (g(1, 4) ? '/' : ' ') << ' ' << (g(0, 3) ? '\\' : ' ') << (g(1, 3) ? '|' : ' ') << (g(2, 4) ? '.' : ' ') << (g(2, 3) ? '/' : ' ') << '\n' << ' ' << ' ';
			cout << (g(4, 5) ? '\\' : ' ') << (g(0, 4) ? '|' : g(1, 4) ? '/' : ' ') << ' ' << (g(2, 4) || g(3, 5) ? '.' : ' ') << ' ' << (g(1, 3) ? '|' : g(0, 3) ? '\\' : ' ') << (g(2, 3) ? '/' : ' ') << '\n';
			cout << "   o" << (g(3, 4) ? "---" : "   ") << "o\n";
			break;
		default:
			cout << "Cannot print graph of this size";
		}
	}


	/// @brief Generate all subgraphs of given graph that have at least [minEdges] edges and at most [maxEdges] edges
	template<int n>
	std::vector<Graph<n>> generateSubgraphs(const Graph<n>& graph, int minEdges, int maxEdges) {
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

			constexpr Q::Graph<numVertices> toGraph() const {
				Q::Graph<numVertices> graph;
				graph.adjacencyMatrix = adjacencyMatrix.toMatrix();
				return graph;
			}

			constexpr const AdjacencyMatrix& getAdjacencyMatrix() const { return adjacencyMatrix; }
		};


	}

}
