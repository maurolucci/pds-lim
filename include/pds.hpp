#ifndef PDS_HPP
#define PDS_HPP
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <range/v3/all.hpp>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/range/algorithm.hpp"

namespace pds {

struct Bus {
  std::string name;
  long id;
  bool zero_injection;
};

using PowerGrid =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Bus>;
using Vertex = PowerGrid::vertex_descriptor;
using Edge = std::pair<Vertex, Vertex>;
using VertexList = std::vector<Vertex>;
using EdgeList = std::list<Edge>;

struct LabelledVertex {
  Vertex label;
};

using PrecedenceDigraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
                          LabelledVertex>;
using Node = PrecedenceDigraph::vertex_descriptor;

enum class SolveState { Optimal, Timeout, Infeasible, Other };

struct SolveResult {
  int variables;
  int constraints;
  double lower;
  double upper;
  double gap;
  double nodes;
  SolveState state;
};

class Pds {
private:
  PowerGrid graph;
  size_t n_channels;

public:
  Pds();
  explicit Pds(PowerGrid &&graph, size_t n_channels);
  explicit Pds(const PowerGrid &graph, size_t n_channels);

  [[nodiscard]] inline const PowerGrid &get_graph() const { return graph; }
  [[nodiscard]] inline size_t get_n_channels() const { return n_channels; }

  [[nodiscard]] inline bool
  isZeroInjection(PowerGrid::vertex_descriptor v) const {
    return graph[v].zero_injection;
  }

  [[nodiscard]] inline std::ptrdiff_t numZeroInjection() const {
    return boost::range::count_if(
        vertices(get_graph()), [this](auto v) { return isZeroInjection(v); });
  }

  VertexList get_monitored_set(std::map<Vertex, double> &s,
                               std::map<Edge, double> &w);

  [[nodiscard]] inline bool isFeasible(VertexList &mS) const {
    return boost::range::count_if(vertices(graph), [mS](auto u) {
             return mS[u];
           }) == static_cast<std::ptrdiff_t>(boost::num_vertices(graph));
  }

}; // end of class Pds

} // end of namespace pds

#endif // PDS_HPP