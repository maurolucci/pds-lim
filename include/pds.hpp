#ifndef PDS_HPP
#define PDS_HPP
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <range/v3/all.hpp>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/range/adaptor/filtered.hpp"
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
  size_t totalCallback;
  size_t totalCallbackTime;
  size_t totalLazy;
};

class Pds {
private:
  PowerGrid graph;
  size_t n_channels;
  VertexList activated;
  VertexList monitoredSet;
  size_t n_monitored;
  std::vector<size_t> n_adj_monitored;
  std::vector<std::vector<bool>> observed_by; // observed_by[u][v]: is u observed by v?
  std::vector<size_t> n_observers;
  std::vector<int> propagates;
  std::vector<int> propagator;
  PrecedenceDigraph digraph;

  void activate_neighbor(Vertex from, Vertex to, std::list<Vertex> &turnedOn, std::list<Vertex> &turnedOff);
  void deactivate_neighbor(Vertex from, Vertex to, std::list<Vertex> &turnedOff);

  void despropagate_to(Vertex to, std::list<Vertex> &turnedOff);
  void despropagate_from(Vertex v, std::list<Vertex> &turnedOff);
  void despropagate(Vertex from, Vertex to, std::list<Vertex> &turnedOff);

  [[nodiscard]] bool try_propagation_to(Vertex v, std::list<Vertex> &turnedOn);
  [[nodiscard]] bool try_propagation_from(Vertex v, std::list<Vertex> &turnedOn);
  [[nodiscard]] bool check_propagation(Vertex from, Vertex to);
  void propagate(Vertex from, Vertex to, std::list<Vertex> &turnedOn);

public:
  Pds();
  explicit Pds(PowerGrid &&graph, size_t n_channels);
  explicit Pds(const PowerGrid &graph, size_t n_channels);

  [[nodiscard]] inline const PowerGrid &get_graph() const { return graph; }
  [[nodiscard]] inline size_t get_n_channels() const { return n_channels; }

  [[nodiscard]] inline bool isZeroInjection(Vertex v) const {
    return graph[v].zero_injection;
  }

  [[nodiscard]] inline std::ptrdiff_t numZeroInjection() const {
    return boost::range::count_if(
        vertices(get_graph()), [this](auto v) { return isZeroInjection(v); });
  }

  [[nodiscard]] inline bool isMonitored(Vertex v) const {
    return monitoredSet[v];
  }

  [[nodiscard]] inline bool isActivated(Vertex v) const {
    return activated[v];
  }

  VertexList get_monitored_set(std::map<Vertex, double> &s,
                               std::map<Edge, double> &w);

  inline void get_unmonitored_set(std::vector<Vertex> &vec) const {
    boost::copy(vertices(graph) |
                    boost::adaptors::filtered(
                        [this](auto v) { return !monitoredSet[v]; }),
                std::back_inserter(vec));
  }

  [[nodiscard]] inline bool isFeasible(VertexList &mS) const {
    return boost::range::count_if(vertices(graph), [mS](auto u) {
             return mS[u];
           }) == static_cast<std::ptrdiff_t>(boost::num_vertices(graph));
  }

  [[nodiscard]] inline bool isFeasible() const {
    return n_monitored == boost::num_vertices(graph);
  }

  void activate(Vertex v, std::vector<Vertex> &neighbors, 
    std::list<Vertex> &turnedOn, std::list<Vertex> &turnedOff);
  void deactivate(Vertex v, std::list<Vertex> &turnedOff);

  void propagate_to(std::list<Vertex> &candidates, std::list<Vertex> &turnedOn);
  void propagate_from(std::list<Vertex> &candidates, std::list<Vertex> &turnedOn);

  [[nodiscard]] bool check_get_monitored_set(std::map<Vertex, double> &s,
                               std::map<Edge, double> &w);

}; // end of class Pds

} // end of namespace pds

#endif // PDS_HPP
