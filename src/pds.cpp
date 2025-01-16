#include "pds.hpp"

namespace pds {

Pds::Pds() : Pds(PowerGrid{}, 0) {}

Pds::Pds(const pds::PowerGrid &graph, size_t n_channels)
    : Pds(PowerGrid{graph}, n_channels) {}

Pds::Pds(PowerGrid &&graph, size_t n_channels)
    : graph(graph), n_channels(n_channels) {}

VertexList Pds::get_monitored_set(std::map<Vertex, double> &s,
                                  std::map<Edge, double> &w) {

  VertexList monitored(num_vertices(graph), false);

  // Domiation rule
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    if (s.at(v) < 0.5)
      continue;
    monitored[v] = true;
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
      if (w.at(std::make_pair(v, u)) < 0.5)
        continue;
      monitored[u] = true;
    }
  }

  // Neighborhood-propagation rule
  /*
  bool stop = false;
  while (!stop) {
    stop = true;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (!monitored[v] || !isZeroInjection(v))
        continue;
      size_t count =
          boost::range::count_if(boost::adjacent_vertices(v, graph),
                                 [monitored](auto u) { return monitored[u]; });
      if (boost::degree(v, graph) - count != 1)
        continue;
      auto it_u =
          boost::range::find_if(boost::adjacent_vertices(v, graph),
                                [monitored](auto u) { return !monitored[u]; });
      monitored[*it_u] = true;
      stop = false;
    }
  }*/

  std::list<Vertex> candidates;
  for (auto v : boost::make_iterator_range(vertices(graph)))
    if (monitored[v] && isZeroInjection(v))
      candidates.push_back(v);

  while (!candidates.empty()) {
    Vertex v = candidates.front();
    candidates.pop_front();
    size_t count =
        boost::range::count_if(boost::adjacent_vertices(v, graph),
                                [monitored](auto u) { return monitored[u]; });
    if (boost::degree(v, graph) - count != 1)
      continue;
    auto it_u =
        boost::range::find_if(boost::adjacent_vertices(v, graph),
                              [monitored](auto u) { return !monitored[u]; });
    monitored[*it_u] = true;
    if (isZeroInjection(*it_u))
      candidates.push_back(*it_u);
    for (auto y: boost::make_iterator_range(adjacent_vertices(*it_u, graph)))
      if (y != v && isZeroInjection(y)) 
        candidates.push_back(y);
  }

  return monitored;
}

} // end of namespace pds