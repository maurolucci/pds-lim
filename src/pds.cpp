#include "pds.hpp"

namespace pds {

Pds::Pds() : Pds(PowerGrid{}, 0) {}

Pds::Pds(const pds::PowerGrid &graph, size_t n_channels)
    : Pds(PowerGrid{graph}, n_channels) {}

Pds::Pds(PowerGrid &&graph, size_t n_channels)
    : graph(graph),
      n_channels(n_channels),
      activated(num_vertices(graph), false),
      monitoredSet(num_vertices(graph), false),
      n_monitored(0),
      n_monitored_neighbors(num_vertices(graph), 0),
      observers(num_vertices(graph), std::set<Vertex>()) {
  for (auto v : boost::make_iterator_range(vertices(graph)))
    add_vertex(LabelledVertex{.label = v}, digraph);
}

VertexList Pds::get_monitored_set(std::map<Vertex, double> &s,
                                  std::map<Edge, double> &w) {
  VertexList monitored(num_vertices(graph), false);

  // Domiation rule
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    if (s.at(v) < 0.5) continue;
    monitored[v] = true;
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (degree(v, graph) <= n_channels || w.at(std::make_pair(v, u)) > 0.5)
        monitored[u] = true;
  }

  // Neighborhood-propagation rule
  std::list<Vertex> candidates;
  for (auto v : boost::make_iterator_range(vertices(graph)))
    if (monitored[v] && isZeroInjection(v)) candidates.push_back(v);

  while (!candidates.empty()) {
    Vertex v = candidates.front();
    candidates.pop_front();
    size_t count =
        boost::range::count_if(boost::adjacent_vertices(v, graph),
                               [monitored](auto u) { return monitored[u]; });
    if (boost::degree(v, graph) - count != 1) continue;
    auto it_u =
        boost::range::find_if(boost::adjacent_vertices(v, graph),
                              [monitored](auto u) { return !monitored[u]; });
    monitored[*it_u] = true;
    if (isZeroInjection(*it_u)) candidates.push_back(*it_u);
    for (auto y : boost::make_iterator_range(adjacent_vertices(*it_u, graph)))
      if (y != v && monitored[y] && isZeroInjection(y)) candidates.push_back(y);
  }

  return monitored;
}

void Pds::activate(Vertex v, std::vector<bool> &dominate) {
  bool act = activated[v];
  std::list<Vertex> turnedOn;
  std::list<Vertex> turnedOff;

  // Activate v
  if (!activated[v]) {
    activated[v] = true;
    if (observers[v].empty() && !propagator.contains(v)) {
      assert(!monitoredSet[v]);
      monitoredSet[v] = true;
      n_monitored++;
      for (Vertex y : boost::make_iterator_range(adjacent_vertices(v, graph)))
        n_monitored_neighbors[y]++;
      turnedOn.push_back(v);
      observers[v].insert(v);  // observe before despropagating
      despropagate_to(v, turnedOff);
    } else
      observers[v].insert(v);
  }

  // Dominate or desdominate neighbors
  size_t i = 0;
  for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
    // Dominate
    if (dominate[i] && !observers[u].contains(v)) {
      if (observers[u].empty() && !propagator.contains(u)) {
        assert(!monitoredSet[u]);
        monitoredSet[u] = true;
        n_monitored++;
        for (Vertex y : boost::make_iterator_range(adjacent_vertices(u, graph)))
          n_monitored_neighbors[y]++;
        turnedOn.push_back(u);
        observers[u].insert(v);  // observe before despropagating
        despropagate_to(u, turnedOff);
      } else
        observers[u].insert(v);
    }
    // Desdominate
    else if (act && !dominate[i] && observers[u].contains(v)) {
      observers[u].erase(v);
      if (observers[u].empty() && !propagator.contains(u)) {
        assert(monitoredSet[u]);
        monitoredSet[u] = false;
        n_monitored--;
        for (Vertex y : boost::make_iterator_range(adjacent_vertices(u, graph)))
          n_monitored_neighbors[y]--;
        turnedOff.push_back(u);
        despropagate_from(u, turnedOff);
      }
    }
    ++i;
  }

  // Try propagations to turned off vertices
  propagate_to(turnedOff, turnedOn);

  // Try propagations from turned on vertices, or their neighbors,
  // as long as they are propagating and monitored
  std::list<Vertex> candidates;
  for (auto u : turnedOn) {
    if (isZeroInjection(u)) candidates.push_back(u);
    for (auto y : boost::make_iterator_range(adjacent_vertices(u, graph)))
      if (isZeroInjection(y) && isMonitored(y)) candidates.push_back(y);
  }
  propagate_from(candidates, turnedOn);
}

void Pds::deactivate(Vertex v) {
  if (!activated[v]) return;

  std::list<Vertex> turnedOn;
  std::list<Vertex> turnedOff;

  // Deactivate v
  activated[v] = false;
  observers[v].erase(v);
  if (observers[v].empty() && !propagator.contains(v)) {
    assert(monitoredSet[v]);
    monitoredSet[v] = false;
    n_monitored--;
    for (Vertex y : boost::make_iterator_range(adjacent_vertices(v, graph)))
      n_monitored_neighbors[y]--;
    turnedOff.push_back(v);
    despropagate_from(v, turnedOff);
  }

  // Desdominate neighbors
  for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
    if (!observers[u].contains(v)) continue;
    observers[u].erase(v);
    if (observers[u].empty() && !propagator.contains(u)) {
      assert(monitoredSet[u]);
      monitoredSet[u] = false;
      n_monitored--;
      for (Vertex y : boost::make_iterator_range(adjacent_vertices(u, graph)))
        n_monitored_neighbors[y]--;
      turnedOff.push_back(u);
      despropagate_from(u, turnedOff);
    }
  }

  // Try propagations to turned off vertices
  propagate_to(turnedOff, turnedOn);
}

void Pds::despropagate_to(Vertex to, std::list<Vertex> &turnedOff) {
  auto it = propagator.find(to);
  if (it == propagator.end()) return;
  despropagate(it->second, to, turnedOff);
}

void Pds::despropagate_from(Vertex v, std::list<Vertex> &turnedOff) {
  std::list<Vertex> candidates;
  candidates.push_back(v);
  while (!candidates.empty()) {
    Vertex u = candidates.front();
    candidates.pop_front();

    // u cannot propagate anymore
    if (propagates.contains(u)) {
      Vertex v = propagates[u];
      despropagate(u, v, turnedOff);
      candidates.push_back(v);
    }

    // u cannot be involved in a propagation
    while (out_degree(u, digraph) > 0) {
      auto e = *out_edges(u, digraph).first;
      Vertex v = target(e, digraph);
      // v can no longer be propagated
      despropagate_to(v, turnedOff);
      candidates.push_back(v);
    }
  }
}

void Pds::despropagate(Vertex from, Vertex to, std::list<Vertex> &turnedOff) {
  if (observers[to].empty()) {
    assert(monitoredSet[to]);
    monitoredSet[to] = false;
    n_monitored--;
    for (Vertex y : boost::make_iterator_range(adjacent_vertices(to, graph)))
      n_monitored_neighbors[y]--;
    turnedOff.push_back(to);
  }
  propagates.erase(from);
  propagator.erase(to);
  remove_edge(from, to, digraph);
  for (auto y : boost::make_iterator_range(adjacent_vertices(from, graph)))
    if (y != to) remove_edge(y, to, digraph);
}

bool Pds::try_propagation_to(Vertex v, std::list<Vertex> &turnedOn) {
  for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
    if (check_propagation(u, v)) {
      propagate(u, v, turnedOn);
      return true;
    }
  return false;
}

bool Pds::try_propagation_from(Vertex v, std::list<Vertex> &turnedOn) {
  for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
    if (check_propagation(v, u)) {
      propagate(v, u, turnedOn);
      return true;
    }
  return false;
}

bool Pds::check_propagation(Vertex from, Vertex to) {
  if (monitoredSet[to] || !monitoredSet[from] || !isZeroInjection(from) ||
      propagates.contains(from) ||
      n_monitored_neighbors[from] != degree(from, graph) - 1)
    return false;
  return true;
}

void Pds::propagate_to(std::list<Vertex> &candidates,
                       std::list<Vertex> &turnedOn) {
  bool keepGoing = true;
  while (keepGoing) {
    keepGoing = false;
    std::set<Vertex> propagated;
    for (auto v : candidates) {
      if (!try_propagation_to(v, turnedOn)) continue;
      keepGoing = true;
      propagated.insert(v);
      break;
    }
    for (auto v : propagated) candidates.remove(v);
  }
}

void Pds::propagate_from(std::list<Vertex> &candidates,
                         std::list<Vertex> &turnedOn) {
  while (!candidates.empty()) {
    Vertex v = candidates.front();
    candidates.pop_front();
    if (!try_propagation_from(v, turnedOn)) continue;
    Vertex u = propagates[v];
    if (isZeroInjection(u)) candidates.push_back(u);
    for (auto y : boost::make_iterator_range(adjacent_vertices(u, graph)))
      if (y != v && isZeroInjection(y) && monitoredSet[y])
        candidates.push_back(y);
  }
}

void Pds::propagate(Vertex from, Vertex to, std::list<Vertex> &turnedOn) {
  assert(!monitoredSet[to]);
  monitoredSet[to] = true;
  n_monitored++;
  for (Vertex y : boost::make_iterator_range(adjacent_vertices(to, graph)))
    n_monitored_neighbors[y]++;
  turnedOn.push_back(to);
  propagates[from] = to;
  propagator[to] = from;
  add_edge(from, to, digraph);
  for (auto v : boost::make_iterator_range(adjacent_vertices(from, graph)))
    if (v != to) add_edge(v, to, digraph);
}

bool Pds::check_get_monitored_set(std::map<Vertex, double> &s,
                                  std::map<Edge, double> &w) {
  VertexList mS = get_monitored_set(s, w);
  return (mS == monitoredSet);
}

}  // end of namespace pds