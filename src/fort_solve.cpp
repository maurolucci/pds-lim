#include "fort_solve.hpp"
#include "gurobi_common.hpp"

#include <boost/range/adaptor/filtered.hpp>
#include <gurobi_c++.h>

namespace pds {

using Fort = std::pair<std::set<Vertex>, std::set<Edge>>;

namespace {

struct LazyFortCB : public GRBCallback {

  MIPModel mipmodel;
  GRBModel &model;
  std::map<pds::Vertex, GRBVar> &s;
  std::map<Vertex, double> sValue;
  std::set<Vertex> pmuSet;
  std::map<Edge, GRBVar> &w;
  std::map<Edge, double> wValue;
  std::map<Edge, GRBVar> y;
  Pds &input;
  const PowerGrid &graph;
  std::map<Edge, EdgeList> translate;
  std::ostream &cbFile, &solFile;
  size_t n_channels;
  size_t lazyLimit;

  LazyFortCB(Pds &input, std::ostream &callbackFile, std::ostream &solutionFile,
             size_t lzLimit)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()), cbFile(callbackFile),
        solFile(solutionFile), n_channels(input.get_n_channels()),
        lazyLimit(lzLimit), observers(num_vertices(graph), std::set<Vertex>()) {

    model.setCallback(this);

    // Add variables
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      s.try_emplace(
          v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
        w.try_emplace(std::make_pair(v, u),
                      model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                   fmt::format("w_{}_{}", v, u)));
      }
    }

    // (2) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      GRBLinExpr constr3 = 0;
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
        constr3 += w.at(std::make_pair(v, u));
      }
      model.addConstr(constr3 <=
                      std::min(degree(v, graph), (n_channels - 1)) *
                          s.at(v));
    }

    // Initialize precedense digraph
    for (auto v : boost::make_iterator_range(vertices(graph)))
      add_vertex(LabelledVertex{.label = v}, digraph);

  }

  SolveResult solve(boost::optional<std::string> logPath, double timeLimit) {
    if (logPath.has_value()) {
      model.set(GRB_IntParam_LogToConsole, false);
    }
    model.set(GRB_IntParam_LazyConstraints, 1);
    return solveMIP(input, mipmodel, logPath, solFile, timeLimit);
  }

  void callback() override {

    switch (where) {

    // MIP solution callback
    // Integer solution founded (it does not necessarily improve the
    // incumbent)
    case GRB_CB_MIPSOL:

      // Recover variable values
      pmuSet.clear();
      for (auto v : boost::make_iterator_range(vertices(graph))) {
        sValue[v] = getSolution(s.at(v));
        if (sValue[v] > 0.5)
          pmuSet.insert(v);
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          wValue[std::make_pair(v, u)] =
              getSolution(w.at(std::make_pair(v, u)));
      }

      // Feasibility check
      VertexList mS = input.get_monitored_set(sValue, wValue);
      if (!input.isFeasible(mS)) {

        // Find violated cycles
        std::set<Fort> forts = violatedForts(mS, lazyLimit);
        std::pair<double, double> avg = addLazyForts(forts);

        // Report to callback file
        cbFile << fmt::format("# forts: {}, avg. vertex size: {:.2f}, avg. "
                              "edge size: {:.2f}",
                              forts.size(), avg.first, avg.second)
               << std::endl;
      }

      break;
    }
  }

private:

  std::vector<std::set<Vertex>> observers;
  std::map<Vertex, Vertex> propagates;
  std::map<Vertex, Vertex> propagator;
  PrecedenceDigraph digraph;

  std::set<Fort> violatedForts(const VertexList &monitoredSet,
                               size_t fortsLimit) {

    std::set<Fort> forts;

    // Get unmonitored set
    VertexList unmonitoredSet;
    boost::copy(vertices(graph) |
                    boost::adaptors::filtered(
                        [monitoredSet](auto v) { return !monitoredSet[v]; }),
                std::back_inserter(unmonitoredSet));

    // Shuffle unmonitored set
    boost::range::random_shuffle(unmonitoredSet);

    // Preselect random neighbors for unmonitored vertices
    for (Vertex v : unmonitoredSet) {
      VertexList neighbors; 
      boost::copy(adjacent_vertices(v, graph) |
                      boost::adaptors::filtered([unmonitoredSet](auto v) {
                        return boost::range::find(unmonitoredSet, v) !=
                               unmonitoredSet.end();
                      }),
                  std::back_inserter(neighbors));
      size_t k = std::min(neighbors.size(), n_channels - 1);
      if (k < neighbors.size())
        random_unique(neighbors.begin(), neighbors.end(), k);
      for (size_t i = 0; i < k; ++i)
        wValue.at(std::make_pair(v, neighbors[i])) = 1.0;      
    }

    // Copy of monitored set
    VertexList mS (monitoredSet);

    // Activate all unmonitored vertices
    // (propagation is unnecessary here)
    for (Vertex v : unmonitoredSet)
      (void) activate(v, mS);

    // MINIMISE FEASIBLE SOLUTION
    for (Vertex v : unmonitoredSet) {

      if (forts.size() >= fortsLimit)
        break;

      // Deactivate v
      std::list<Vertex> changed = deactivate(v, mS);
      
      // Changed vertices can no longer propagate
      std::set<Vertex> changed2 = despropagate(changed, mS);

      // Try propagations to changed vertices
      propagate_to(changed2, mS);

      // Feasibility check
      VertexList mS2 = input.get_monitored_set(sValue, wValue);
      std::cout << fmt::format("mS: {},\t mS2: {}\n", mS, mS2);
      assert(mS == mS2);
      if (!input.isFeasible(mS)) {
        
        // Find and insert fort
        forts.insert(findFort(mS));

        // Reactivate v
        std::list<Vertex> changed3 = activate(v, mS);

        // Changed vertices no longer need to be propagated
        for (auto u: changed3)
          (void) try_despropagate(u, mS);

        // Try propagations from changed vertices or their neighbors
        std::set<Vertex> changed4; 
        for (auto u: changed3) {
          changed4.insert(u);
          for (auto y: boost::make_iterator_range(adjacent_vertices(u, graph)))
            changed4.insert(y);
        }
        std::list<Vertex> changed5 (changed4.begin(), changed4.end());
        propagate_from(changed5, mS);

      }

    }

    return forts;
  }

  std::list<Vertex> activate(Vertex v, VertexList &monitoredSet) {
    std::list<Vertex> changed; // Vertices that are no longer unmonitored
    sValue.at(v) = 1.0;
    if (observers[v].empty()) {
      monitoredSet[v] = true;
      changed.push_back(v);
    }
    observers[v].insert(v);
    for (auto u: boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (wValue.at(std::make_pair(v,u)) > 0.5) {
        if (observers[u].empty()) {
          monitoredSet[u] = true;
          changed.push_back(u);
        }
        observers[u].insert(v);
      }
    return changed;
  }

  std::list<Vertex> deactivate(Vertex v, VertexList &monitoredSet) {
    std::list<Vertex> changed; // Vertices that are no longer monitored
    // Deactivate v
    sValue.at(v) = 0.0;
    observers[v].erase(v);
    if (observers[v].empty()) {
      monitoredSet[v] = false;
      changed.push_back(v);
    }
    // Deactivate neighbors of v
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (wValue.at(std::make_pair(v, u)) > 0.5) {
        observers[u].erase(v);
        if (observers[u].empty()) {
          monitoredSet[u] = false;
          changed.push_back(u);
        }
      }
    // Quick fix
    std::list<Vertex> changed2;
    for (auto v: changed)
      if (!try_propagation_to(v, monitoredSet))
        changed2.push_back(v);
    return changed2;
  }

  bool try_propagation_to(Vertex v, VertexList &monitoredSet) {
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (check_propagation(u, v, monitoredSet)) {
        propagate(u, v, monitoredSet);
        return true;
      }
    return false;
  }

  bool try_propagation_from(Vertex v, VertexList &monitoredSet) {
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (check_propagation(v, u, monitoredSet)) {
        propagate(v, u, monitoredSet);
        return true;
      }
    return false;
  }

  bool check_propagation(Vertex from, Vertex to, VertexList &monitoredSet) {
    if (!monitoredSet[from] || !input.isZeroInjection(from))
      return false;
    size_t count = boost::range::count_if(
        adjacent_vertices(from, graph),
        [monitoredSet, to](auto v) { return v != to && monitoredSet[v]; });
    return (degree(from, graph) - count == 1);
  }

  void propagate_to(std::set<Vertex> &candidates, VertexList &monitoredSet) {
    bool keepGoing = true;
    while (keepGoing) {
      keepGoing = false;
      std::set<Vertex> propagated;
      for (auto v : candidates) {
        if (!try_propagation_to(v, monitoredSet))
          continue;
        keepGoing = true;
        propagated.insert(v);
        break;
      }
      for (auto v : propagated)
        candidates.erase(v);
    }
  }

  void propagate_from(std::list<Vertex> &candidates, VertexList &monitoredSet) {
    VertexList revised (num_vertices(graph), false);
    while (!candidates.empty()) {
      Vertex v = candidates.front();
      candidates.pop_front();
      if (revised[v])
        continue;
      revised[v] = true;
      if (!try_propagation_from(v, monitoredSet))
        continue;
      Vertex u = propagates[v];
      candidates.push_back(u);
      revised[u] = false;
      for (auto y: boost::make_iterator_range(adjacent_vertices(u, graph)))
        if (y != v) {
          candidates.push_back(y);
          revised[y] = false;
        }
    }
  }

  void propagate(Vertex from, Vertex to, VertexList &monitoredSet) {
    monitoredSet[to] = true;
    propagates[from] = to;
    propagator[to] = from;
    add_edge(from, to, digraph);
    for (auto v : boost::make_iterator_range(adjacent_vertices(from, graph)))
      if (v != to)
        add_edge(v, to, digraph);
  }

  bool try_despropagate(Vertex v, VertexList &monitoredSet) {
    if (propagator.contains(v)) {
      despropagate(propagator[v], v, monitoredSet);
      return true;
    }
    return false;
  }

  std::set<Vertex> despropagate(std::list<Vertex> &candidates, VertexList &monitoredSet) {
    std::set<Vertex> changed;
    while (!candidates.empty()) {
      Vertex u = candidates.front();
      candidates.pop_front();
      changed.insert(u);

      // u cannot propagate anymore
      if (propagates.contains(u)) {
        Vertex v = propagates[u];
        despropagate(u, v, monitoredSet);
        candidates.push_back(v);     
      }

      // u cannot be involved in a propagation
      for (auto e : boost::make_iterator_range(out_edges(u, digraph))) {
        Vertex v = target(e, digraph);
        // v can no longer be propagated
        Vertex y = propagator[v];
        despropagate(y, v, monitoredSet);
        candidates.push_back(v);
      }
    }
    return changed;
  }

  void despropagate(Vertex from, Vertex to, VertexList &monitoredSet) {
    if (observers[to].empty())
      monitoredSet[to] = false;
    propagates.erase(from);
    propagator.erase(to);
    remove_edge(from, to, digraph);
    for (auto y : boost::make_iterator_range(adjacent_vertices(from, graph)))
      if (y != to)
        remove_edge(y, to, digraph);
  }

  template <class bidiiter>
  bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
      bidiiter r = begin;
      std::advance(r, rand() % left);
      std::swap(*begin, *r);
      ++begin;
      --left;
    }
    return begin;
  }

  Fort findFort(const VertexList &monitoredSet) {
    Fort fort;
    for (auto u : boost::make_iterator_range(vertices(graph))) {
      if (monitoredSet[u])
        continue;
      fort.first.insert(u);
      for (auto z : boost::make_iterator_range(adjacent_vertices(u, graph))) {
        if (!monitoredSet[z])
          continue;
        if (!pmuSet.contains(z))
          fort.first.insert(z);
        else
          fort.second.insert(std::make_pair(z, u));
      }
    }
    return fort;
  }

  std::pair<double, double> addLazyForts(std::set<Fort> forts) {

    size_t accumVertices = 0;
    size_t accumEdges = 0;

    for (auto &f : forts) {
      GRBLinExpr fortSum;
      for (auto v : f.first) {
        fortSum += s.at(v);
        accumVertices++;
      }
      for (auto e : f.second) {
        fortSum += w.at(e);
        accumEdges++;
      }
      addLazy(fortSum >= 1);
    }

    return std::make_pair(static_cast<double>(accumVertices) / forts.size(),
                          static_cast<double>(accumEdges) / forts.size());
  }
};
} // namespace

SolveResult solveLazyForts(Pds &input, boost::optional<std::string> logPath,
                           std::ostream &callbackFile, std::ostream &solFile,
                           double timeLimit, size_t lazyLimit) {
  LazyFortCB lazyForts(input, callbackFile, solFile, lazyLimit);
  return lazyForts.solve(logPath, timeLimit);
}

} // end of namespace pds