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
        lazyLimit(lzLimit) {

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
      for (auto v : boost::make_iterator_range(vertices(graph))) {
        sValue[v] = getSolution(s.at(v));
        if (sValue[v] > 0.5)
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          wValue[std::make_pair(v, u)] =
              getSolution(w.at(std::make_pair(v, u)));
      }

      // Update solution

      // First deactivate vertices
      std::list<Vertex> turnedOff, turnedOn;
      for (auto v : boost::make_iterator_range(vertices(graph)))
        if (getSolution(s.at(v)) < 0.5)
          input.deactivate(v, turnedOff);
      // Try propagation to turned off vertices
      input.propagate_to(turnedOff, turnedOn);

      // Second activate vertices
      turnedOff.clear();
      for (auto v : boost::make_iterator_range(vertices(graph)))
        if (getSolution(s.at(v)) > 0.5) {
          std::vector<Vertex> neighbors;
          for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
            if (getSolution(w.at(std::make_pair(v, u))) > 0.5)
              neighbors.push_back(u);
          input.activate(v, neighbors, turnedOn, turnedOff);
        }
      // Try propagation to turned off vertices
      input.propagate_to(turnedOff, turnedOn);
      // Try propagation from turned on vertices or their neighbors
      std::list<Vertex> candidates; 
      for (auto u: turnedOn) { 
        if (input.isZeroInjection(u))
          candidates.push_back(u);
        for (auto y: boost::make_iterator_range(adjacent_vertices(u, graph)))
          if (input.isZeroInjection(y) && input.isMonitored(y))
            candidates.push_back(y);
      }
      input.propagate_from(candidates, turnedOn);

      // Feasibility check
      //assert(input.check_get_monitored_set(sValue, wValue));
      
      if (!input.isFeasible()) {

        // Find violated cycles
        std::set<Fort> forts = violatedForts(lazyLimit);
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

  std::set<Fort> violatedForts(size_t fortsLimit) {

    std::set<Fort> forts;

    // Copy solution
    Pds newSolution (input);

    // Get unmonitored set
    std::vector<Vertex> unmonitoredSet;
    newSolution.get_unmonitored_set(unmonitoredSet);

    // Shuffle unmonitored set
    boost::range::random_shuffle(unmonitoredSet);

    // Preselect random neighbors for unmonitored vertices
    std::map<Vertex, std::vector<Vertex>> neighbors; 
    for (Vertex v : unmonitoredSet) {
      neighbors[v] = std::vector<Vertex> ();
      boost::copy(adjacent_vertices(v, graph) |
                      boost::adaptors::filtered([unmonitoredSet](auto v) {
                        return boost::range::find(unmonitoredSet, v) !=
                               unmonitoredSet.end();
                      }),
                  std::back_inserter(neighbors[v]));
      size_t k = std::min(neighbors[v].size(), n_channels - 1);
      if (k < neighbors[v].size()) {
        neighbors[v].resize(k);
        random_unique(neighbors[v].begin(), neighbors[v].end(), k);
      }
    }

    // Initialize components
    /*
    observers = std::vector<std::set<Vertex>> (num_vertices(graph), std::set<Vertex> ());
    propagates.clear();
    propagator.clear();
    digraph.clear();
    for (auto v : boost::make_iterator_range(vertices(graph)))
      add_vertex(LabelledVertex{.label = v}, digraph);
    */

    // Copy of monitored set
    // VertexList mS (monitoredSet);

    // Activate all unmonitored vertices
    // (propagation is unnecessary here)
    std::list<Vertex> trash;
    for (Vertex v : unmonitoredSet)
      newSolution.activate(v, neighbors[v], trash, trash);

    // MINIMISE FEASIBLE SOLUTION
    for (Vertex v : unmonitoredSet) {

      if (forts.size() >= fortsLimit)
        break;

      // Deactivate v
      std::list<Vertex> turnedOff;
      newSolution.deactivate(v, turnedOff);

      // Try propagations to changed vertices
      newSolution.propagate_to(turnedOff, trash);

      // Feasibility check
      //VertexList mS2 = input.get_monitored_set(sValue, wValue);
      //std::cout << fmt::format("mS: {},\t mS2: {}\n", mS, mS2);
      //assert(mS == mS2);
      if (!newSolution.isFeasible()) {
        
        // Find and insert fort
        forts.insert(findFort(newSolution));

        // Reactivate v
        std::list<Vertex> turnedOn;
        newSolution.activate(v, neighbors[v], turnedOn, trash);

        // Try propagations from changed vertices or their neighbors
        std::list<Vertex> candidates; 
        for (auto u: turnedOn) { 
          if (newSolution.isZeroInjection(u))
            candidates.push_back(u);
          for (auto y: boost::make_iterator_range(adjacent_vertices(u, graph)))
            if (newSolution.isZeroInjection(y) && newSolution.isMonitored(y))
              candidates.push_back(y);
        }
        newSolution.propagate_from(candidates, trash);

      }

    }

    return forts;
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

  Fort findFort(Pds &input) {
    Fort fort;
    for (auto u : boost::make_iterator_range(vertices(graph))) {
      if (input.isMonitored(u))
        continue;
      fort.first.insert(u);
      for (auto z : boost::make_iterator_range(adjacent_vertices(u, graph))) {
        if (!input.isMonitored(z))
          continue;
        if (!input.isActivated(z))
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
