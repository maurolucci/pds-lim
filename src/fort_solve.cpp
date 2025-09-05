#include "fort_solve.hpp"
#include "gurobi_common.hpp"

#include <boost/range/adaptor/filtered.hpp>
#include <chrono>
#include <gurobi_c++.h>

namespace pds {

using Fort = std::tuple<std::set<Vertex>, std::set<Vertex>, std::set<Edge>>;

namespace {

struct LazyFortCB : public GRBCallback {

  MIPModel mipmodel;
  GRBModel &model;
  std::map<pds::Vertex, GRBVar> &s;
  std::map<Edge, GRBVar> &w;
  std::map<Edge, GRBVar> y;
  Pds &input;
  const PowerGrid &graph;
  std::map<Edge, EdgeList> translate;
  std::ostream &cbFile, &solFile;
  size_t n_channels;
  size_t lazyLimit;

  size_t &totalCallback;
  size_t &totalCallbackTime;
  size_t &totalLazy;

  LazyFortCB(Pds &input, std::ostream &callbackFile, std::ostream &solutionFile,
             size_t lzLimit)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()), cbFile(callbackFile),
        solFile(solutionFile), n_channels(input.get_n_channels()),
        lazyLimit(lzLimit), totalCallback(mipmodel.totalCallback),
        totalCallbackTime(mipmodel.totalCallbackTime),
        totalLazy(mipmodel.totalLazy) {

    model.setCallback(this);

    // Add variables
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      s.try_emplace(
          v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
      if (degree(v, graph) > n_channels - 1) {
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          w.try_emplace(std::make_pair(v, u),
                        model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     fmt::format("w_{}_{}", v, u)));
      }
    }

    // Add constraints
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      // (2) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V_2
      if (degree(v, graph) > n_channels - 1) {
        GRBLinExpr constr = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr += w.at(std::make_pair(v, u));
        model.addConstr(constr ==
                        std::min(degree(v, graph), (n_channels - 1)) * s.at(v));
      }
    }

    // Turn-off presolve
    model.set(GRB_IntParam_MinRelNodes, 0);
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

      auto t0 = std::chrono::high_resolution_clock::now();
      totalCallback++;

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
            if (degree(v, graph) <= input.get_n_channels() - 1 ||
                getSolution(w.at(std::make_pair(v, u))) > 0.5)
              neighbors.push_back(u);
          input.activate(v, neighbors, turnedOn, turnedOff);
        }
      // Try propagation to turned off vertices
      input.propagate_to(turnedOff, turnedOn);
      // Try propagation from turned on vertices or their neighbors
      std::list<Vertex> candidates;
      for (auto u : turnedOn) {
        if (input.isZeroInjection(u))
          candidates.push_back(u);
        for (auto y : boost::make_iterator_range(adjacent_vertices(u, graph)))
          if (input.isZeroInjection(y) && input.isMonitored(y))
            candidates.push_back(y);
      }
      input.propagate_from(candidates, turnedOn);

      // Feasibility check
      if (!input.isFeasible()) {

        // Find violated forts
        std::set<Fort> forts = violatedForts(lazyLimit);
        std::pair<double, double> avg = addLazyForts(forts);
        totalLazy += forts.size();

        // Report to callback file
        cbFile << fmt::format("# forts: {}, avg. vertex size: {:.2f}, avg. "
                              "edge size: {:.2f}",
                              forts.size(), avg.first, avg.second)
               << std::endl;
      }

      auto t1 = std::chrono::high_resolution_clock::now();
      totalCallbackTime +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();

      break;
    }
  }

private:
  std::set<Fort> violatedForts(size_t fortsLimit) {

    std::set<Fort> forts;

    // Copy solution
    Pds newSolution(input);

    // Get unactivated set
    std::vector<Vertex> unActivatedSet;
    newSolution.get_unactivated_set(unActivatedSet);

    // Shuffle unactivated set
    boost::range::random_shuffle(unActivatedSet);

    // Preselect random neighbors for unactivated vertices
    std::map<Vertex, std::vector<Vertex>> neighbors;
    for (Vertex v : unActivatedSet) {
      neighbors[v] = std::vector<Vertex>();
      boost::copy(adjacent_vertices(v, graph) |
                      boost::adaptors::filtered([unActivatedSet](auto v) {
                        return boost::range::find(unActivatedSet, v) !=
                               unActivatedSet.end();
                      }),
                  std::back_inserter(neighbors[v]));
      size_t k = std::min(neighbors[v].size(), n_channels - 1);
      if (k < neighbors[v].size()) {
        neighbors[v].resize(k);
        random_unique(neighbors[v].begin(), neighbors[v].end(), k);
      }
    }

    // Activate all unactivated vertices
    // (propagation is unnecessary here)
    std::list<Vertex> trash;
    for (Vertex v : unActivatedSet)
      newSolution.activate(v, neighbors[v], trash, trash);

    // MINIMISE FEASIBLE SOLUTION
    for (Vertex v : unActivatedSet) {

      if (forts.size() >= fortsLimit)
        break;

      // Deactivate v
      std::list<Vertex> turnedOff;
      newSolution.deactivate(v, turnedOff);

      // Try propagations to changed vertices
      newSolution.propagate_to(turnedOff, trash);

      // Feasibility check
      if (!newSolution.isFeasible()) {

        // Find and insert fort
        forts.insert(findFort(newSolution));

        // Reactivate v
        std::list<Vertex> turnedOn;
        newSolution.activate(v, neighbors[v], turnedOn, trash);

        // Try propagations from changed vertices or their neighbors
        std::list<Vertex> candidates;
        for (auto u : turnedOn) {
          if (newSolution.isZeroInjection(u))
            candidates.push_back(u);
          for (auto y : boost::make_iterator_range(adjacent_vertices(u, graph)))
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
      std::get<0>(fort).insert(u);
      for (auto z : boost::make_iterator_range(adjacent_vertices(u, graph))) {
        if (!input.isMonitored(z))
          continue;
        if (degree(z, graph) <= input.get_n_channels() - 1)
          std::get<1>(fort).insert(z);
        else
          std::get<2>(fort).insert(std::make_pair(z, u));
      }
    }
    return fort;
  }

  std::pair<double, double> addLazyForts(std::set<Fort> forts) {

    size_t accumVertices = 0;
    size_t accumEdges = 0;

    for (auto &f : forts) {
      GRBLinExpr fortSum;
      for (auto v : std::get<0>(f)) {
        fortSum += s.at(v);
        accumVertices++;
      }
      for (auto v : std::get<1>(f)) {
        fortSum += s.at(v);
        accumVertices++;
      }
      for (auto e : std::get<2>(f)) {
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
