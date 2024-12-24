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
        lazyLimit(lzLimit) {

    model.setCallback(this);

    // Add variables
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      s.try_emplace(
          v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        w.try_emplace(std::make_pair(v, u),
                      model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                   fmt::format("w_{}_{}", v, u)));
      }
    }

    // (2) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      GRBLinExpr constr3 = 0;
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        constr3 += w.at(std::make_pair(v, u));
      }
      model.addConstr(constr3 <=
                      std::min(boost::degree(v, graph), (n_channels - 1)) *
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
  std::set<Fort> violatedForts(const VertexList &monitoredSet,
                               size_t fortsLimit) {

    std::set<Fort> forts;

    // Shuffle unmonitored vertices
    VertexList unmonitoredSet;
    boost::copy(vertices(graph) |
                    boost::adaptors::filtered(
                        [monitoredSet](auto v) { return !monitoredSet[v]; }),
                std::back_inserter(unmonitoredSet));
    boost::range::random_shuffle(unmonitoredSet);

    // Complete feasible solution
    for (Vertex v : unmonitoredSet) {
      sValue.at(v) = 1.0;
      // Apply random domination rules
      size_t k = std::min(degree(v, graph), n_channels - 1);
      VertexList neighbors(degree(v, graph));
      boost::copy(adjacent_vertices(v, graph), neighbors.begin());
      if (k < degree(v, graph))
        random_unique(neighbors.begin(), neighbors.end(), k);
      for (size_t i = 0; i < k; ++i)
        wValue.at(std::make_pair(v, neighbors[i])) = 1.0;
    }

    // Minimize feasible solution
    for (Vertex v : unmonitoredSet) {

      if (forts.size() >= fortsLimit)
        break;

      sValue.at(v) = 0.0;
      VertexList mS = input.get_monitored_set(sValue, wValue);
      if (!input.isFeasible(mS)) {
        forts.insert(findFort(mS));
        sValue.at(v) = 1.0;
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