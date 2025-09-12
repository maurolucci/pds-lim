#include "fort_solve.hpp"

#include <gurobi_c++.h>

#include <boost/range/adaptor/filtered.hpp>
#include <chrono>

#include "gurobi_common.hpp"

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
      : mipmodel(),
        model(*mipmodel.model),
        s(mipmodel.s),
        w(mipmodel.w),
        y(),
        input(input),
        graph(input.get_graph()),
        cbFile(callbackFile),
        solFile(solutionFile),
        n_channels(input.get_n_channels()),
        lazyLimit(lzLimit),
        totalCallback(mipmodel.totalCallback),
        totalCallbackTime(mipmodel.totalCallbackTime),
        totalLazy(mipmodel.totalLazy) {
    model.setCallback(this);

    // Add variables
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      s.try_emplace(
          v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
      if (degree(v, graph) > n_channels) {
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          w.try_emplace(std::make_pair(v, u),
                        model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     fmt::format("w_{}_{}", v, u)));
      }
    }

    // Add constraints
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      // (2) sum_{u \in N(v)} w_v_u <= omega_v * s_v, \forall v \in V_2
      if (degree(v, graph) > n_channels) {
        GRBLinExpr constr = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr += w.at(std::make_pair(v, u));
        model.addConstr(constr == n_channels * s.at(v));
      }
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

        auto t0 = std::chrono::high_resolution_clock::now();
        totalCallback++;

        // Update solution

        for (auto v : boost::make_iterator_range(vertices(graph))) {
          if (getSolution(s.at(v)) < 0.5)
            input.deactivate(v);
          else {
            std::vector<bool> dominate(degree(v, graph), false);
            size_t i = 0;
            for (auto u :
                 boost::make_iterator_range(adjacent_vertices(v, graph))) {
              if (degree(v, graph) <= input.get_n_channels() ||
                  getSolution(w.at(std::make_pair(v, u))) > 0.5)
                dominate[i] = true;
              ++i;
            }
            input.activate(v, dominate);
          }
        }

        // Feasibility check
        if (!input.isFeasible()) {
          // Find violated forts
          std::set<Fort> forts = violatedForts(lazyLimit);
          std::pair<double, double> avg = addLazyForts(forts);
          totalLazy += forts.size();

          // Report to callback file
          cbFile << fmt::format(
                        "# forts: {}, avg. vertex size: {:.2f}, avg. "
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

    // Get blank vertices
    // Unmonitored or unactivated?
    std::vector<Vertex> blank;
    // newSolution.get_unmonitored_set(blank);
    newSolution.get_unactivated_set(blank);

    // Shuffle blank set
    boost::range::random_shuffle(blank);

    // Activate blank vertices
    // Choose random neighbors to dominate
    std::map<Vertex, std::vector<bool>> dominate;
    for (Vertex v : blank) {
      dominate.emplace(v, std::vector<bool>(degree(v, graph), false));
      std::vector<size_t> indices(degree(v, graph));
      std::iota(indices.begin(), indices.end(), 0);
      size_t k = std::min(degree(v, graph), n_channels);
      if (k < degree(v, graph)) {
        for (size_t len = degree(v, graph), i = 0; i < k; ++i, --len) {
          int r = rand() % len;
          int j = indices[r];
          indices[r] = indices[i];
          indices[i] = j;
          dominate[v][indices[i]] = true;
        }
      }
      newSolution.activate(v, dominate[v]);
    }

    // Desactivate some vertices
    for (Vertex v : blank) {
      if (forts.size() >= fortsLimit) break;

      // Deactivate v
      newSolution.deactivate(v);

      // Feasibility check
      if (!newSolution.isFeasible()) {
        // Find and insert fort
        forts.insert(findFort(newSolution));

        // Reactivate v
        newSolution.activate(v, dominate[v]);
      }
    }

    return forts;
  }

  Fort findFort(Pds &input) {
    Fort fort;
    for (auto u : boost::make_iterator_range(vertices(graph))) {
      if (input.isMonitored(u)) continue;
      std::get<0>(fort).insert(u);
      for (auto z : boost::make_iterator_range(adjacent_vertices(u, graph))) {
        if (!input.isMonitored(z)) continue;
        if (degree(z, graph) <= input.get_n_channels())
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
}  // namespace

SolveResult solveLazyForts(Pds &input, boost::optional<std::string> logPath,
                           std::ostream &callbackFile, std::ostream &solFile,
                           double timeLimit, size_t lazyLimit) {
  LazyFortCB lazyForts(input, callbackFile, solFile, lazyLimit);
  return lazyForts.solve(logPath, timeLimit);
}

}  // end of namespace pds
