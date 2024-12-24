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
    std::vector<std::set<Vertex>> observedBy(num_vertices(graph),
                                             std::set<Vertex>());
    for (Vertex v : unmonitoredSet) {
      sValue.at(v) = 1.0;
      observedBy[v].insert(v);
      // Apply random domination rules
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
      for (size_t i = 0; i < k; ++i) {
        wValue.at(std::make_pair(v, neighbors[i])) = 1.0;
        observedBy[neighbors[i]].insert(v);
      }
    }
    VertexList mS(num_vertices(graph), true);

    // Minimize feasible solution
    std::vector<std::set<Vertex>> inNeighbors(num_vertices(graph),
                                              std::set<Vertex>());
    std::vector<std::set<Vertex>> outNeighbors(num_vertices(graph),
                                               std::set<Vertex>());
    std::map<Vertex, Vertex> propagates;
    std::map<Vertex, Vertex> propagatedBy;
    for (Vertex v : unmonitoredSet) {

      if (forts.size() >= fortsLimit)
        break;

      sValue.at(v) = 0.0;

      std::list<Vertex> changes;
      observedBy[v].erase(v);
      changes.push_back(v);
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        if (wValue.at(std::make_pair(v, u)) > 0.5) {
          observedBy[u].erase(v);
          changes.push_back(u);
        }

      std::set<Vertex> needsPropagation;
      while (!changes.empty()) {
        Vertex u = changes.front();
        changes.pop_front();
        if (observedBy[u].empty()) {
          // u is no longer monitored
          mS[u] = false;
          needsPropagation.insert(u);
          // u cannot propagate anymore
          if (propagates.contains(u)) {
            Vertex x = propagates[u];
            propagates.erase(u);
            propagatedBy.erase(x);
            mS[x] = false;
            needsPropagation.insert(x);
            inNeighbors[x].erase(u);
            outNeighbors[u].erase(x);
            for (auto y :
                 boost::make_iterator_range(adjacent_vertices(u, graph))) {
              if (y == x)
                continue;
              inNeighbors[x].erase(y);
              outNeighbors[y].erase(x);
            }
            changes.push_back(x);
          }
          // u cannot join a propagation
          std::set<Vertex> uOut(outNeighbors[u]);
          for (Vertex x : uOut) {
            // x can no longer be propagated
            Vertex y = propagatedBy[x];
            propagates.erase(y);
            propagatedBy.erase(x);
            mS[x] = false;
            needsPropagation.insert(x);
            inNeighbors[x].erase(y);
            outNeighbors[y].erase(x);
            for (auto z :
                 boost::make_iterator_range(adjacent_vertices(y, graph))) {
              if (z == x)
                continue;
              inNeighbors[x].erase(z);
              outNeighbors[z].erase(x);
            }
            changes.push_back(x);
          }
        }
      }

      bool keepGoing = true;
      while (keepGoing) {
        keepGoing = false;
        std::set<Vertex> propagated;
        for (auto u : needsPropagation) {
          for (auto x :
               boost::make_iterator_range(adjacent_vertices(u, graph))) {
            if (!mS[x] || !input.isZeroInjection(x))
              continue;
            size_t count = boost::range::count_if(
                boost::adjacent_vertices(x, graph),
                [mS, u](auto y) { return y != u && mS[y]; });
            if (boost::degree(x, graph) - count != 1)
              continue;
            mS[u] = true;
            propagates[x] = u;
            propagatedBy[u] = x;
            inNeighbors[u].insert(x);
            outNeighbors[x].insert(u);
            for (auto y :
                 boost::make_iterator_range(adjacent_vertices(x, graph))) {
              if (y == u)
                continue;
              inNeighbors[u].insert(y);
              outNeighbors[y].insert(u);
            }
            keepGoing = true;
            propagated.insert(u);
            break;
          }
        }
        for (auto u : propagated)
          needsPropagation.erase(u);
      }

      VertexList mS2 = input.get_monitored_set(sValue, wValue);
      assert(mS == mS2);
      if (!input.isFeasible(mS2)) {
        forts.insert(findFort(mS2));
        sValue.at(v) = 1.0;

        // Acomodate everything
        std::list<Vertex> changes;
        // v ahora se monitorea solo
        if (observedBy[v].empty()) {
          changes.push_back(v);
          // Ya no necesita que lo propaguen
          if (propagatedBy.contains(v)) {
            Vertex u = propagatedBy[v];
            propagatedBy.erase(v);
            propagates.erase(u);
            inNeighbors[v].erase(u);
            outNeighbors[u].erase(v);
            for (auto y :
                 boost::make_iterator_range(adjacent_vertices(u, graph))) {
              if (y == v)
                continue;
              inNeighbors[v].erase(y);
              outNeighbors[y].erase(v);
            }
          }
        }
        observedBy[v].insert(v);
        mS[v] = true;
        // v monitorea a sus vecinos
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
          if (wValue.at(std::make_pair(v, u)) > 0.5) {
            if (observedBy[u].empty()) {
              changes.push_back(u);
              // u ya no necesita que lo propaguen
              if (propagatedBy.contains(u)) {
                Vertex y = propagatedBy[u];
                propagatedBy.erase(u);
                propagates.erase(y);
                inNeighbors[u].erase(y);
                outNeighbors[y].erase(u);
                for (auto z :
                     boost::make_iterator_range(adjacent_vertices(y, graph))) {
                  if (z == u)
                    continue;
                  inNeighbors[u].erase(z);
                  outNeighbors[z].erase(u);
                }
              }
            }
            observedBy[u].insert(v);
            mS[u] = true;
            for (auto z :
                 boost::make_iterator_range(adjacent_vertices(u, graph)))
              changes.push_back(z);
          }
          changes.push_back(u);
        }

        // Propagate!
        while (!changes.empty()) {
          auto u = changes.front();
          changes.pop_front();
          // u can propagate?
          if (input.isZeroInjection(u)) {
            size_t count =
                boost::range::count_if(boost::adjacent_vertices(u, graph),
                                       [mS](auto y) { return mS[y]; });
            if (boost::degree(u, graph) - count != 1)
              continue;
            auto it_y =
                boost::range::find_if(boost::adjacent_vertices(u, graph),
                                      [mS](auto y) { return !mS[y]; });
            mS[*it_y] = true;
            propagates[u] = *it_y;
            propagatedBy[*it_y] = u;
            inNeighbors[*it_y].insert(u);
            outNeighbors[u].insert(*it_y);
            changes.push_back(*it_y);
            for (auto z :
                 boost::make_iterator_range(adjacent_vertices(*it_y, graph)))
              changes.push_back(z);
            for (auto z :
                 boost::make_iterator_range(adjacent_vertices(u, graph))) {
              if (z == *it_y)
                continue;
              inNeighbors[*it_y].insert(z);
              outNeighbors[z].insert(*it_y);
            }
          }
        }
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