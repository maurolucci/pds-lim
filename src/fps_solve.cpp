#include "fps_solve.hpp"
#include "gurobi_common.hpp"

#include <fstream>
#include <gurobi_c++.h>

namespace pds {

namespace {
struct LazyFpsCB : public GRBCallback {

  MIPModel mipmodel;
  GRBModel &model;
  std::map<pds::Vertex, GRBVar> &s;
  std::map<Edge, GRBVar> &w;
  std::map<Edge, GRBVar> y;
  Pds &input;
  const PowerGrid &graph;
  std::map<Edge, EdgeList> translate;
  std::ostream &cbFile, &solFile;
  size_t lazyLimit;

  size_t &totalCallback;
  size_t &totalCallbackTime;
  size_t &totalLazy;

  LazyFpsCB(Pds &input, std::ostream &callbackFile, std::ostream &solutionFile,
            size_t lzLimit)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()), cbFile(callbackFile),
        solFile(solutionFile), lazyLimit(lzLimit),
        totalCallback(mipmodel.totalCallback),
        totalCallbackTime(mipmodel.totalCallbackTime),
        totalLazy(mipmodel.totalLazy) {

    size_t n_channels = input.get_n_channels();

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
        if (input.isZeroInjection(v)) {
          y.try_emplace(std::make_pair(v, u),
                        model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     fmt::format("y_{}_{}", v, u)));
        }
      }
    }

    // Add constraints
    for (auto v : boost::make_iterator_range(vertices(graph))) {

      // (1) s_v + sum_{u \in N(v)} w_u_v + sum_{u \in N(v) \cap V_Z} y_u_v ==
      // 1, \forall v \in V
      GRBLinExpr constr1 = 0;
      constr1 += s.at(v);
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        constr1 += w.at(std::make_pair(u, v));
        if (input.isZeroInjection(u))
          constr1 += y.at(std::make_pair(u, v));
      }
      model.addConstr(constr1 >= 1);

      // (3) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V
      GRBLinExpr constr3 = 0;
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        constr3 += w.at(std::make_pair(v, u));
      }
      model.addConstr(constr3 ==
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
            if (getSolution(w.at(std::make_pair(v, u))) > 0.5)
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

        // Find violated fpss
        PrecedenceDigraph digraph = build_precedence_digraph();
        std::set<VertexList> fpss = violatedFpss(digraph, lazyLimit);
        double avg = addLazyFpss(fpss);
        totalLazy += fpss.size();

        // Report to callback file
        cbFile << fmt::format("# fps: {}, avg. size: {:.2f}", fpss.size(), avg)
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
  PrecedenceDigraph build_precedence_digraph() {

    translate.clear();

    PrecedenceDigraph digraph;
    std::map<Vertex, Vertex> name;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (input.isMonitored(v))
        continue;
      if (!name.contains(v))
        name[v] = boost::add_vertex(LabelledVertex{.label = v}, digraph);
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        if (!input.isZeroInjection(u))
          continue;
        if (getSolution(y.at(std::make_pair(u, v))) < 0.5)
          continue;
        if (!input.isMonitored(u)) {
          if (!name.contains(u))
            name[u] = boost::add_vertex(LabelledVertex{.label = u}, digraph);
          boost::add_edge(name[u], name[v], digraph);
          translate[std::make_pair(u, v)].push_back(std::make_pair(u, v));
        }
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(u, graph))) {
          if (w == v || input.isMonitored(w))
            continue;
          if (!name.contains(w))
            name[w] = boost::add_vertex(LabelledVertex{.label = w}, digraph);
          boost::add_edge(name[w], name[v], digraph);
          translate[std::make_pair(w, v)].push_back(std::make_pair(u, v));
        }
      }
    }

    return digraph;
  }

  std::set<VertexList> violatedFpss(PrecedenceDigraph &digraph,
                                    size_t fpssLimit) {
    std::set<VertexList> fpss;
    // auto randomVertices = boost::range::random_shuffle(vertices(digraph));

    for (auto v : boost::make_iterator_range(vertices(digraph))) {
      if (fpss.size() >= fpssLimit)
        break;
      VertexList fps = findFps(digraph, v);
      if (fps.empty())
        continue;
      fpss.insert(fps);
    }
    return fpss;
  }

  VertexList findFps(PrecedenceDigraph &digraph, Node v) {
    VertexList fps;
    std::map<Node, std::pair<Node, size_t>> precededBy;
    Node lastVertex = v;

    for (int i = 1; !precededBy.contains(lastVertex); ++i) {
      if (in_degree(lastVertex, digraph) == 0)
        return fps; // Return empty Fps

      // Check if lastVertex has an already visited in-neighbor
      // Select the last visited in-neighbor
      Node next;
      size_t max_i = 0;
      for (auto e : make_iterator_range(in_edges(lastVertex, digraph))) {
        Node u = source(e, digraph);
        if (!precededBy.contains(u))
          continue;
        auto p = precededBy.at(u);
        if (p.second <= max_i) {
          continue;
        }
        next = u;
        max_i = p.second;
      }

      if (max_i == 0) {
        // Choose a random in-neighbor
        auto e = std::next(in_edges(lastVertex, digraph).first,
                           rand() % in_degree(lastVertex, digraph));
        next = source(*e, digraph);
      }
      precededBy.emplace(lastVertex, std::make_pair(next, i));
      lastVertex = next;
    }

    // Traverse de cycle
    auto u = lastVertex;
    do {
      fps.push_back(digraph[u].label);
      u = precededBy.at(u).first;
    } while (u != lastVertex);

    // Rotate the fps so the minium element is in the front
    boost::range::rotate(fps, boost::range::min_element(fps));

    return fps;
  }

  double addLazyFpss(std::set<VertexList> &fpss) {

    size_t accumFps = 0;

    for (const VertexList &fps : fpss) {
      EdgeList translated;
      for (auto it = fps.rbegin(); it != fps.rend();) {
        Vertex v = *it++;
        int u = it != fps.rend() ? *it : fps.back();
        auto e = std::next(translate.at(std::make_pair(v, u)).begin(),
                           rand() % translate.at(std::make_pair(v, u)).size());
        translated.push_back(*e);
        accumFps++;
      }
      GRBLinExpr pathSum;
      for (auto [u, v] : translated)
        pathSum += y.at(std::make_pair(u, v));
      addLazy(pathSum <= fps.size() - 1);
    }

    return static_cast<double>(accumFps) / fpss.size();
  }
};
} // end of namespace

SolveResult solveLazyFpss(Pds &input, boost::optional<std::string> logPath,
                          std::ostream &callbackFile, std::ostream &solFile,
                          double timeLimit, size_t lazyLimit) {
  LazyFpsCB lazyFpss(input, callbackFile, solFile, lazyLimit);
  return lazyFpss.solve(logPath, timeLimit);
}

} // end of namespace pds
