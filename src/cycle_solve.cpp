#include "gurobi_common.hpp"

#include <gurobi_c++.h>

namespace pds {

namespace {
struct LazyCycleCB : public GRBCallback {

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

  LazyCycleCB(Pds &input)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()) {

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
      model.addConstr(constr3 <=
                      std::min(boost::degree(v, graph), (n_channels - 1)) *
                          s.at(v));
    }

    // Build translation function
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (!input.isZeroInjection(v))
        continue;
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        translate[std::make_pair(v, u)].push_back(std::make_pair(v, u));
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
          if (w == u)
            continue;
          translate[std::make_pair(w, u)].push_back(std::make_pair(v, u));
        }
      }
    }
  }

  SolveResult solve(boost::optional<std::string> outPath, double timeLimit) {
    model.set(GRB_IntParam_LazyConstraints, 1);
    return solveMIP(input, mipmodel, outPath, timeLimit);
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
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          wValue[std::make_pair(v, u)] =
              getSolution(w.at(std::make_pair(v, u)));
      }

      // Feasibility check
      VertexList mS = input.get_monitored_set(sValue, wValue);
      if (!input.isFeasible(mS)) {

        // Find violated cycles
        PrecedenceDigraph digraph = build_precedence_digraph(mS);
        std::list<VertexList> cycles = violatedCycles(digraph, 100);
        addLazyCycles(cycles);
      }
      break;
    }
  }

private:
  PrecedenceDigraph build_precedence_digraph(VertexList &isMonitored) {

    PrecedenceDigraph digraph;
    std::map<Vertex, Vertex> name;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (isMonitored[v])
        continue;
      if (!name.contains(v))
        name[v] = boost::add_vertex(LabelledVertex{.label = v}, digraph);
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        if (!input.isZeroInjection(u))
          continue;
        if (getSolution(y.at(std::make_pair(u, v))) < 0.5)
          continue;
        if (!isMonitored[u]) {
          if (!name.contains(u))
            name[u] = boost::add_vertex(LabelledVertex{.label = u}, digraph);
          boost::add_edge(name[u], name[v], digraph);
        }
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(u, graph))) {
          if (w == v || isMonitored[w])
            continue;
          if (!name.contains(w))
            name[w] = boost::add_vertex(LabelledVertex{.label = w}, digraph);
          boost::add_edge(name[w], name[v], digraph);
        }
      }
    }

    return digraph;
  }

  std::list<VertexList> violatedCycles(PrecedenceDigraph &digraph,
                                       size_t cyclesLimit) {
    std::list<VertexList> cycles;
    // auto randomVertices = boost::range::random_shuffle(vertices(digraph));

    for (auto v : boost::make_iterator_range(vertices(digraph))) {
      if (cycles.size() >= cyclesLimit)
        break;
      VertexList cycle = findCycle(digraph, v);
      if (cycle.empty())
        continue;
      cycles.push_back(cycle);
    }
    return cycles;
  }

  VertexList findCycle(PrecedenceDigraph &digraph, Node v) {
    VertexList cycle;
    std::map<Node, std::pair<Node, size_t>> precededBy;
    Node lastVertex = v;

    for (int i = 1; !precededBy.contains(lastVertex); ++i) {
      if (in_degree(lastVertex, digraph) == 0)
        return cycle; // Return empty cycle

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
      cycle.push_back(digraph[u].label);
      u = precededBy.at(u).first;
    } while (u != lastVertex);

    return cycle;
  }

  void addLazyCycles(std::list<VertexList> &cycles) {

    for (VertexList &cycle : cycles) {
      EdgeList translated;
      for (auto it = cycle.rbegin(); it != cycle.rend();) {
        Vertex v = *it++;
        int u = it != cycle.rend() ? *it : cycle.back();
        for (auto e : translate.at(std::make_pair(v, u)))
          translated.push_back(e);
      }
      GRBLinExpr pathSum;
      for (auto [u, v] : translated)
        pathSum += y.at(std::make_pair(u, v));
      addLazy(pathSum <= cycle.size() - 1);
    }
  }
};
} // end of namespace

SolveResult solveLazyCycles(Pds &input, boost::optional<std::string> outPath,
                            double timeLimit) {
  LazyCycleCB lazyCycles(input);
  return lazyCycles.solve(outPath, timeLimit);
}

} // end of namespace pds
