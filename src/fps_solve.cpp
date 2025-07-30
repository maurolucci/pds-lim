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
  // Map from precedences to propagation
  std::map<Edge, Edge> prec2props;
  std::ostream &cbFile, &solFile;
  size_t valIneq;
  size_t lazyLimit;

  size_t &totalCallback;
  size_t &totalCallbackTime;
  size_t &totalLazy;

  LazyFpsCB(Pds &input, std::ostream &callbackFile, std::ostream &solutionFile,
            size_t valIneq, size_t lzLimit)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()), cbFile(callbackFile),
        solFile(solutionFile), valIneq(valIneq), lazyLimit(lzLimit),
        totalCallback(mipmodel.totalCallback),
        totalCallbackTime(mipmodel.totalCallbackTime),
        totalLazy(mipmodel.totalLazy) {

    size_t n_channels = input.get_n_channels();

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
      if (input.isZeroInjection(v)) {
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          y.try_emplace(std::make_pair(v, u),
                        model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     fmt::format("y_{}_{}", v, u)));
      }
    }

    // Add constraints
    for (auto v : boost::make_iterator_range(vertices(graph))) {

      // (1) s_v + sum_{u \in N(v) \cap V_1} s_u + sum_{u \in N(v) \cap V_2}
      // w_u_v + sum_{u \in N(v) \cap V_Z} y_u_v == 1, \forall v \in V
      GRBLinExpr constr1 = 0;
      constr1 += s.at(v);
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        if (degree(u, graph) <= n_channels - 1)
          constr1 += s.at(u);
        else
          constr1 += w.at(std::make_pair(u, v));
        if (input.isZeroInjection(u))
          constr1 += y.at(std::make_pair(u, v));
      }
      model.addConstr(constr1 >= 1);

      // (3) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V_2
      if (degree(v, graph) > n_channels - 1) {
        GRBLinExpr constr3 = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr3 += w.at(std::make_pair(v, u));
        model.addConstr(constr3 ==
                        std::min(degree(v, graph), (n_channels - 1)) * s.at(v));
      }

      // Limitation of outgoing propagations
      // (4.1) sum_{u \in N(v)} y_{vu} <= 1 - s_v, \forall v \in V_Z \cap V_1
      // (4.1) sum_{u \in N(v)} y_{vu} <= 1, \forall v \in V_Z \cap V_2
      if ((valIneq == 1 || valIneq == 3) && input.isZeroInjection(v)) {
        GRBLinExpr constr4 = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr4 += y.at(std::make_pair(v, u));
        if (degree(v, graph) <= n_channels - 1)
          constr4 += s.at(v);
        model.addConstr(constr4 <= 1);
      }

      // Limitation of incomming propagations
      // (5.1) sum_{u \in N(v) \cap V_Z} y_{uv} <= 1, \forall v \in V
      if (valIneq == 2 || valIneq == 3) {
        GRBLinExpr constr5 = 0;
        size_t nlhs = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          if (input.isZeroInjection(u)) {
            constr5 += y.at(std::make_pair(u, v));
            nlhs++;
          }
        if (nlhs > 0)
          model.addConstr(constr5 <= 1);
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
    // Integer solution found (it does not necessarily improve the incumbent)
    case GRB_CB_MIPSOL:

      auto t0 = std::chrono::high_resolution_clock::now();
      totalCallback++;

      // Update solution

      // First deactivate vertices that became turned off
      std::list<Vertex> turnedOff, turnedOn;
      for (auto v : boost::make_iterator_range(vertices(graph)))
        if (getSolution(s.at(v)) < 0.5)
          input.deactivate(v, turnedOff);
      // Try propagation to turned off vertices
      input.propagate_to(turnedOff, turnedOn);

      // Second activate vertices that became turned on
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

        // Find violated fpss
        PrecedenceDigraph digraph = build_precedence_digraph();
        std::set<EdgeList> fpss = find_fpss(digraph, lazyLimit);
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
  // Function to build the precedence digraph according to the current
  // propagations. ATTENTION: redundant propagations are ignored and only
  // the unmonitored vertices are considered
  PrecedenceDigraph build_precedence_digraph() {

    prec2props.clear();

    PrecedenceDigraph digraph;
    std::map<Vertex, Vertex> name;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (input.isMonitored(v))
        continue;
      if (!name.contains(v))
        name[v] = boost::add_vertex(LabelledVertex{.label = v}, digraph);
      bool redudant = false;
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        if (!input.isZeroInjection(u))
          continue;
        if (getSolution(y.at(std::make_pair(u, v))) < 0.5)
          continue;
        if (redudant)
          break;
        redudant = true;
        if (!input.isMonitored(u)) {
          if (!name.contains(u))
            name[u] = boost::add_vertex(LabelledVertex{.label = u}, digraph);
          boost::add_edge(name[u], name[v], digraph);
          prec2props[std::make_pair(u, v)] = std::make_pair(u, v);
        }
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(u, graph))) {
          if (w == v || input.isMonitored(w))
            continue;
          if (!name.contains(w))
            name[w] = boost::add_vertex(LabelledVertex{.label = w}, digraph);
          boost::add_edge(name[w], name[v], digraph);
          prec2props[std::make_pair(w, v)] = std::make_pair(u, v);
        }
      }
    }

    return digraph;
  }

  // Function that gets the nearest visited predecessor of v, if any
  std::pair<Node, bool> get_nearest_visited_predecessor(
      PrecedenceDigraph &digraph,
      std::map<Node, std::pair<Node, size_t>> &preceded_by, Node v) {
    Node pred;
    size_t max_level = 0;
    for (auto e : make_iterator_range(in_edges(v, digraph))) {
      Node u = source(e, digraph);
      if (!preceded_by.contains(u))
        continue;
      auto p = preceded_by.at(u);
      if (p.second <= max_level) {
        continue;
      }
      pred = u;
      max_level = p.second;
    }
    return std::make_pair(pred, max_level != 0);
  }

  // Function that gets the farthest visited successor of v different from u,
  // if any
  std::pair<Node, bool> get_farthest_visited_successor(
      PrecedenceDigraph &digraph,
      std::map<Node, std::pair<Node, size_t>> &preceded_by, Node v, Node u) {
    Node succ;
    size_t min_level = std::numeric_limits<size_t>::max();
    for (auto e : make_iterator_range(out_edges(v, digraph))) {
      Node w = target(e, digraph);
      if (!preceded_by.contains(w) || w == u)
        continue;
      auto p = preceded_by.at(w);
      if (p.second >= min_level) {
        continue;
      }
      succ = w;
      min_level = p.second;
    }
    return std::make_pair(succ,
                          min_level != std::numeric_limits<size_t>::max());
  }

  // Functions that cuts a cycle from vertex u to vertex v
  void cut_cycle(std::map<Node, std::pair<Node, size_t>> &preceded_by, Node u,
                 Node v) {

    Node remove = u;
    Node pred;

    while (remove != v) {
      pred = preceded_by[remove].first;
      preceded_by.erase(remove);
      remove = pred;
    }

    return;
  }

  // Functions that cuts a cycle from vertex u to vertex v
  // and makes u preceded by v
  void cut_and_join_cycle(std::map<Node, std::pair<Node, size_t>> &preceded_by,
                          Node u, Node v) {

    size_t level = preceded_by[u].second;
    cut_cycle(preceded_by, u, v);
    preceded_by[u] = std::make_pair(v, level);

    return;
  }

  // Function that takes a cycle and makes it chordless
  void
  make_cycle_chordless(PrecedenceDigraph &digraph,
                       std::map<Node, std::pair<Node, size_t>> &preceded_by,
                       Node v) {
    Node pred = preceded_by[v].first;
    Node succ = v;
    // Check if pred has already visited successors (different from succ).
    // Select the farthest successor, i.e. with the lowest level.
    while (pred != v) {
      auto [u, ok] =
          get_farthest_visited_successor(digraph, preceded_by, pred, succ);
      if (ok)
        cut_and_join_cycle(preceded_by, u, pred);
      succ = pred;
      pred = preceded_by[succ].first;
    }
  }

  // Function that finds a chordless cycle in the precedence digraph, starting
  // from vertex v (v is not guaranteed to be in the cycle).  We already know
  // that such a cycle exists and that every vertex has at least one
  // predecessor.
  VertexList find_chordless_cycle(PrecedenceDigraph &digraph, Node v) {
    // Map from vertex to predecessor and level (height)
    std::map<Node, std::pair<Node, size_t>> preceded_by;
    Node lastVertex = v;

    for (int i = 1; !preceded_by.contains(lastVertex); ++i) {
      assert(in_degree(lastVertex, digraph) > 0);

      // Check if lastVertex has already visited predecessors.
      // Select the nearest predecessor, i.e. with the highest level.
      // This strategy avoid chords in one direction, while the cycle must be
      // postprocessed to avoid chords in the other direction.
      auto [next, ok] =
          get_nearest_visited_predecessor(digraph, preceded_by, lastVertex);
      if (!ok) {
        // Choose a random in-neighbor
        auto e = std::next(in_edges(lastVertex, digraph).first,
                           rand() % in_degree(lastVertex, digraph));
        next = source(*e, digraph);
      }
      preceded_by.emplace(lastVertex, std::make_pair(next, i));
      lastVertex = next;
    }

    // Cut every vertex but the cycle
    cut_cycle(preceded_by, v, lastVertex);

    // Cut the cycle until getting a chordless cycle
    make_cycle_chordless(digraph, preceded_by, lastVertex);

    // Get cycle in terms of the original vertices
    VertexList cycle;
    auto u = lastVertex;
    do {
      cycle.push_back(digraph[u].label);
      u = preceded_by.at(u).first;
    } while (u != lastVertex);

    // Rotate the fps so the minium element is in the front
    boost::range::rotate(fps, boost::range::min_element(fps));

    return cycle;
  }

  // Function thats maps a cycle to an FPS
  EdgeList get_fps(VertexList &cycle) {
    EdgeList fps;
    for (auto it = cycle.rbegin(); it != cycle.rend();) {
      Vertex v = *it++;
      int u = it != cycle.rend() ? *it : cycle.back();
      auto e = prec2props.at(std::make_pair(v, u));
      fps.push_back(e);
    }
    return fps;
  }

  // Function that finds a set of minimal FPSs (associated to chordless graphs)
  std::set<EdgeList> find_fpss(PrecedenceDigraph &digraph, size_t fpssLimit) {
    std::set<EdgeList> fpss;
    // Iterate over the vertices
    for (auto v : boost::make_iterator_range(vertices(digraph))) {
      // Check if the maximum number of FPSs has been found
      if (fpss.size() >= fpssLimit)
        break;
      // Find a chordless cycle (its existence is already guaranteed)
      VertexList cycle = find_chordless_cycle(digraph, v);
      // Map cycle to FPS
      EdgeList fps = get_fps(cycle);
      fpss.insert(fps);
    }
    return fpss;
  }

  double addLazyFpss(std::set<EdgeList> &fpss) {

    size_t accumFps = 0;

    for (const EdgeList &fps : fpss) {
      accumFps += fps.size();
      GRBLinExpr pathSum;
      for (auto [u, v] : fps)
        pathSum += y.at(std::make_pair(u, v));
      addLazy(pathSum <= fps.size() - 1);
    }

    return static_cast<double>(accumFps) / fpss.size();
  }
};
} // end of namespace

SolveResult solveLazyFpss(Pds &input, boost::optional<std::string> logPath,
                          std::ostream &callbackFile, std::ostream &solFile,
                          double timeLimit, size_t valIneq, size_t lazyLimit) {
  LazyFpsCB lazyFpss(input, callbackFile, solFile, valIneq, lazyLimit);
  return lazyFpss.solve(logPath, timeLimit);
}

} // end of namespace pds
