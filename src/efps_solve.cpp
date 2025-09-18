#include "efps_solve.hpp"

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <fstream>
#include <gurobi_c++.h>

#include "gurobi_common.hpp"

#define EPSILON 1e-6

namespace pds {

namespace {
struct LazyEfpsCB : public GRBCallback {
  MIPModel mipmodel;
  GRBModel &model;
  std::map<pds::Vertex, GRBVar> &s;
  std::map<Edge, GRBVar> &w;
  std::map<Edge, GRBVar> y;
  Pds &input;
  const PowerGrid &graph;
  // Map from precedence to list of propagations
  std::map<Edge, EdgeList> prec2props;
  std::ostream &cbFile, &solFile;
  size_t lazyLimit;
  size_t cutLimit;

  size_t &totalLazyCBCalls;
  size_t &totalLazyCBTime;
  size_t &totalLazyAdded;

  size_t &totalCutCBCalls;
  size_t &totalCutCBTime;
  size_t &totalCutAdded;

  LazyEfpsCB(Pds &input, std::ostream &callbackFile, std::ostream &solutionFile,
             bool inProp, bool outProp, bool initEFPS, size_t lzLimit,
             size_t cutLimit)
      : mipmodel(), model(*mipmodel.model), s(mipmodel.s), w(mipmodel.w), y(),
        input(input), graph(input.get_graph()), cbFile(callbackFile),
        solFile(solutionFile), lazyLimit(lzLimit), cutLimit(cutLimit),
        totalLazyCBCalls(mipmodel.totalLazyCBCalls),
        totalLazyCBTime(mipmodel.totalLazyCBTime),
        totalLazyAdded(mipmodel.totalLazyAdded),
        totalCutCBCalls(mipmodel.totalCutCBCalls),
        totalCutCBTime(mipmodel.totalCutCBTime),
        totalCutAdded(mipmodel.totalCutAdded) {
    size_t n_channels = input.get_n_channels();

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
        if (degree(u, graph) <= n_channels)
          constr1 += s.at(u);
        else
          constr1 += w.at(std::make_pair(u, v));
        if (input.isZeroInjection(u))
          constr1 += y.at(std::make_pair(u, v));
      }
      model.addConstr(constr1 >= 1);

      // (3) sum_{u \in N(v)} w_v_u <= omega_v * s_v, \forall v \in V_2
      if (degree(v, graph) > n_channels) {
        GRBLinExpr constr3 = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr3 += w.at(std::make_pair(v, u));
        model.addConstr(constr3 == n_channels * s.at(v));
      }

      // Limitation of outgoing propagations
      // (4.1) sum_{u \in N(v)} y_{vu} <= 1 - s_v, \forall v \in V_Z \cap V_1
      // (4.1) sum_{u \in N(v)} y_{vu} <= 1, \forall v \in V_Z \cap V_2
      if (outProp && input.isZeroInjection(v)) {
        GRBLinExpr constr4 = 0;
        for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
          constr4 += y.at(std::make_pair(v, u));
        if (degree(v, graph) <= n_channels)
          constr4 += s.at(v);
        model.addConstr(constr4 <= 1);
      }

      // Limitation of incomming propagations
      // (5.1) sum_{u \in N(v) \cap V_Z} y_{uv} <= 1, \forall v \in V
      if (inProp) {
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

    build_prec2props_map();

    // Initial FPS constraints associated to cycles of length 2
    if (initEFPS) {
      // Find cyclic precedences of length 2
      for (auto &[prec, props1] : prec2props) {
        auto [u, v] = prec;
        // Avoid symmetries
        if (u >= v)
          continue;
        // Check if there is an opposite precedence
        if (!prec2props.contains(std::make_pair(v, u)))
          continue;
        auto &props2 = prec2props[std::make_pair(v, u)];
        GRBLinExpr constr6;
        for (auto &p1 : props1)
          constr6 += y.at(p1);
        for (auto &p2 : props2)
          constr6 += y.at(p2);
        model.addConstr(constr6 <= 1);
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
    case GRB_CB_MIPSOL: {
      auto t0 = std::chrono::high_resolution_clock::now();
      totalLazyCBCalls++;

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
        // Find violated cycles
        PrecedenceDigraph digraph = build_precedence_digraph();
        std::set<std::pair<EdgeList, size_t>> efpss =
            find_efps_lazy_constraints(digraph);
        std::pair<double, double> avg = addLazyEfpss(efpss);
        totalLazyAdded += efpss.size();

        // Report to callback file
        cbFile
            << fmt::format(
                   "# lazy efps: {}, avg. size: {:.2f}, avg. ext. size: {:.2f}",
                   efpss.size(), avg.first, avg.second)
            << std::endl;
      }

      auto t1 = std::chrono::high_resolution_clock::now();
      totalLazyCBTime +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();

      break;
    }
    // Node callback
    case GRB_CB_MIPNODE: {
      if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {

        auto t0 = std::chrono::high_resolution_clock::now();
        totalCutCBCalls++;

        // Build weighted precedence digraph
        WeightedPrecedenceDigraph digraph = build_weighted_precedence_digraph();
        // Find violated cycles
        std::set<std::pair<EdgeList, size_t>> efpss = find_efps_cuts(digraph);
        if (efpss.size() == 0)
          break;
        std::pair<double, double> avg = addCutEfpss(efpss);
        totalCutAdded += efpss.size();

        // Report to callback file
        cbFile
            << fmt::format(
                   "# cut efps: {}, avg. size: {:.2f}, avg. ext. size: {:.2f}",
                   efpss.size(), avg.first, avg.second)
            << std::endl;

        // Time spent
        auto t1 = std::chrono::high_resolution_clock::now();
        totalCutCBTime +=
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
                .count();
      }

      break;
    }
    }
  }

private:
  // Function to build the map from a precedence to the list of propagations
  // that impose it
  void build_prec2props_map() {
    prec2props.clear();
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      if (!input.isZeroInjection(v))
        continue;
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        prec2props[std::make_pair(v, u)].push_back(std::make_pair(v, u));
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
          if (w == u)
            continue;
          prec2props[std::make_pair(w, u)].push_back(std::make_pair(v, u));
        }
      }
    }
  }

  // Function to build the precedence digraph according to the current
  // propagations. ATTENTION: redundant propagations are ignored and only
  // the unmonitored vertices are considered
  PrecedenceDigraph build_precedence_digraph() {
    PrecedenceDigraph digraph;
    std::map<Vertex, Node> name;
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
          if (!boost::edge(name[u], name[v], digraph).second)
            boost::add_edge(name[u], name[v], digraph);
        }
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(u, graph))) {
          if (w == v || input.isMonitored(w))
            continue;
          if (!name.contains(w))
            name[w] = boost::add_vertex(LabelledVertex{.label = w}, digraph);
          if (!boost::edge(name[w], name[v], digraph).second)
            boost::add_edge(name[w], name[v], digraph);
        }
      }
    }
    return digraph;
  }

  using EdgeWeightMap = boost::adj_list_edge_property_map<
      boost::bidirectional_tag, double, double &, pds::WeightedNode,
      boost::property<boost::edge_weight_t, double>, boost::edge_weight_t>;

  // Function to add an edge (u,v) to the weighted precedence digraph if not
  // already present, or update its weight otherwise
  void add_edge_if_needed(Vertex u, Vertex v, double weight,
                          WeightedPrecedenceDigraph &digraph,
                          std::map<Vertex, WeightedNode> &name,
                          EdgeWeightMap &weightMap) {
    if (!name.contains(v))
      name[v] = boost::add_vertex(LabelledVertex{.label = v}, digraph);
    if (!name.contains(u))
      name[u] = boost::add_vertex(LabelledVertex{.label = u}, digraph);
    if (!boost::edge(name[u], name[v], digraph).second) {
      auto e = boost::add_edge(name[u], name[v], digraph).first;
      boost::put(weightMap, e, std::max(0.0, 1.0 - weight));
    } else {
      auto e = boost::edge(name[u], name[v], digraph).first;
      boost::put(weightMap, e,
                 std::max(0.0, boost::get(weightMap, e) - weight));
    }
  }

  // Function to build the weighted precedence digraph according to the
  // current propagations
  WeightedPrecedenceDigraph build_weighted_precedence_digraph() {
    WeightedPrecedenceDigraph digraph;
    EdgeWeightMap weightMap = boost::get(boost::edge_weight, digraph);
    std::map<Vertex, WeightedNode> name;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        if (!input.isZeroInjection(u))
          continue;
        double weight = getNodeRel(y.at(std::make_pair(u, v)));
        if (weight < EPSILON)
          continue;
        add_edge_if_needed(u, v, weight, digraph, name, weightMap);
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(u, graph))) {
          if (w == v)
            continue;
          add_edge_if_needed(w, v, weight, digraph, name, weightMap);
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
    boost::range::rotate(cycle, boost::range::min_element(cycle));

    return cycle;
  }

  // Function thats maps a cycle to an EFPS
  EdgeList get_efps(VertexList &cycle) {
    EdgeList efps;
    for (auto it = cycle.rbegin(); it != cycle.rend();) {
      Vertex v = *it++;
      int u = it != cycle.rend() ? *it : cycle.back();
      for (auto &e : prec2props.at(std::make_pair(v, u)))
        efps.push_back(e);
    }
    return efps;
  }

  // Function that finds a set of EFPS constraints to be used as lazy
  // constraints. The size of the cycle is saved in the 2nd components.
  std::set<std::pair<EdgeList, size_t>>
  find_efps_lazy_constraints(PrecedenceDigraph &digraph) {
    std::set<std::pair<EdgeList, size_t>> efpss;
    // Iterate over the vertices
    for (auto v : boost::make_iterator_range(vertices(digraph))) {
      // Check if the maximum number of FPSs has been found
      if (efpss.size() >= lazyLimit)
        break;
      // Find a chordless cycle (its existence is already guaranteed)
      VertexList cycle = find_chordless_cycle(digraph, v);
      // Map cycle to EFPS
      efpss.emplace(std::make_pair(get_efps(cycle), cycle.size()));
    }
    return efpss;
  }

  std::pair<double, double>
  addLazyEfpss(std::set<std::pair<EdgeList, size_t>> &efpss) {
    size_t accumCycle = 0;
    size_t accumExt = 0;

    for (auto &[efps, size] : efpss) {
      accumCycle += size;
      accumExt += efps.size();
      GRBLinExpr pathSum;
      for (auto [u, v] : efps)
        pathSum += y.at(std::make_pair(u, v));
      addLazy(pathSum <= size - 1);
    }

    return std::make_pair(static_cast<double>(accumCycle) / efpss.size(),
                          static_cast<double>(accumExt) / efpss.size());
  }

  // Find the minimum weight cycle that contains the edge e
  // Returns a tuple with the cycle, its weight, and a boolean that is true
  // if a cycle has been found, false otherwise
  std::tuple<VertexList, double, bool>
  find_min_weight_cycle(WeightedPrecedenceDigraph &digraph, WeightedArc e) {

    // Find the shortest path from v to u
    auto v = target(e, digraph);
    auto u = source(e, digraph);
    std::vector<double> distances(num_vertices(digraph));
    std::vector<WeightedNode> predecessors(num_vertices(digraph));
    auto weight_map = boost::get(boost::edge_weight, digraph);
    dijkstra_shortest_paths(digraph, v,
                            boost::distance_map(&distances[0])
                                .predecessor_map(&predecessors[0])
                                .weight_map(weight_map));

    // If there is not path from v to u, return false
    if (distances[u] == std::numeric_limits<double>::max())
      return std::make_tuple(VertexList(), 0.0, false);

    // If the weight of the cycle is greater or equal than 1, return false
    if (distances[u] + boost::get(weight_map, e) >= 1.0 - EPSILON)
      return std::make_tuple(VertexList(), 0.0, false);

    // Reconstruct the shortest path from v to u
    VertexList path;
    for (auto x = u; x != v; x = predecessors[x])
      path.push_back(digraph[x].label);
    path.push_back(digraph[v].label);
    // Rotate the cycle so the minimum element is in the front
    boost::range::rotate(path, boost::range::min_element(path));

    return std::make_tuple(path, distances[u] + boost::get(weight_map, e),
                           true);
  }

  // Function that finds a set of EFPS constraints to be used as cuts.
  // The size of the cycle is saved in the 2nd components.
  std::set<std::pair<EdgeList, size_t>>
  find_efps_cuts(WeightedPrecedenceDigraph &digraph) {
    std::set<std::pair<EdgeList, size_t>> efpss;
    // Iterate through edges
    for (auto e : boost::make_iterator_range(edges(digraph))) {
      // Check if the maximum number of FPSs has been found
      if (efpss.size() >= cutLimit)
        break;
      // Find minimum weight cycle containing e
      auto [cycle, weight, ok] = find_min_weight_cycle(digraph, e);
      if (!ok)
        continue;
      // Map cycle to EFPS
      auto x = std::make_pair(get_efps(cycle), cycle.size());
      if (efpss.contains(x)) {
        std::cout << "Duplicate EFPS of size " << cycle.size() << " found."
                  << std::endl;
        continue;
      }
      auto ret = efpss.emplace(x);

      if (ret.second) {
        // PRint value and size of last efpss
        std::cout << "Found EFPS of size " << cycle.size() << " with value ";
        double val = 0.0;
        for (auto [u, v] : efpss.rbegin()->first)
          val += getNodeRel(y.at(std::make_pair(u, v)));
        std::cout << " = " << val << " ";
        std::cout << "with elements ";
        for (auto [u, v] : efpss.rbegin()->first)
          std::cout << "(" << u << "," << v << ") ";
        std::cout << std::endl;
      }
    }
    return efpss;
  }

  // Function that adds a set of EFPS as cuts
  std::pair<double, double>
  addCutEfpss(std::set<std::pair<EdgeList, size_t>> &efpss) {
    size_t accumCycle = 0;
    size_t accumExt = 0;

    for (auto &[efps, size] : efpss) {
      accumCycle += size;
      accumExt += efps.size();
      GRBLinExpr pathSum;
      for (auto [u, v] : efps)
        pathSum += y.at(std::make_pair(u, v));
      addCut(pathSum <= size - 1);
    }

    return std::make_pair(static_cast<double>(accumCycle) / efpss.size(),
                          static_cast<double>(accumExt) / efpss.size());
  }
};
} // end of namespace

SolveResult solveLazyEfpss(Pds &input, boost::optional<std::string> logPath,
                           std::ostream &callbackFile, std::ostream &solFile,
                           double timeLimit, bool inProp, bool outProp,
                           bool initEFPS, size_t lazyLimit, size_t cutLimit) {
  LazyEfpsCB lazyEfpss(input, callbackFile, solFile, inProp, outProp, initEFPS,
                       lazyLimit, cutLimit);
  return lazyEfpss.solve(logPath, timeLimit);
}

} // end of namespace pds
