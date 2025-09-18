#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {

GRBEnv &getEnv() {
  static thread_local GRBEnv env;
  return env;
}

MIPModel::MIPModel()
    : model(std::make_unique<GRBModel>(getEnv())), s(), w(),
      totalLazyCBCalls(0), totalLazyCBTime(0), totalLazyAdded(0),
      totalCutCBCalls(0), totalCutCBTime(0), totalCutAdded(0) {}
MIPModel::~MIPModel() {}

SolveResult solveMIP(Pds &input, MIPModel &mipmodel,
                     boost::optional<std::string> outPath,
                     std::ostream &solFile, double timeLimit) {
  auto &model = *mipmodel.model;
  if (outPath.has_value()) {
    model.set(GRB_IntParam_LogToConsole, false);
    model.set(GRB_StringParam_LogFile, outPath.get());
  } else {
    model.set(GRB_IntParam_LogToConsole, true);
  }
  model.set(GRB_DoubleParam_TimeLimit, timeLimit);

  model.optimize();
  SolveResult result = {model.get(GRB_IntAttr_NumVars),
                        model.get(GRB_IntAttr_NumConstrs),
                        -1.0,
                        -1.0,
                        100.0,
                        model.get(GRB_DoubleAttr_NodeCount),
                        SolveState::Other,
                        mipmodel.totalLazyCBCalls,
                        mipmodel.totalLazyCBTime,
                        mipmodel.totalLazyAdded,
                        mipmodel.totalCutCBCalls,
                        mipmodel.totalCutCBTime,
                        mipmodel.totalCutAdded};

  switch (model.get(GRB_IntAttr_Status)) {
  case GRB_INFEASIBLE:
    result.state = SolveState::Infeasible;
    break;
  case GRB_OPTIMAL:
    result.lower =
        static_cast<size_t>(round(model.get(GRB_DoubleAttr_ObjBound)));
    result.upper = static_cast<size_t>(round(model.get(GRB_DoubleAttr_ObjVal)));
    result.gap = model.get(GRB_DoubleAttr_MIPGap);
    result.state = SolveState::Optimal;
    break;
  case GRB_TIME_LIMIT:
    result.lower = model.get(GRB_DoubleAttr_ObjBound);
    result.state = SolveState::Timeout;
    if (model.get(GRB_IntAttr_SolCount) > 0) {
      result.upper =
          static_cast<size_t>(round(model.get(GRB_DoubleAttr_ObjVal)));
      result.gap = model.get(GRB_DoubleAttr_MIPGap);
    }
    break;
  }

  // Check correctness of solution
  if (model.get(GRB_IntAttr_SolCount) > 0) {

    // Recover variable values
    const PowerGrid &graph = input.get_graph();
    std::map<Vertex, double> sValue;
    std::map<Edge, double> wValue;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      sValue[v] = mipmodel.s.at(v).get(GRB_DoubleAttr_X);
      if (degree(v, graph) <= input.get_n_channels())
        continue;
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
        wValue[std::make_pair(v, u)] =
            mipmodel.w.at(std::make_pair(v, u)).get(GRB_DoubleAttr_X);
      }
    }

    VertexList mS = input.get_monitored_set(sValue, wValue);
    if (!input.isFeasible(mS))
      throw std::runtime_error("Error: The solution IS NOT FEASIBLE");

    // Print solution
    mipmodel.write_sol(solFile);
  }

  return result;
}

MIPModel brimkovModel(Pds &input) {
  MIPModel mipmodel;
  auto &model = *mipmodel.model;

  auto &s = mipmodel.s;
  auto &w = mipmodel.w;
  std::map<Vertex, GRBVar> x;
  std::map<Edge, GRBVar> y;

  auto &graph = static_cast<const Pds &>(input).get_graph();
  auto n_channels = input.get_n_channels();
  size_t T = boost::num_vertices(graph);

  // Add variables
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    s.try_emplace(
        v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
    x.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER,
                                  fmt::format("x_{}", v)));
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

    // (1) s_v + sum_{u \in N(v) \cap V_1} s_u + sum_{u \in N(v) \cap V_2} w_u_v
    // + sum_{u \in N(v) \cap V_Z} y_u_v == 1, \forall v \in V
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

    // (2) x_w - x_u + (T+1)y_v_u <= T, \forall (v,u) \in A_Z, w \in N[v] -
    // u
    if (input.isZeroInjection(v)) {
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
        GRBLinExpr constr2 = 0;
        constr2 += x.at(v) - x.at(u) + (T + 1) * y.at(std::make_pair(v, u));
        model.addConstr(constr2 <= T);
        for (auto w : boost::make_iterator_range(adjacent_vertices(v, graph))) {
          if (w == u)
            continue;
          GRBLinExpr constr2n = 0;
          constr2n += x.at(w) - x.at(u) + (T + 1) * y.at(std::make_pair(v, u));
          model.addConstr(constr2n <= T);
        }
      }
    }

    // (3) sum_{u \in N(v)} w_v_u <= omega_v * s_v, \forall v \in V2
    if (degree(v, graph) > n_channels) {
      GRBLinExpr constr3 = 0;
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        constr3 += w.at(std::make_pair(v, u));
      model.addConstr(constr3 == n_channels * s.at(v));
    }
  }

  return mipmodel;
}

MIPModel jovanovicModel(Pds &input) {
  MIPModel mipmodel;
  auto &model = *mipmodel.model;

  auto &s = mipmodel.s;
  auto &w = mipmodel.w;
  std::map<Vertex, GRBVar> x;
  std::map<Edge, GRBVar> y;

  auto &graph = static_cast<const Pds &>(input).get_graph();
  auto n_channels = input.get_n_channels();
  size_t M = boost::num_vertices(graph);

  // Add variables
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    s.try_emplace(
        v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
    x.try_emplace(v, model.addVar(1.0, static_cast<double>(M), 0.0, GRB_INTEGER,
                                  fmt::format("x_{}", v)));
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

    // (1) x_v >= 1, \forall v \in V
    // Already implied by the lower bound of the variable

    // (2) x_v <= s_v + M(1-s_v), \forall v in V
    GRBLinExpr constr2 = 0;
    constr2 += x.at(v) - s.at(v) - M * (1 - s.at(v));
    model.addConstr(constr2 <= 0);

    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
      // (3.1) x_v <= s_u + M(1-s_u), \forall v in V, u \in N(v) \cap V1
      if (degree(u, graph) <= n_channels) {
        GRBLinExpr constr3 = 0;
        constr3 += x.at(v) - s.at(u) - M * (1 - s.at(u));
        model.addConstr(constr3 <= 0);
      }
      // (3.2) x_v <= w_u_v + M(1-w_u_v), \forall v in V, u \in N(v) \cap V2
      else {
        GRBLinExpr constr3 = 0;
        constr3 += x.at(v) - w.at(std::make_pair(u, v)) -
                   M * (1 - w.at(std::make_pair(u, v)));
        model.addConstr(constr3 <= 0);
      }
    }

    // (4) x_v <= M(s_v + sum_{u \in N(v) \cap V1} s_u + sum_{u \in N(v) \cap
    // V2} w_u_v + sum_{u \in N(v) \cap V_Z} y_u_v), \forall v \in V
    GRBLinExpr constr4 = 0;
    constr4 += x.at(v) - M * s.at(v);
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
      if (degree(u, graph) <= n_channels)
        constr4 -= M * s.at(u);
      else
        constr4 -= M * w.at(std::make_pair(u, v));
      if (input.isZeroInjection(u))
        constr4 -= M * y.at(std::make_pair(u, v));
    }
    model.addConstr(constr4 <= 0);

    // (5) sum_{u \in N(v)} y_v_u <= 1, \forall v \in V_Z
    if (input.isZeroInjection(v)) {
      GRBLinExpr constr5 = 0;
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        constr5 += y.at(std::make_pair(v, u));
      model.addConstr(constr5 <= 1);
    }

    // (6) sum_{u \in N(v) \cap V_Z} y_u_v <= 1, \forall v \in V
    GRBLinExpr constr6 = 0;
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      if (input.isZeroInjection(u))
        constr6 += y.at(std::make_pair(u, v));
    model.addConstr(constr6 <= 1);

    // (7) y_u_v + y_v_u <= 1, \forall (u,v),(v,u) in A_Z
    if (input.isZeroInjection(v))
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        if (input.isZeroInjection(u)) {
          GRBLinExpr constr7 = 0;
          constr7 += y.at(std::make_pair(u, v)) + y.at(std::make_pair(v, u));
          model.addConstr(constr7 <= 1);
        }

    // (8) x_u >= x_w + 1 - M(1-y_v_u), \forall (v,u) \in A_Z, w \in N[v]-u
    if (input.isZeroInjection(v)) {
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
        GRBLinExpr constr8 = 0;
        constr8 += x.at(u) - x.at(v) - 1 + M * (1 - y.at(std::make_pair(v, u)));
        model.addConstr(constr8 >= 0);
        for (auto w : boost::make_iterator_range(adjacent_vertices(v, graph))) {
          if (w == u)
            continue;
          GRBLinExpr constr8n = 0;
          constr8n +=
              x.at(u) - x.at(w) - 1 + M * (1 - y.at(std::make_pair(v, u)));
          model.addConstr(constr8n >= 0);
        }
      }
    }

    // (9) sum_{u \in N(v)} w_v_u <= omega_v * s_v, \forall v \in V2
    if (degree(v, graph) > n_channels) {
      GRBLinExpr constr9 = 0;
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        constr9 += w.at(std::make_pair(v, u));
      model.addConstr(constr9 == n_channels * s.at(v));
    }
  }

  return mipmodel;
}

} // end of namespace pds
