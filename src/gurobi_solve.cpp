#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {

GRBEnv &getEnv() {
  static thread_local GRBEnv env;
  return env;
}

MIPModel::MIPModel() : model(std::make_unique<GRBModel>(getEnv())), s(), w(),
  totalCallback(0), totalCallbackTime(0), totalLazy(0) {}
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
                        mipmodel.totalCallback, 
                        mipmodel.totalCallbackTime,
                        mipmodel.totalLazy};

  switch (model.get(GRB_IntAttr_Status)) {
  case GRB_INFEASIBLE:
    result.state = SolveState::Infeasible;
    break;
  case GRB_OPTIMAL:
    result.lower = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound));
    result.upper = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal));
    result.gap = model.get(GRB_DoubleAttr_MIPGap);
    result.state = SolveState::Optimal;
    break;
  case GRB_TIME_LIMIT:
    result.lower = model.get(GRB_DoubleAttr_ObjBound);
    result.state = SolveState::Timeout;
    if (model.get(GRB_IntAttr_SolCount) > 0) {
      result.upper = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal));
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
      for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
        wValue[std::make_pair(v, u)] =
            mipmodel.w.at(std::make_pair(v, u)).get(GRB_DoubleAttr_X);
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
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph))) {
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

    // (2) x_w - x_v + (T+1)y_u_v <= T, \forall (u,v) \in A_Z, w \in N[u] -
    // v
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

    // (3) sum_{u \in N(v)} w_v_u <= (omega_v - 1) s_v, \forall v \in V
    GRBLinExpr constr3 = 0;
    for (auto u : boost::make_iterator_range(adjacent_vertices(v, graph)))
      constr3 += w.at(std::make_pair(v, u));
    model.addConstr(constr3 ==
                    std::min(degree(v, graph), (n_channels - 1)) * s.at(v));
  }

  return mipmodel;
}

} // end of namespace pds
