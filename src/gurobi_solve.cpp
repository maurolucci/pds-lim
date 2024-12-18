#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {

GRBEnv &getEnv() {
  static thread_local GRBEnv env;
  return env;
}

MIPModel::MIPModel() : model(std::make_unique<GRBModel>(getEnv())), s(), w() {}
MIPModel::~MIPModel() {}

SolveResult solveMIP(const Pds &state, MIPModel &mipmodel,
                     boost::optional<std::string> outPath, double timeLimit) {
  auto &model = *mipmodel.model;
  if (outPath.has_value()) {
    model.set(GRB_StringParam_LogFile, outPath.get());
    model.set(GRB_IntParam_LogToConsole, false);
  } else {
    model.set(GRB_IntParam_LogToConsole, true);
  }
  model.set(GRB_DoubleParam_TimeLimit, timeLimit);

  model.optimize();

  switch (model.get(GRB_IntAttr_Status)) {
  case GRB_INFEASIBLE:
    return {model.get(GRB_IntAttr_NumVars),
            model.get(GRB_IntAttr_NumConstrs),
            1.0,
            0.0,
            model.get(GRB_DoubleAttr_MIPGap),
            model.get(GRB_DoubleAttr_NodeCount),
            SolveState::Infeasible};
  case GRB_OPTIMAL:
    return {model.get(GRB_IntAttr_NumVars),
            model.get(GRB_IntAttr_NumConstrs),
            model.get(GRB_DoubleAttr_ObjBound),
            model.get(GRB_DoubleAttr_ObjVal),
            model.get(GRB_DoubleAttr_MIPGap),
            model.get(GRB_DoubleAttr_NodeCount),
            SolveState::Optimal};
  case GRB_TIME_LIMIT:
    return {model.get(GRB_IntAttr_NumVars),
            model.get(GRB_IntAttr_NumConstrs),
            model.get(GRB_DoubleAttr_ObjBound),
            model.get(GRB_DoubleAttr_ObjVal),
            model.get(GRB_DoubleAttr_MIPGap),
            model.get(GRB_DoubleAttr_NodeCount),
            SolveState::Timeout};
  default:
    return {model.get(GRB_IntAttr_NumVars),
            model.get(GRB_IntAttr_NumConstrs),
            1.0,
            0.0,
            model.get(GRB_DoubleAttr_MIPGap),
            model.get(GRB_DoubleAttr_NodeCount),
            SolveState::Other};
  }
}

MIPModel brimkovModel(Pds &input) {
  MIPModel mipmodel;
  auto &model = *mipmodel.model;

  auto &s = mipmodel.s;
  auto &w = mipmodel.w;
  std::map<PowerGrid::vertex_descriptor, GRBVar> x;
  std::map<
      std::pair<PowerGrid::vertex_descriptor, PowerGrid::vertex_descriptor>,
      GRBVar>
      y;

  auto &graph = static_cast<const Pds &>(input).get_graph();
  auto n_channels = input.get_n_channels();
  size_t T = boost::num_vertices(graph);

  // Add variables
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    s.try_emplace(
        v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
    x.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER,
                                  fmt::format("x_{}", v)));
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
    model.addConstr(constr1 == 1);

    // (2) x_w - x_v + (T+1)y_u_v <= T, \forall (u,v) \in A_Z, w \in N[u] -
    // v
    if (input.isZeroInjection(v)) {
      for (auto u :
           boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
        GRBLinExpr constr2 = 0;
        constr2 += x.at(v) - x.at(u) + (T + 1) * y.at(std::make_pair(v, u));
        model.addConstr(constr2 <= T);
        for (auto w :
             boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
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
    for (auto u :
         boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
      constr3 += w.at(std::make_pair(v, u));
    }
    model.addConstr(constr3 <=
                    std::min(boost::degree(v, graph), (n_channels - 1)) *
                        s.at(v));
  }

  return mipmodel;
}

} // end of namespace pds