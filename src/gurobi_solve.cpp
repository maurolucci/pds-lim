#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {

GRBEnv &getEnv() {
    static thread_local GRBEnv env;
    return env;
}

MIPModel::MIPModel() : model(std::make_unique<GRBModel>(getEnv())), s(), w() { }
MIPModel::~MIPModel() { }

SolveResult solveMIP(const Pds& state, MIPModel & mipmodel, bool output, double timeLimit) {
    auto& model = *mipmodel.model;
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);

    model.optimize();

    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_INFEASIBLE:
            return { size_t{1}, size_t{0}, SolveState::Infeasible };
        case GRB_OPTIMAL:
            return { size_t(model.get(GRB_DoubleAttr_ObjBound)), size_t(model.get(GRB_DoubleAttr_ObjVal)), SolveState::Optimal };
        case GRB_TIME_LIMIT:
            return { size_t(model.get(GRB_DoubleAttr_ObjBound)), size_t(model.get(GRB_DoubleAttr_ObjVal)), SolveState::Timeout };
        default:
            return { size_t{1}, size_t{0}, SolveState::Other };
    }
}

MIPModel brimkovModel(Pds& input) {
    try {

        MIPModel mipmodel;
        auto& model = *mipmodel.model;

        auto& s = mipmodel.s;
        auto& w = mipmodel.w;
        std::map<PowerGrid::vertex_descriptor, GRBVar> x;
        std::map<std::pair<PowerGrid::vertex_descriptor, PowerGrid::vertex_descriptor>, GRBVar> y;

        auto &graph = static_cast<const Pds&>(input).get_graph();
        size_t T = boost::num_vertices(graph);
        
        // Add variables
        for (auto v: boost::make_iterator_range(vertices(graph))) {
            std::cout << v << std::endl;
        }

/*
        for (auto v: state.graph().vertices()) {
            xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, "x"));
            si.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER, "s"));
            for (auto e: state.graph().outEdges(v)) {
                assert(!ye.contains(e));
                ye.try_emplace(e, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
                if (!state.isZeroInjection(v)) {
                    model.addConstr(ye.at(e) <= xi.at(v));
                }
            }
            if (state.isActive(v)) {
                model.addConstr(xi.at(v) == 1.0);
            }
            else if (state.isInactive(v)) {
                model.addConstr(xi.at(v) == 0.0);
            }
        }
        for (auto u: state.graph().vertices()) {
            // x_u + \sum_{v \in N(u)} y_{vu} = 1 (3)
            GRBLinExpr observers = 0;
            observers += xi.at(u);
            for (auto e: state.graph().inEdges(u)) {
                observers += ye.at(e);
            }
            model.addConstr(observers == 1);
            for (auto e: state.graph().outEdges(u)) {
                auto v = state.graph().target(e);
                // s_u - s_v + (T + 1) y_{uv} <= T (4) e@(u,v) \in E'
                model.addConstr(si.at(u) - si.at(v) + (T + 1) * ye.at(e) <= T);
                for (auto w: state.graph().neighbors(u)) {
                    if (v != w) {
                        unused(u, w, v, e);
                        // s_w - s_v + (T + 1) y_e <= T + (T+1) x_u , e@(u,v) \in E', w \in N(u) - v
                        model.addConstr(si.at(w) - si.at(v) + (T + 1) * ye.at(e) <= T + (T+1) * xi.at(u));
                    }
                }
            }
        }
*/


        return mipmodel;
    } catch (GRBException ex) {
        fmt::print(stderr, "Gurobi Exception {}: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }
}

} // end of namespace pds