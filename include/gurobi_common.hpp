#ifndef GUROBI_COMMON_HPP
#define GUROBI_COMMON_HPP

#include <chrono>
#include <iostream>
#include <map>

#include "gurobi_c++.h"
#include "pds.hpp"

namespace pds {

struct MIPModel {
  std::unique_ptr<GRBModel> model;
  std::map<Vertex, GRBVar> s;
  std::map<Edge, GRBVar> w;
  MIPModel();
  MIPModel(MIPModel &&other) = default;
  virtual ~MIPModel();

  void write_sol(std::ostream &solFile) {
    for (auto &[v, var] : s)
      if (var.get(GRB_DoubleAttr_X) > 0.5) {
        solFile << v << ": ";
        for (auto &[e, var2] : w)
          if (e.first == v && var2.get(GRB_DoubleAttr_X) > 0.5)
            solFile << e.second << " ";
        solFile << std::endl;
      }
  }
};

void preloadMIPSolver();
GRBEnv &getEnv();

void relaxMIPModel(MIPModel &);

SolveResult solveMIP(Pds &, MIPModel &, boost::optional<std::string>,
                     std::ostream &solFile, double);

} // namespace pds

#endif // GUROBI_COMMON_HPP
