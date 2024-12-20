#ifndef GUROBI_COMMON_HPP
#define GUROBI_COMMON_HPP

#include "pds.hpp"
#include <map>

struct GRBModel;
struct GRBVar;
struct GRBEnv;

namespace pds {

struct MIPModel {
  std::unique_ptr<GRBModel> model;
  std::map<Vertex, GRBVar> s;
  std::map<Edge, GRBVar> w;
  MIPModel();
  MIPModel(MIPModel &&other) = default;
  virtual ~MIPModel();
};

void preloadMIPSolver();
GRBEnv &getEnv();

void relaxMIPModel(MIPModel &);

SolveResult solveMIP(Pds &, MIPModel &, boost::optional<std::string>, double);

} // namespace pds

#endif // GUROBI_COMMON_HPP
