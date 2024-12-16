#ifndef GUROBI_SOLVE_HPP
#define GUROBI_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {

inline const double TIME_LIMIT = 3600;

MIPModel brimkovModel(Pds &input);

SolveResult solveMIP(const Pds &input, MIPModel &model, bool output = false,
                     double timeLimit = TIME_LIMIT);

} // namespace pds

#endif // GUROBI_SOLVE_HPP