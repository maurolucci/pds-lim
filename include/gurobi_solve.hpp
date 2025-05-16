#ifndef GUROBI_SOLVE_HPP
#define GUROBI_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

#include <boost/optional.hpp>

namespace pds {

MIPModel brimkovModel(Pds &input);
MIPModel brimkovModel2(Pds &input);

MIPModel jovanovicModel(Pds &inputs);
MIPModel jovanovicModel2(Pds &inputs);

SolveResult solveMIP(const Pds &input, MIPModel &model,
                     boost::optional<std::string> output, double timeLimit);

} // namespace pds

#endif // GUROBI_SOLVE_HPP
