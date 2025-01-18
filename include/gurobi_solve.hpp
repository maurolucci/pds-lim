#ifndef GUROBI_SOLVE_HPP
#define GUROBI_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

#include <boost/optional.hpp>
#include <chrono>

namespace pds {

auto now();
template <typename T> auto µs(T time);

MIPModel brimkovModel(Pds &input);

SolveResult solveMIP(const Pds &input, MIPModel &model,
                     boost::optional<std::string> output, double timeLimit);

} // namespace pds

#endif // GUROBI_SOLVE_HPP