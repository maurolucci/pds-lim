#ifndef CYCLE_SOLVE_HPP
#define CYCLE_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {
SolveResult solveLazyCycles(Pds &, boost::optional<std::string>, std::ostream &,
                            double, size_t);
} // end of namespace pds

#endif // CYCLE_SOLVE_HPP
