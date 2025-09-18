#ifndef EFPS_SOLVE_HPP
#define EFPS_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {
SolveResult solveLazyEfpss(Pds &, boost::optional<std::string>, std::ostream &,
                           std::ostream &, double, bool, bool, bool, size_t,
                           bool, size_t);
} // end of namespace pds

#endif // EFPS_SOLVE_HPP
