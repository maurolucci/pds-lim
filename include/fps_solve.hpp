#ifndef FPS_SOLVE_HPP
#define FPS_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {
SolveResult solveLazyFpss(Pds &, boost::optional<std::string>, std::ostream &,
                          std::ostream &, double, bool, bool, bool, bool, bool,
                          size_t);
} // end of namespace pds

#endif // FPS_SOLVE_HPP
