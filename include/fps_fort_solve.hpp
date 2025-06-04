#ifndef FPS_FORT_SOLVE_HPP
#define FPS_FORT_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {
SolveResult solveLazyFpssForts(Pds &, boost::optional<std::string>,
                               std::ostream &, std::ostream &, double, size_t);
} // end of namespace pds

#endif // FPS2_SOLVE_HPP
