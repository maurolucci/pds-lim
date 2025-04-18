#ifndef FORT_SOLVE_HPP
#define FORT_SOLVE_HPP

#include "gurobi_common.hpp"
#include "pds.hpp"

namespace pds {
SolveResult solveLazyForts(Pds &, boost::optional<std::string>, std::ostream &,
                           std::ostream &, double, size_t);
} // end of namespace pds

#endif // FORT_SOLVE_HPP
