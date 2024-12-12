#ifndef GRAPHIO_HPP
#define GRAPHIO_HPP

#include "pds.hpp"

namespace pds {

PowerGrid readGraphML(const std::string& filename, bool allZeroInjection = false);

}

#endif // GRAPHIO_HPP