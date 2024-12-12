#include "pds.hpp"

namespace pds {

Pds::Pds() : Pds(PowerGrid{}) { }

Pds::Pds(const pds::PowerGrid &graph) : Pds(PowerGrid{graph}) { }

Pds::Pds(PowerGrid&& graph) : graph(graph) { }

} // end of namespace pds