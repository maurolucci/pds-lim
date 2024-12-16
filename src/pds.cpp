#include "pds.hpp"

namespace pds {

Pds::Pds() : Pds(PowerGrid{}, 0) {}

Pds::Pds(const pds::PowerGrid &graph, size_t n_channels)
    : Pds(PowerGrid{graph}, n_channels) {}

Pds::Pds(PowerGrid &&graph, size_t n_channels)
    : graph(graph), n_channels(n_channels) {}

} // end of namespace pds