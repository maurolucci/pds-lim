#ifndef PDS_HPP
#define PDS_HPP

#include <fmt/core.h>
#include <fmt/ranges.h>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"

namespace pds {

struct Bus {
    std::string name;
    long id;
    bool zero_injection;
};

using PowerGrid = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Bus>;

enum class SolveState {
    Optimal,
    Timeout,
    Infeasible,
    Other
};

struct SolveResult {
    size_t lower;
    size_t upper;
    SolveState state;
};

class Pds {
private:
    PowerGrid graph;
public:
    Pds();
    explicit Pds(PowerGrid&& graph);
    explicit Pds(const PowerGrid& graph);

    [[nodiscard]] inline const PowerGrid& get_graph() const { return graph; }
}; // end of class Pds

} // end of namespace pds

#endif // PDS_HPP