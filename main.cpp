#include "pds.hpp"
#include "graphio.hpp"
#include <iostream>

int main () {

    using namespace pds;
    std::string instance_path = "inputs/case1888rte.graphml";
    Pds pds(readGraphML(instance_path, false));

    std::cout << boost::num_vertices(pds.get_graph()) << std::endl;
}