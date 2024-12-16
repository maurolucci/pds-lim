#include "graphio.hpp"
#include "gurobi_solve.hpp"
#include "pds.hpp"

#include <iostream>

int main () {

    using namespace pds;
    std::string instance_path = "inputs/case5.graphml";
    Pds pds(readGraphML(instance_path, false));

    auto model = brimkovModel(pds);
    solveMIP(pds, model);

    return 0;
}