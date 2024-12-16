#include "graphio.hpp"
#include "gurobi_solve.hpp"
#include "pds.hpp"

#include <iostream>

int main() {

  using namespace pds;
  std::string instance_path = "inputs/case5.graphml";
  Pds pds(readGraphML(instance_path, false), 2);

  MIPModel model = brimkovModel(pds);
  SolveResult ret = solveMIP(pds, model);

  std::string state;
  switch (ret.state) {
  case SolveState::Optimal:
    state = "Optimal";
    break;
  case SolveState::Timeout:
    state = "Timeout";
    break;
  case SolveState::Infeasible:
    state = "Infeasible";
    break;
  default:
    state = "Other";
    break;
  }
  std::cout << "Lower: " << ret.lower << ", Upper: " << ret.upper << ", "
            << state << std::endl;

  return 0;
}