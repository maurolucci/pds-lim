#include "graphio.hpp"
#include "gurobi_solve.hpp"
#include "pds.hpp"

#include <boost/program_options.hpp>
#include <chrono>
#include <filesystem>
#include <fmt/format.h>
#include <iostream>
#include <limits>

namespace fs = std::filesystem;
namespace po = boost::program_options;
using namespace pds;
using Solver = std::function<SolveResult(Pds &, double)>;

std::string format_solve_state(SolveState state) {
  std::string name = "unknown";
  switch (state) {
  case pds::SolveState::Optimal:
    name = "Optimal";
    break;
  case pds::SolveState::Other:
    name = "Other";
    break;
  case pds::SolveState::Infeasible:
    name = "Infeasible";
    break;
  case pds::SolveState::Timeout:
    name = "Timeout";
    break;
  }
  return name;
}

auto now() { return std::chrono::high_resolution_clock::now(); }

template <typename T> auto µs(T time) {
  return std::chrono::duration_cast<std::chrono::microseconds>(time).count();
}

auto getModel(const std::string &name) {
  if (name == "brimkov") {
    return brimkovModel;
  } else {
    throw std::invalid_argument("unknown model " + name);
  }
}

auto getSolver(po::variables_map &vm) {
  std::string solverName = vm["solver"].as<std::string>();
  try {
    return Solver{[model = getModel(solverName)](auto &input, double timeout) {
      auto mip = model(input);
      auto result = solveMIP(input, mip, true, timeout);
      return result;
    }};
  } catch (std::invalid_argument &ex) {
    fmt::print(stderr, "{}", ex.what());
    throw ex;
  }
}

int main(int argc, const char **argv) {

  // Parse arguments
  po::options_description desc(argv[0]);
  desc.add_options()("help,h", "show this help");
  desc.add_options()("solver,s", po::value<std::string>()->required(),
                     "solver, can be any of [brimkov,]");
  desc.add_options()("n-channels,w", po::value<size_t>()->required(),
                     "number of channels");
  desc.add_options()(
      "graph,f",
      po::value<std::vector<std::string>>()->required()->multitoken(),
      "input files");
  desc.add_options()("all-zi,z", "consider all nodes zero-innjection");
  desc.add_options()("repeat,n",
                     po::value<size_t>()->default_value(1)->implicit_value(5),
                     "number of experiment repetitions");
  desc.add_options()("timeout,t", po::value<double>()->default_value(3600.0),
                     "gurobi time limit (seconds)");
  po::positional_options_description pos;
  pos.add("graph", -1);
  po::variables_map vm;

  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    po::notify(vm);
  } catch (po::error const &e) {
    std::cerr << "Invalid arguments. Please, write './pds-lim -h' for help.\n";
    std::cerr << e.what() << std::endl;
    return 2;
  }

  // Read arguments
  if (vm.count("help")) {
    desc.print(std::cout);
    return 1;
  }
  std::string solver = vm["solver"].as<std::string>();
  bool allZeroInjection = vm.count("all-zi");
  size_t repetitions = vm["repeat"].as<size_t>();
  double timeout = vm["timeout"].as<double>();
  size_t n_channels = vm["n-channels"].as<size_t>();
  std::vector<std::string> inputs;
  if (vm.count("graph")) {
    inputs = vm["graph"].as<std::vector<std::string>>();
  } else {
    fmt::print(stderr, "no input\n");
    return 2;
  }

  // Read inputs
  for (const std::string &filename : inputs) {

    std::string currentName = fs::path(filename).filename().string();
    currentName = currentName.substr(0, currentName.rfind('.'));
    Pds input(readGraphML(filename, allZeroInjection), n_channels);
    for (size_t run = 0; run < repetitions; ++run) {

      auto &graph = static_cast<const Pds &>(input).get_graph();
      size_t n = boost::num_vertices(graph);
      size_t m = boost::num_edges(graph);
      size_t zi = input.numZeroInjection();

      auto t0 = now();
      Solver solve = getSolver(vm);
      auto result = solve(input, timeout);
      auto t1 = now();

      using namespace fmt::literals;
      fmt::print("{solver},{name},{n},{m},{zi},{channels},{variables},{"
                 "constraints},{run},{lower_bound},{upper_bound},{gap},{result}"
                 ",{nodes},{t_solver}\n",
                 "solver"_a = solver, "name"_a = filename, "n"_a = n, "m"_a = m,
                 "zi"_a = zi, "channels"_a = n_channels,
                 "variables"_a = result.variables,
                 "constraints"_a = result.constraints, "run"_a = run,
                 "lower_bound"_a = result.lower, "upper_bound"_a = result.upper,
                 "gap"_a = result.gap,
                 "result"_a = format_solve_state(result.state),
                 "nodes"_a = result.nodes, "t_solver"_a = µs(t1 - t0));
    }
  }

  return 0;
}