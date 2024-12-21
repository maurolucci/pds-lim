#include "cycle_solve.hpp"
#include "graphio.hpp"
#include "gurobi_solve.hpp"
#include "pds.hpp"

#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <gurobi_c++.h>
#include <iostream>
#include <limits>
#include <map>

namespace fs = std::filesystem;
namespace po = boost::program_options;
using namespace pds;
using Solver =
    std::function<SolveResult(Pds &, boost::optional<std::string>, double)>;

std::map<std::string, fs::path> outDirs = {{"log", fs::path()},
                                           {"stat", fs::path()},
                                           {"sol", fs::path()},
                                           {"cb", fs::path()}};

class outPut {
public:
  std::ofstream cbFileAux;
  std::ostream &cbFile;
  outPut(std::ostream &cbStream = std::cout) : cbFileAux(), cbFile(cbStream) {}
  outPut(std::string sbPath)
      : cbFileAux(sbPath, std::ofstream::app), cbFile(this->cbFileAux) {}
};

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
    return Solver{
        [model = getModel(solverName)](
            auto &input, boost::optional<std::string> logPath, double timeout) {
          auto mip = model(input);
          auto result = solveMIP(input, mip, logPath, timeout);
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
  desc.add_options()("repeat,n", po::value<size_t>()->default_value(1),
                     "number of experiment repetitions");
  desc.add_options()("timeout,t", po::value<double>()->default_value(3600.0),
                     "gurobi time limit (seconds)");
  desc.add_options()("outdir,o", po::value<std::string>(),
                     "write outputs to the specified directory");
  po::positional_options_description pos;
  pos.add("graph", -1);
  po::variables_map vm;

  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    po::notify(vm);
  } catch (po::error const &e) {
    std::cerr << e.what() << std::endl;
    desc.print(std::cout);
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
  fs::path outDir;
  if (vm.count("outdir")) {
    outDir = fs::path(vm["outdir"].as<std::string>());
    // Create output directories
    if (!fs::is_directory(outDir)) {
      fs::create_directory(outDir);
    }
    for (auto [key, _] : outDirs) {
      outDirs[key] = outDir / key / "";
      if (!fs::is_directory(outDirs[key])) {
        fs::create_directories(outDirs[key]);
      }
    }
  }

  // Read inputs
  for (const std::string &filename : inputs) {
    Pds input(readGraphML(filename, allZeroInjection), n_channels);
    for (size_t run = 0; run < repetitions; ++run) {

      fs::path currentName(fs::path(filename).stem().string() +
                           fmt::format("-{}-{}-{}", solver, n_channels, run));
      std::cout << "Solving Instance " << currentName << " ..." << std::endl;
      if (vm.count("outdir")) {
        for (auto [key, _] : outDirs) {
          outDirs[key].replace_filename(currentName);
          outDirs[key].replace_extension(key);
        }

        // Resume check
        if (fs::is_regular_file(outDirs["stat"])) {
          std::ifstream stats;
          stats.open(outDirs["stat"]);
          std::string line;
          while (getline(stats, line)) {
            std::cout << line << std::endl;
          }
          continue;
        }
      }

      auto &graph = static_cast<const Pds &>(input).get_graph();
      size_t n = boost::num_vertices(graph);
      size_t m = boost::num_edges(graph);
      size_t zi = input.numZeroInjection();

      // Prepare output files
      boost::optional<std::string> logPath =
          vm.count("outdir")
              ? boost::optional<std::string>(outDirs["log"].string())
              : boost::none;
      outPut output =
          (vm.count("outdir")) ? outPut(outDirs["cb"].string()) : outPut();

      // Solve
      auto t0 = now();
      SolveResult result;
      std::string solverName = vm["solver"].as<std::string>();
      try {
        if (solverName == "cycles") {
          result = solveLazyCycles(input, logPath, output.cbFile, timeout);
        } else {
          Solver solve = getSolver(vm);
          result = solve(input, logPath, timeout);
        }
      } catch (GRBException ex) {
        fmt::print(stderr, "Gurobi Exception {}: {}\n", ex.getErrorCode(),
                   ex.getMessage());
        throw ex;
      }
      auto t1 = now();

      // Write stats
      using namespace fmt::literals;
      std::string stats(fmt::format(
          "{solver},{name},{n},{m},{zi},{channels},{variables},{"
          "constraints},{run},{lower_bound},{upper_bound},{gap},{result}"
          ",{nodes},{t_solver}\n",
          "solver"_a = solver, "name"_a = filename, "n"_a = n, "m"_a = m,
          "zi"_a = zi, "channels"_a = n_channels,
          "variables"_a = result.variables,
          "constraints"_a = result.constraints, "run"_a = run,
          "lower_bound"_a = result.lower, "upper_bound"_a = result.upper,
          "gap"_a = result.gap, "result"_a = format_solve_state(result.state),
          "nodes"_a = result.nodes, "t_solver"_a = µs(t1 - t0)));
      if (vm.count("outdir")) {
        std::ofstream statFile(outDirs["stat"], std::ofstream::out);
        statFile << stats;
        statFile.close();
      }
      std::cout << stats;
    }
  }

  return 0;
}