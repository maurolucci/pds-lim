#include "fort_solve.hpp"
#include "fps_fort_solve.hpp"
#include "fps_solve.hpp"
#include "graphio.hpp"
#include "gurobi_common.hpp"
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
using Solver = std::function<SolveResult(Pds &, boost::optional<std::string>,
                                         std::ostream &, double)>;

std::map<std::string, fs::path> outDirs = {{"log", fs::path()},
                                           {"stat", fs::path()},
                                           {"sol", fs::path()},
                                           {"cb", fs::path()}};

class outPut {
public:
  std::ofstream cbFileAux, solFileAux;
  std::ostream &cbFile, &solFile;
  outPut(std::ostream &cbStream = std::cout,
         std::ostream &solStream = std::cout)
      : cbFileAux(), solFileAux(), cbFile(cbStream), solFile(solStream) {}
  outPut(std::string sbPath, std::string solPath)
      : cbFileAux(sbPath, std::ofstream::app),
        solFileAux(solPath, std::ofstream::app), cbFile(this->cbFileAux),
        solFile(this->solFileAux) {}
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

auto getModel(const std::string &name) {
  if (name == "brimkov") {
    return brimkovModel;
  } else if (name == "jovanovic") {
    return jovanovicModel;
  } else {
    throw std::invalid_argument("unknown model " + name);
  }
}

auto getSolver(po::variables_map &vm) {
  std::string solverName = vm["solver"].as<std::string>();
  try {
    return Solver{[model = getModel(solverName)](
                      auto &input, boost::optional<std::string> logPath,
                      std::ostream &solFile, double timeout) {
      auto mip = model(input);
      auto result = solveMIP(input, mip, logPath, solFile, timeout);
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
  desc.add_options()(
      "solver,s", po::value<std::string>()->required(),
      "solver, can be any of "
      "[brimkov,brimkov2,jovanovic,jovanovic2,fpss,fpss2,cycles,forts,forts2]");
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
  desc.add_options()("in-prop", "consider incoming propagation constraints");
  desc.add_options()("out-prop", "consider outgoing propagation constraints");
  desc.add_options()("init-fps-1", "consider initial FPS constraints type 1");
  desc.add_options()("init-fps-2", "consider initial FPS constraints type 2");
  desc.add_options()("init-fps-3", "consider initial FPS constraints type 3");
  desc.add_options()(
      "lazy-limit",
      po::value<size_t>()->default_value(std::numeric_limits<size_t>::max()),
      "maximum number of lazy contraints added per callback");
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
  bool inProp = vm.count("in-prop");
  bool outProp = vm.count("out-prop");
  bool initFPS1 = vm.count("init-fps-1");
  bool initFPS2 = vm.count("init-fps-2");
  bool initFPS3 = vm.count("init-fps-3");
  size_t lazyLimit = vm["lazy-limit"].as<size_t>();
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

  // Update solver name, if necessary
  if (inProp)
    solver.append("-inp");
  if (outProp)
    solver.append("-outp");
  if (initFPS1)
    solver.append("-fps1");
  if (initFPS2)
    solver.append("-fps2");
  if (initFPS3)
    solver.append("-fps3");

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
      outPut output = (vm.count("outdir")) ? outPut(outDirs["cb"].string(),
                                                    outDirs["sol"].string())
                                           : outPut();

      // Solve
      auto t0 = std::chrono::high_resolution_clock::now();
      ;
      SolveResult result;
      std::string solverName = vm["solver"].as<std::string>();
      try {
        if (solverName == "efpss") {
          result = solveLazyEfpss(input, logPath, output.cbFile, output.solFile,
                                  timeout, inProp, outProp, lazyLimit);
        } else if (solverName == "fpss") {
          result = solveLazyFpss(input, logPath, output.cbFile, output.solFile,
                                 timeout, inProp, outProp, initFPS1, initFPS2,
                                 initFPS3, lazyLimit);
        } else if (solverName == "forts") {
          result = solveLazyForts(input, logPath, output.cbFile, output.solFile,
                                  timeout, lazyLimit);
        } else {
          Solver solve = getSolver(vm);
          result = solve(input, logPath, output.solFile, timeout);
        }
      } catch (GRBException ex) {
        fmt::print(stderr, "Gurobi Exception {}: {}\n", ex.getErrorCode(),
                   ex.getMessage());
        throw ex;
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      ;

      // Write stats
      using namespace fmt::literals;
      std::string stats(fmt::format(
          "{solver},{name},{n},{m},{zi},{channels},{variables},{"
          "constraints},{run},{lower_bound},{upper_bound},{gap},{result}"
          ",{nodes},{t_solver},{callback},{t_callback},{lazy}\n",
          "solver"_a = solver, "name"_a = filename, "n"_a = n, "m"_a = m,
          "zi"_a = zi, "channels"_a = n_channels,
          "variables"_a = result.variables,
          "constraints"_a = result.constraints, "run"_a = run,
          "lower_bound"_a = result.lower, "upper_bound"_a = result.upper,
          "gap"_a = result.gap, "result"_a = format_solve_state(result.state),
          "nodes"_a = result.nodes,
          "t_solver"_a =
              std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
                  .count(),
          "callback"_a = result.totalCallback,
          "t_callback"_a = result.totalCallbackTime,
          "lazy"_a = result.totalLazy));
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
