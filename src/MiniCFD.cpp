#include <iostream>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include "Simulation.h"

int main(int argc, char** argv) {
    // Command line parsing
    cxxopts::Options options("MiniCFD",
                             "Simple CFD simulation for didactic purposes.");

    // clang-format off
    options.add_options()
        ("l,log-level", "Log level (trace, debug, info, warn, err, critical or off)", cxxopts::value<std::string>()->default_value("info"))
        ("d,domain-size", "Number of the simulation cells along all three axes", cxxopts::value<size_t>()->default_value("10"))
        ("c,cell-size", "Size of each simulation cell", cxxopts::value<ScalarT>()->default_value("1.0"))
        ("e,end-time", "Simulation duration (seconds)", cxxopts::value<double>()->default_value("1.0"))
        ("s,step-size", "Simulation step size (seconds)", cxxopts::value<double>()->default_value("0.05"))
        ("o,output-prefix", "Output file prefix", cxxopts::value<std::string>()->default_value("fields"))
        ("p,preconditioner", "Preconditioner type (none, jacobi, dic)", cxxopts::value<std::string>()->default_value("dic"))
        ("h,help", "Print usage")
    ;
    // clang-format on

    auto args = options.parse(argc, argv);
    if (args.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // read log level from CLI options
    const std::string logLevel = args["log-level"].as<std::string>();
    if (logLevel == "trace") {
        spdlog::set_level(spdlog::level::trace);
    } else if (logLevel == "debug") {
        spdlog::set_level(spdlog::level::debug);
    } else if (logLevel == "info") {
        spdlog::set_level(spdlog::level::info);
    } else if (logLevel == "warn") {
        spdlog::set_level(spdlog::level::warn);
    } else if (logLevel == "err") {
        spdlog::set_level(spdlog::level::err);
    } else if (logLevel == "critical") {
        spdlog::set_level(spdlog::level::critical);
    } else if (logLevel == "off") {
        spdlog::set_level(spdlog::level::off);
    } else {
        spdlog::critical("Log level '{}' not recognized!");
        exit(-1);
    }

    spdlog::info("Welcome to MiniCFD!");

    // read preconditioner type from CLI options
    PreconditionerType preconditionerType{PreconditionerType::DIC};
    const std::string precondTypeStr = args["preconditioner"].as<std::string>();
    if (precondTypeStr == "none") {
        preconditionerType = PreconditionerType::NONE;
    } else if (precondTypeStr == "jacobi") {
        preconditionerType = PreconditionerType::JACOBI;
    } else if (precondTypeStr == "dic") {
        preconditionerType = PreconditionerType::DIC;
    } else {
        spdlog::critical("Invalid preconditioner '{}'", precondTypeStr);
        exit(-1);
    }

    // set up simulation config
    auto N = args["domain-size"].as<size_t>();
    const SimulationConfig cfg = {
        .width = N,
        .height = N,
        .depth = N,
        .endTime = args["end-time"].as<double>(),
        .stepSize = args["step-size"].as<double>(),
        .cellSize = args["cell-size"].as<ScalarT>(),
        .outputPrefix = args["output-prefix"].as<std::string>(),
        .preconditionerType = preconditionerType,
        .disableFileOutput = false,
    };

    // run simulation
    simulate(cfg);

    return 0;
}
