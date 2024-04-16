#ifndef SIMULATION_H
#define SIMULATION_H

#include <cstddef>
#include <string>

using ScalarT = double;

enum class PreconditionerType { NONE, JACOBI, DIC };

struct SimulationConfig {
    size_t width, height, depth;
    double endTime, stepSize;
    ScalarT cellSize;
    std::string outputPrefix;
    PreconditionerType preconditionerType;
};

/**
    Main simulation routine
*/
void simulate(const SimulationConfig& cfg);

#endif
