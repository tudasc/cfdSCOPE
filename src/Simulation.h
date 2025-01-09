#ifndef SIMULATION_H
#define SIMULATION_H

#include "Grid.h"
#include <cstddef>
#include <memory>
#include <string>

using ScalarT = double;

enum class PreconditionerType { NONE, JACOBI, DIC };

struct SimulationConfig {
    size_t width, height, depth;
    double endTime, stepSize;
    ScalarT cellSize;
    ScalarT lidSpeed;
    std::string outputPrefix;
    PreconditionerType preconditionerType;
    bool disableFileOutput;
};

struct SimulationOutput {
    std::unique_ptr<VelocityField<ScalarT>> velocityField;
    std::unique_ptr<PressureField<ScalarT>> pressureField;
};

/**
    Main simulation routine
*/
SimulationOutput simulate(const SimulationConfig& cfg);

#endif
