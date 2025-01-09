#ifndef SIMULATION_H
#define SIMULATION_H

#include <cstddef>
#include <memory>
#include <string>

#include "Config.h"
#include "Grid.h"

struct SimulationOutput {
    std::unique_ptr<VelocityField<ScalarT>> velocityField;
    std::unique_ptr<PressureField<ScalarT>> pressureField;
};

/**
    Main simulation routine
*/
SimulationOutput simulate(const SimulationConfig& cfg);

#endif
