#ifndef CONFIG_H
#define CONFIG_H

#include <cstddef>
#include <string>

using ScalarT = double;

enum class PreconditionerType { NONE, JACOBI, DIC };

enum class FileFormat {CSV, RAW};

struct SimulationConfig {
    size_t width, height, depth;
    double endTime, stepSize;
    ScalarT cellSize;
    ScalarT lidSpeed;
    std::string outputPrefix;
    PreconditionerType preconditionerType;
    FileFormat fileFormat;
    bool disableFileOutput;
};

#endif // CONFIG_H
