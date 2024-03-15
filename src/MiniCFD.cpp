
#include "IO.h"

int main(int argc, char** argv) {
    // Command line parsing (later)

    auto pair = read_from_file<double>("sample_io");

    // Set up grid

    // Time loop
    // - Solve advection: implicit euler
    // - Pressure corection
    // - Write field (later)

    // - Validierung

    // - Visualisierung

    return 0;
}
