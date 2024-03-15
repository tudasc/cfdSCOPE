
#include "IO.h"

int main(int argc, char** argv) {
    // Command line parsing (later)

    auto pair = read_from_file<double>("/home/tim/peng-prak-code-mini-cfd/sample_io");
    write_to_file(pair.first,pair.second,"/home/tim/peng-prak-code-mini-cfd/sample_io_out");

    // Set up grid

    // Time loop
    // - Solve advection: implicit euler
    // - Pressure corection
    // - Write field (later)

    // - Validierung

    // - Visualisierung

    return 0;
}
