# cfdSCOPE

## Build instructions
```bash
mkdir build
cd build
cmake ..
make
```
To enable trace logging, pass `-DENABLE_TRACE_LOG=ON` to CMake. (Note, that you also need to set the correct logging level when calling the executable.)

## Run instructions
```bash
./cfdscope -h # print help about command-line arguments
```

### Visualization runs 
A good size for trying out the simulation and generating data for visualization is 64x64x64 grid cells.
```bash
OMP_NUM_THREADS=8 ./cfdscope -d 64 -e 5 -s 0.4
```

### Benchmarking runs
To conduct performance measurements, use a larger simulation setup:
```bash
time OMP_NUM_THREADS=16 ./cfdscope -d 100 -e 6 -s 0.4
```

## Visualization
The `visualize.py` script can be used to render a specific frame of the simulation with [Paraview](https://www.paraview.org/).
To use it, first install Paraview via your system's package manager.
See `python visualize.py -h` for available parameters.

### Run visualization
```bash
pvpython visualize.py "<simulation_dir>/fields.csv.*" --frame 20 --size 64
```
Make sure to match the `-size` parameter to the grid size (`-d`) from the simulation call.

## Tests
To run the tests, first enable their compilation at configure time by passing `-DBUILD_TESTS=ON` to CMake.
Then, run:
```bash
tests/unittests
tests/integrationtests
```
