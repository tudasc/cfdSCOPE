# MiniCFD

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
./minicfd -h # print help about command-line arguments
```

### Visualization runs 
The visualization script tends to take a while for larger simulations.
Thus, in order to just try out the simulation, use the call:
```bash
OMP_NUM_THREADS=8 ./minicfd -d 20 -e 5 -s 0.4
```

### Benchmarking runs
To conduct performance measurements, use a larger simulation setup:
```bash
time OMP_NUM_THREADS=16 ./minicfd -d 100 -e 6 -s 0.4
```

## Visualization
The `visualize.py` script can visualize the force field with an animation. Refer to `python visualize.py --help` for usage instructions.

### Python environment
Steps necessary to set up the Python virtual environment for the visualization (or refer to `requirements.txt`):
```bash
# setup venv
python -m venv venv
source venv/bin/activate
# install required packages
pip install --upgrade pip
pip install numpy scipy matplotlib tqdm
```

### Run visualization
```bash
python ../visualize.py ./ --frames_per_step 12 --time_step=0.4 --num_particles=1000
```
Make sure to match the `--time_step` parameter to the step size (`-s`) from the simulation call.
By default, the script outputs a file called `visu.html` and a directory `visu_frames` that contains all the frames.

## Tests
To run the tests, first enable their compilation at configure time by passing `-DBUILD_TESTS=ON` to CMake.
Then, run:
```bash
tests/unittests
tests/integrationtests
```
