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
./minicfd --log-level debug --domain-size 5 --end-time 2.0
```

## Visualization
The `visualize.py` script can visualize the force field with an animation. Refer to `python visualize.py --help` for usage instructions.

### Python environment
Steps necessary to set up the python environment for the visualization (or refer to `requirements.txt`):
```bash
# setup venv
python -m venv venv
source venv/bin/activate
# install required packages
pip install --upgrade pip
pip install numpy scipy matplotlib tqdm
```
