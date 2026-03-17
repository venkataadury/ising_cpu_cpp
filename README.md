# ising_cpu_cpp
Rapidly simulate a 2D ising model. Useful to quickly run and validate tests. Saves the entire trajectory to simple text.

## Running the program
Two simple scripts are provided to run the Ising Model simulation (on a single thread - safe but slow, multithreaded - may rarely crash due to race conditions. You can rerun it)\\
``> simulate_ising.sh``
``> simulate_ising_parallel.sh``

## Editing parameters
You can easily edit the following parameters by changing the heaader (*#define lines*) in the two `*.cpp` files
1. Jx (Default: 1) - Coupling constant in the x-direction
2. Jy (Default: 1) - Coupling constant in the y-direction
3. GRID\_DIM (Default: 8) - NxN grid is simulated. This is the value of N

