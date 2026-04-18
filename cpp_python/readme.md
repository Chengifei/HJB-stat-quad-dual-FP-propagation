# HJB-Stat-Quad-Dual Solver

This is a solver for Hamilton-Jacobi-Bellman equations of statical control problems. The solver is especially efficient for problems with low dimensional nonlinearities, as it exploits the domain of nonlinearity and converts the control problem into an equivalent stat-quad dual formulation.

The core loops are written in C++ for performance using Eigen3 and Boost, and it exports native bindings for Python and MATLAB (deprecated).


## Architecture

The implementation and interface are modular and consist of the following components:

1. **Problem Definition:** The state space, coefficients, and nonlinearity specifications (e.g., $l(x)$ and $f(x)$) are encapsulated in the `Ntilde` and `PDE_info` objects.
2. **Numerical Integration:** System dynamics (state/costate trajectories) are functors compatible with `boost::numeric::odeint`.
3. **Iterative Sweeps:** The final staticization is solved by iteratively solving coupled $p$ and $q$ systems using custom Runge-Kutta 4 steppers and interpolating state transformation matrices ($R$ and $S$, Bernoulli substitution).

Files:
* **`pde_info.h`**: The central data structure holding the PDE coefficients ($A, C, M_0, L_0, \hat{C}, \Gamma$) and the `Ntilde` class.
* **`N_tilde.cpp`**: Implements the nonlinear mappings and their exact first and second-order derivatives (Hessians).
* **`solver.h`** & **`solver.cpp`**: Orchestrates the time grids, executes ODE integrations, handles error tolerance loops, and reconstructs the value/cost.
* **`ode.h`**: Defines the ODE dynamics and a custom `pre_discretized_rk4_stepper`.
* **`solvers.h`**: Implements templated iterative solvers (fixed-point iteration and Newton-Raphson).
* **`transform.h`**: Manages state transformations interpolations by precomputing exponential maps of the Hamiltonian blocks.

### Interfaces & Examples
* **`wrapper.cpp`**: Bridges C++ to external environments. Compiles as a MATLAB MEX function or a Python extension module depending on the build flags.
* **`test.py`**: A Python testing and visualization driver that instantiates the solver.

## Dependencies

* **C++20** (only tested with gcc-14 and libstdc++)
* **CMake** (>= 3.10)
* **Eigen3**
* **Boost** (`boost::numeric::odeint`)
* **Python 3 & NumPy** (For Python bindings)
* **MATLAB** (for MEX bindings)

## Build Instructions

To build the Python shared object library (`HJB_solve.so`), run the following commands:

```bash
mkdir build
cd build
cmake ..
make
```
