# HJB solver via stat-quad duality
This is a MATLAB implementation of the numerical procedure for solving HJB PDEs
equations associated with mildly nonlinear optimal control problems, especially
for systems where the nonlinearity is defined on lower dimensional
subvarieties.

## Project structure

### Main files, interface
* `test.m`: the primary driver script. Defines the test problem, calls the solver and does some plotting.
* `main.m`: The main solver routine. Caches shared variable and dispatches the solver.

### Utilities
* `eqn_P_Psi.m`: Defines the Riccati equation for $P$.
* `D_sigma_s`: Computes the integration kernel $\breve{\mathcal{D}}_{\cdot, \cdot}$.
* `propalpnlpsiTbndim.m`: Propagates the $\alpha$ process.
* `Getalphafuncnl.m`: Implements the fixed-point iteration for $\alpha$.
* `newtonforeta.m`: Evaluates $\eta^{-1}$ using Newton's method. (This is only for speed)
* `etainverse.m`: Evaluates $\eta^{-1}$ using fixed point iterations.
* `qrintegrals.m`: Performs the integrals involved for $q$ and $r$ processes given $\alpha$.
* `gradtheta.m`/`gettheta.m`/`gradeta.m`: Wrappers for $\nabla \Theta$, $\Theta$, $\nabla \eta$.
