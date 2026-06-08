<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

The Heat Equation App is a Python-based solver for one-dimensional heat diffusion problems. The program provides an interactive Graphical User Interface (GUI) that allows users to choose between two numerical methods, specify the boundary conditions, and define the physical and numerical parameters of the problem, including the spatial domain, simulation time, and discretization steps.

The available numerical methods are the Forward-Time Centered-Space (FTCS) scheme and the Crank–Nicolson method. Boundary conditions can be specified as Dirichlet, Neumann, or Robin, with support for both steady-state and time-dependent conditions. Once the problem is defined, the numerical solution is computed and visualized through an animated representation of the temperature distribution.

</div>

## Project Structure

The application is divided into two modules:

- **Mini_1D_Heat_App.py**: Tkinter-based GUI responsible for user interaction and problem setup.
- **HeatLib.py**: Numerical library implementing the FTCS and Crank–Nicolson schemes, boundary-condition treatment, the Thomas algorithm for tridiagonal systems, and solution visualization.

The GUI collects the user-defined parameters and passes them to HeatLib, where the numerical solution of the heat equation is performed.

Mathematical Model

The Heat Equation App solves the one-dimensional heat equation

a\frac{\partial^2 u}{\partial x^2},
]

where (u(x,t)) denotes the temperature distribution and (a) is the thermal diffusivity. The boundary conditions are written in the general form

[
k\frac{\partial u}{\partial x}=f-au,
]

which allows the implementation of Dirichlet, Neumann, and Robin boundary conditions by appropriate choices of the coefficients (a) and (k).

Numerical Methods

Two finite-difference schemes are implemented.

Forward-Time Centered-Space (FTCS): an explicit method combining a forward Euler discretization in time with a second-order central difference approximation in space. The spatial discretization is implemented using NumPy slicing for computational efficiency.
Crank–Nicolson: an implicit second-order accurate method obtained by averaging the spatial operator between consecutive time levels. At each time step, the method leads to the solution of a tridiagonal linear system.
Tridiagonal System

The Crank–Nicolson discretization produces, for each node,

Q_i,
]

where

(A_W) denotes the lower diagonal coefficients,
(A_P) the main diagonal coefficients,
(A_E) the upper diagonal coefficients,
(Q) the forcing vector.

Consequently, the resulting coefficient matrix is tridiagonal.

Thomas Algorithm

The tridiagonal system is solved using the Thomas algorithm, which is a specialized form of Gaussian elimination for tridiagonal matrices. Instead of storing and manipulating the full coefficient matrix, the algorithm operates directly on the three diagonals.

During the forward sweep, modified coefficients are computed recursively,

\frac{A_E^i}
{A_P^i-A_W^iA_E^{*(i-1)}},
]

and

\frac{
Q_i-A_W^iQ_{i-1}^*
}
{
A_P^i-A_W^iA_E^{*(i-1)}
}.
]

The first row provides the initial values,

\frac{Q_0}{A_P^0}.
]

Once the forward sweep is completed, the solution is recovered through backward substitution,

[
\Phi_N=Q_N^*,
]

and

Q_i^*

A_E^{*i}\Phi_{i+1},
]

for (i=N-1,\ldots,0).

In the implementation, the arrays Ae_prime and Q_prime correspond to the modified upper diagonal coefficients and modified forcing terms generated during the forward sweep, while the backward sweep reconstructs the solution vector. Since the Thomas algorithm exploits the tridiagonal structure of the Crank–Nicolson matrix, it requires only (O(N)) operations and significantly reduces the computational cost compared with general-purpose linear system solvers.



