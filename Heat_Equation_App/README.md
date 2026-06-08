<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

The Heat Equation App is a Python-based solver for one-dimensional heat diffusion problems. The program provides an interactive Graphical User Interface (GUI) that allows users to choose between two numerical methods, specify the boundary conditions, and define the physical and numerical parameters of the problem, including the spatial domain, simulation time, and discretization steps.

The available numerical methods are the Forward-Time Centered-Space (FTCS) scheme and the Crank–Nicolson method. Boundary conditions can be specified as Dirichlet, Neumann, or Robin, with support for both steady-state and time-dependent conditions. Once the problem is defined, the numerical solution is computed and visualized through an animated representation of the temperature distribution.

</div>

## Project Structure

The application is divided into two modules:

- **Mini_1D_Heat_App.py**: Tkinter-based GUI responsible for user interaction and problem setup.
- **HeatLib.py**: Numerical library implementing the FTCS and Crank–Nicolson schemes, boundary-condition treatment, the Thomas algorithm for tridiagonal systems, and solution visualization.

The GUI collects the user-defined parameters and passes them to HeatLib, where the numerical solution of the heat equation is performed.

## Mathematical Model

- The first method provides the simple Second order Centered-space scheme with explict Euler time discritization, written via Numpy Slicing.

- The second method it is the Crank-Nicolson where the Tridiagonal System is sloved with Thomas algorithm.

Both shemes are defined in the begining of the script.


