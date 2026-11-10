<h2 align="center">Heat Equation</h2>

<div style="text-align: justify; text-justify: inter-word; font-size: 16px;">

<p>
This repository contains educational and scientific Python projects developed within the framework of Numerical Methods and Computational Physics. The projects focus on heat transfer and diffusion phenomena, combining classical numerical techniques with scientific computing and interactive visualization. Their objective is to provide both practical simulation tools and accessible implementations of partial differential equation (PDE) solvers.
</p>

<p>
The repository currently consists of three independent projects:
</p>

</div>

---

<h3>Mini 1D Heat Equation App</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 16px;">

<p>
The <b>Mini 1D Heat Equation App</b> is an interactive Python application for solving one-dimensional heat diffusion problems. The project combines a Tkinter-based Graphical User Interface (GUI) with a dedicated numerical library implementing the underlying finite-difference algorithms.
</p>

<p>
Main features include:
</p>

<ul>

<li>Forward-Time Centered-Space (FTCS) discretization.</li>

<li>Crank–Nicolson implicit discretization.</li>

<li>Thomas algorithm for tridiagonal linear systems.</li>

<li>Dirichlet, Neumann, and Robin boundary conditions.</li>

<li>Support for steady-state and time-dependent boundary data.</li>

<li>Animated visualization of the temperature field.</li>

</ul>

<p>
The project serves both as a practical heat equation solver and as an educational implementation of classical finite-difference methods.
</p>

</div>

---

<h3>Initial-Boundary Compatibility, Ghost Points and Animation</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 16px;">

<p>
The <b>Initial-Boundary Compatibility, Ghost Points and Animation</b> notebook is an educational notebook designed to demonstrate three important aspects of numerical methods through two one-dimensional examples of the heat equation: compatibility between initial and boundary conditions, the use of ghost points to impose boundary conditions, and the visualization of numerical solutions through animations.
</p>

<ul>

<li><b>Compatibility:</b> Demonstration of compatible and incompatible initial and boundary conditions, along with a proposed treatment for incompatible cases. Initial-boundary compatibility concerns the regularity of the solution and should not be confused with the consistency of a numerical scheme (sometimes called compatibility) or the well-posedness of the mathematical problem.</li>

<li><b>Ghost Points:</b> Implementation of ghost points for imposing boundary conditions using a simple finite volume method (FVM), where ghost-cell values are explicitly stored in the solution array.</li>

<li><b>Animation:</b> A visualization approach that accounts for the gaps between the boundary cell centers and the physical boundaries in a cell-centered FVM grid.</li>

</ul>

<p>
The notebook provides a practical introduction to these numerical techniques, emphasizing their mathematical interpretation, implementation, and visualization through simple heat diffusion examples.
</p>

</div>

---

<h3>Stochastic Heat Equation</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 16px;">

<p>
The <b>Stochastic Heat Equation</b> is a two-dimensional advection-diffusion simulation originally developed as a numerical methods coursework project and progressively extended into a more comprehensive computational model.
</p>

<p>
The project incorporates:
</p>

<ul>

<li>Multiple localized Gaussian heat sources.</li>

<li>Stochastic source positioning through custom probability distributions.</li>

<li>Discrete-time Poisson source activation.</li>

<li>Finite Volume discretization.</li>

<li>Ghost-cell treatment of boundary conditions through algebraic elimination of ghost values from the discretized equations.</li>

<li>Upwind and Central Differencing schemes.</li>

<li>Numba JIT and parallel acceleration.</li>

<li>Interactive Plotly visualization and animation.</li>

</ul>

<p>
The notebook combines stochastic processes, numerical discretization techniques, and scientific visualization within a unified educational framework.
</p>

</div>

---

<h3>Future Development</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 16px;">

<p>
The long-term objective of this repository is to progressively incorporate more advanced numerical techniques for heat-transfer and diffusion problems while preserving the educational character of the projects.
</p>

<p>
Planned developments include:
</p>

<ul>

<li>Support for mixed boundary conditions in the one-dimensional solver.</li>

<li>Extension to variable thermal diffusivity models.</li>

<li>Implementation and comparison of advanced tridiagonal and sparse linear-system solvers, including Cyclic Reduction and related algorithms.</li>

</ul>

</div>
