<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

The Heat Equation App is a Python-based solver for one-dimensional heat diffusion problems. The program provides an interactive Graphical User Interface (GUI) that allows users to choose between two numerical methods, specify the boundary conditions, and define the physical and numerical parameters of the problem, including the spatial domain, simulation time, and discretization steps.

The available numerical methods are the Forward-Time Centered-Space (FTCS) scheme and the Crank–Nicolson method. Boundary conditions can be specified as Dirichlet, Neumann, or Robin, with support for both steady-state and time-dependent conditions. Once the problem is defined, the numerical solution is computed and visualized through an animated representation of the temperature distribution.

</div>

## Project Structure

The application is divided into two modules:

- **Mini_1D_Heat_App.py**: Tkinter-based GUI responsible for user interaction and problem setup.
- **HeatLib.py**: Numerical library implementing the FTCS and Crank–Nicolson schemes, boundary-condition treatment, the Thomas algorithm for tridiagonal systems, and solution visualization.

The GUI collects the user-defined parameters and passes them to HeatLib, where the numerical solution of the heat equation is performed.

<h2>Mathematical Model</h2>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

<p>
The Heat Equation App solves the one-dimensional heat equation
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\frac{\partial&space;u}{\partial&space;t}=a_{T}\frac{\partial^2&space;u}{\partial&space;x^2}" />
</p>

<p>
where <code>u(x,t)</code> denotes the temperature distribution and <code>a<sub>T</sub></code> is the thermal diffusivity. The boundary conditions are written in the general form
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?k\frac{\partial&space;u}{\partial&space;x}=f-au" />
</p>

<p>
which allows the implementation of Dirichlet, Neumann, and Robin boundary conditions by appropriate choices of the coefficients <code>a</code>, <code>k</code> and <code>f</code>.
</p>

<h3>Numerical Methods</h3>

<p>
Two finite-difference schemes are implemented:
</p>

<ul>
  <li>
    <strong>Forward-Time Centered-Space (FTCS):</strong>
    an explicit method combining a forward Euler discretization in time with a second-order central difference approximation in space. The spatial discretization is implemented using NumPy slicing for computational efficiency.
  </li>
  <li>
    <strong>Crank–Nicolson:</strong>
    an implicit second-order accurate method obtained by averaging the spatial operator between consecutive time levels. At each time step, the method leads to the solution of a tridiagonal linear system.
  </li>
</ul>

<h3>Tridiagonal System</h3>

<p>
The Crank–Nicolson discretization produces, for each node,
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?A_W^i\Phi_{i-1}+A_P^i\Phi_i+A_E^i\Phi_{i+1}=Q_i" />
</p>

<p>
where <code>A<sub>W</sub></code> denotes the lower diagonal coefficient, <code>A<sub>P</sub></code> the main diagonal coefficient, <code>A<sub>E</sub></code> the upper diagonal coefficient, and <code>Q</code> the forcing vector. Consequently, the coefficient matrix is tridiagonal.
</p>

<h3>Thomas Algorithm</h3>

<p>
The tridiagonal system is solved using the Thomas algorithm, a specialized form of Gaussian elimination for tridiagonal matrices. Instead of storing and manipulating the full coefficient matrix, the algorithm operates directly on the three diagonals.
</p>

<p>
During the forward elimination step, the main diagonal and forcing term are modified recursively:
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?A_P^{*i}=A_P^i-\frac{A_W^iA_E^{i-1}}{A_P^{*(i-1)}}" />
  &nbsp;&nbsp;&nbsp;
  which corresponds to <code>Ap_star</code>.
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?Q_i^*=Q_i-\frac{A_W^iQ_{i-1}^*}{A_P^{*(i-1)}}" />
  &nbsp;&nbsp;&nbsp;
  which corresponds to <code>Q_star</code>.
</p>

<p>
The first row is not modified:
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?{A_P^{0*}=A_P^0},\qquad{Q_0^*=Q_0}" />
</p>

<p>
Once the forward elimination is complete, the solution is recovered through backward substitution:
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\Phi_N=\frac{Q_N^*}{A_P^{*N}}" />
</p>

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\Phi_i=\frac{Q_i^*-A_E^i\Phi_{i+1}}{A_P^{i*}},\qquad{i=N-1,\ldots,0}" />
</p>

<p>
In the implementation, <code>Ap_star</code> stores the modified main diagonal and <code>Q_star</code> stores the modified forcing vector generated during forward elimination. The backward substitution then reconstructs the solution vector from the last node to the first.
</p>

<h3>Boundary Condition Handling</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

<p>
The Heat Equation App employs the general linear boundary condition
</p>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?k\frac{\partial&space;u}{\partial&space;x}=f-au"/>
</p>

<p>
where <code>a</code>, <code>k</code>, and <code>f</code> may depend on time. This formulation allows the three most common boundary-condition types to be treated within a unified framework.
</p>

<ul>
<li>
<strong>Dirichlet boundary conditions</strong> are obtained by setting
<code>k=0</code>, giving
</li>
</ul>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?u=\frac{f}{a}"/>
</p>

<ul>
<li>
<strong>Neumann boundary conditions</strong> are obtained by setting
<code>a=0</code>, yielding
</li>
</ul>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?k\frac{\partial&space;u}{\partial&space;x}=f"/>
</p>

<ul>
<li>
<strong>Robin boundary conditions</strong> correspond to the general case where both coefficients are non-zero:
</li>
</ul>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?k\frac{\partial&space;u}{\partial&space;x}=f-au"/>
</p>

<p>
The program automatically determines the appropriate numerical treatment by examining the boundary coefficients supplied by the user:
</p>

<ul>
<li><code>k=0</code> and <code>a≠0</code> &rarr; Dirichlet.</li>
<li><code>a=0</code> and <code>k≠0</code> &rarr; Neumann.</li>
<li><code>a≠0</code> and <code>k≠0</code> &rarr; Robin.</li>
</ul>

<p>
This unified implementation avoids the need for separate solvers for different boundary-condition types and naturally accommodates both steady-state and time-dependent boundary data. At present, mixed boundary-condition configurations, such as Dirichlet on one side and Neumann or Robin on the other, are not yet implemented. This is a planned extension of the boundary-condition handling.
</p>

</div>
