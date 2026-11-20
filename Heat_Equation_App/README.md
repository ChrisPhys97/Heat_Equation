
# Heat Equation App

The **Heat Equation App** is a Python-based numerical solver for one-dimensional heat diffusion problems. The application provides an interactive Graphical User Interface (GUI) that allows users to select a numerical method, define the initial temperature distribution, specify boundary conditions, and configure the physical and numerical parameters of the problem.

Two finite-difference methods are implemented: the **Forward-Time Centered-Space (FTCS)** scheme and the **Crank–Nicolson method**. Both methods support Dirichlet, Neumann, and Robin boundary conditions, including mixed configurations and time-dependent boundary forcing.

An important feature of the application is the separate specification of the initial boundary conditions and the boundary conditions imposed during the subsequent evolution. The user is responsible for calculating the initial boundary values so that they are consistent with the prescribed initial temperature distribution.

Once the problem is defined, the numerical solution is computed and visualized through an animated representation of the temperature distribution.

## Project Structure

The application is divided into two Python modules:

- **`Mini_1D_Heat_App.py`** — Tkinter-based GUI responsible for user interaction, parameter input, and numerical method selection.
- **`HeatLib.py`** — Numerical library implementing the FTCS and Crank–Nicolson schemes, boundary-condition treatment, the Thomas algorithm for tridiagonal systems, and solution visualization.

The GUI collects the user-defined parameters and passes them to the numerical routines in `HeatLib.py`. The numerical library can also be used independently of the graphical interface.

## Mathematical Model

The application solves the one-dimensional heat equation

$$
\frac{\partial u}{\partial t}=\alpha\frac{\partial^2u}{\partial x^2}
$$

where $u(x,t)$ denotes the temperature distribution and $\alpha$ is the constant thermal diffusivity.

The spatial domain and simulation interval are defined by

$$
0\leq x\leq L,\qquad 0\leq t\leq T.
$$

The initial condition is prescribed as

$$
u(x,0)=\phi(x).
$$

The general boundary condition is written as

$$
a_bu+k_b\frac{\partial u}{\partial n}=f_b(t),
$$

where $a_b$ and $k_b$ are boundary coefficients, $f_b(t)$ is the prescribed boundary forcing, and $\partial u/\partial n$ denotes the derivative along the outward normal.

For a one-dimensional domain, the outward-normal derivatives are

$$
\left.\frac{\partial u}{\partial n}\right|_{x=0}=-\left.\frac{\partial u}{\partial x}\right|_{x=0}
$$

and

$$
\left.\frac{\partial u}{\partial n}\right|_{x=L}=+\left.\frac{\partial u}{\partial x}\right|_{x=L}.
$$

This convention allows the same boundary-update formula to be used at both ends of the domain.

## Numerical Methods

### Forward-Time Centered-Space (FTCS)

The FTCS scheme combines a forward Euler approximation in time with a second-order central difference approximation in space.

For an interior node $i$, the discretized equation is

$$
u_i^{n+1}=u_i^n+r\left(u_{i+1}^n-2u_i^n+u_{i-1}^n\right),
$$

where

$$
r=\frac{\alpha\Delta t}{\Delta x^2}.
$$

The spatial discretization is implemented using NumPy slicing, avoiding explicit loops over the interior nodes.

For standard interior-node stability, the FTCS scheme requires

$$
r\leq\frac12.
$$

Neumann and Robin boundary conditions may introduce additional restrictions. The implementation checks a sufficient stability condition based on the boundary coefficients when their ratios are nonnegative.

### Crank–Nicolson Method

The Crank–Nicolson method is an implicit finite-difference scheme obtained by averaging the spatial diffusion operator between consecutive time levels.

For interior nodes,

$$
\begin{aligned}-r u_{i-1}^{n+1}+(2+2r)u_i^{n+1}-r u_{i+1}^{n+1}={}&r u_{i-1}^{n}\\&+(2-2r)u_i^n+r u_{i+1}^n.\end{aligned}
$$

The method is second-order accurate in time and space for sufficiently smooth solutions and boundary data.

Unlike FTCS, Crank–Nicolson is not subject to the standard explicit diffusion stability restriction. However, its numerical accuracy still depends on the chosen spatial and temporal resolutions.

At every time step, the discretization produces a tridiagonal linear system.

## Tridiagonal System

The Crank–Nicolson discretization produces equations of the form

$$
A_W^i u_{i-1}^{n+1}+A_P^i u_i^{n+1}+A_E^i u_{i+1}^{n+1}=Q_i,
$$

where:

- $A_W$ represents the lower diagonal coefficient.
- $A_P$ represents the main diagonal coefficient.
- $A_E$ represents the upper diagonal coefficient.
- $Q$ represents the right-hand-side vector.

The coefficient matrix is therefore tridiagonal.

The current implementation constructs the three diagonals directly and solves the resulting system without assembling a full dense matrix.

Boundary conditions are incorporated into the first and last equations of the system, allowing Dirichlet, Neumann, Robin, and mixed boundary-condition configurations.

## Thomas Algorithm

The tridiagonal system is solved using the **Thomas algorithm**, a specialized form of Gaussian elimination for tridiagonal matrices.

Instead of storing and manipulating the full coefficient matrix, the algorithm operates directly on its three diagonals.

### Forward Elimination

The modified main diagonal is calculated recursively:

$$
A_P^{\ast i}=A_P^i-\frac{A_W^i A_E^{i-1}}{A_P^{\ast(i-1)}}.
$$

The modified right-hand side is

$$
Q_i^{\ast}={Q_i}-\frac{A_W^{i}Q_{i-1}^{\ast}}{A_P^{*(i-1)}}.
$$

The first row remains unchanged:

$$
A_P^{\ast 0}=A_P^0,\qquad Q_0^{\ast}=Q_0.
$$

### Backward Substitution

After forward elimination, the final unknown is obtained from

$$
u_N=\frac{Q_N^{\ast}}{A_P^{\ast N}}.
$$

The remaining values are calculated recursively:

$$
u_i=\frac{Q_i^{\ast}-A_E^i u_{i+1}}{A_P^{\ast i}},\qquad i=N-1,\ldots,0.
$$

In `HeatLib.py`, `Ap_star` stores the modified main diagonal and `Q_star` stores the modified right-hand-side vector.

The Thomas algorithm has linear computational complexity, making it efficient for one-dimensional diffusion problems.

## Boundary Condition Handling

The application supports three types of boundary conditions through the general formulation

$$
a_bu+k_b\frac{\partial u}{\partial n}=f_b(t).
$$

### Dirichlet Boundary Conditions

Dirichlet conditions prescribe the boundary temperature.

Setting $k_b=0$ gives

$$
u_b(t)=\frac{f_b(t)}{a_b},\qquad a_b\neq0.
$$

The prescribed temperature is assigned directly to the boundary node.

### Neumann Boundary Conditions

Neumann conditions prescribe the outward-normal temperature gradient.

Setting $a_b=0$ gives

$$
\frac{\partial u}{\partial n}=\frac{f_b(t)}{k_b},\qquad k_b\neq0.
$$

A second-order central difference with ghost nodes is used to incorporate the boundary gradient into the numerical scheme.

### Robin Boundary Conditions

Robin conditions combine the boundary temperature and its outward-normal gradient:

$$
a_bu+k_b\frac{\partial u}{\partial n}=f_b(t),\qquad a_b\neq0,\quad k_b\neq0.
$$

As with Neumann conditions, ghost-node approximations are used to derive the boundary equations.

### Initial and Prescribed Boundary Values

The user must provide **two boundary values for each side of the domain**:

1. The initial boundary value at $n=0$, calculated from the initial temperature distribution.
2. The new prescribed boundary value or expression for $n\geq1$.

For the initial temperature distribution

$$
u(x,0)=\phi(x),
$$

the required initial boundary values are:

| Condition | Initial boundary value |
|---|---|
| Dirichlet | $f_b^0=a_b\phi(x_b)$ |
| Neumann | $f_b^0=k_b\,\partial\phi/\partial n$ |
| Robin | $f_b^0=a_b\phi(x_b)+k_b\,\partial\phi/\partial n$ |

The initial values are supplied directly by the user. HeatLib does not calculate them automatically.

The boundary forcing arrays are constructed as


$$
\begin{cases}
f_b^n = f_{b,\mathrm{initial}}, & n = 0 \\
\rule{0pt}{24pt}f_b^n = f_{b,\mathrm{new}}(t_n), & n \geq 1
\end{cases}
$$


The initial temperature profile is preserved in `u[0, :]` without modification.

### Numerical Implementation of Boundary Conditions: FTCS

The treatment of boundary conditions is an important part of the numerical solution because the standard central-difference stencil cannot be applied directly at the domain boundaries without additional information. HeatLib handles the boundary nodes separately from the interior nodes, using the general boundary-condition formulation

$$
a_bu+k_b\frac{\partial u}{\partial n}=f_b(t).
$$

The coefficients $a_b$ and $k_b$ determine the type of boundary condition, while $f_b(t)$ represents the prescribed boundary forcing. The latter may be constant or time-dependent and is evaluated at the discrete time levels of the simulation.

#### Boundary Forcing Arrays

For each boundary, the user specifies an initial forcing value and a new forcing expression. These are combined into a time-dependent array using the `boundary_forcing()` function:

```python
def boundary_forcing(t, initial_value, new_expression):

    f = evaluate_expression(new_expression, "t", t)

    f[0] = float(initial_value)

    return f
```

The resulting array is defined by

$$
f_b^n =
\begin{cases}
f_{b,\mathrm{initial}}, & n=0,\\
f_{b,\mathrm{new}}(t_n), & n\geq1.
\end{cases}
$$

The initial value is provided by the user and should be consistent with the initial temperature distribution. The new expression determines the boundary forcing during the subsequent simulation.

The forcing $f_b$ does not necessarily represent a temperature. Depending on the boundary coefficients, it may prescribe a boundary temperature, a temperature gradient, or a linear combination of both.

#### Dirichlet Boundary Conditions

For Dirichlet conditions, $k_b=0$, and the general boundary equation reduces to

$$
u_b(t)=\frac{f_b(t)}{a_b}.
$$

Since the boundary temperature is explicitly prescribed, it is assigned directly at the new time level:

$$
u_b^{n+1}=\frac{f_b^{n+1}}{a_b}.
$$

In HeatLib, the left boundary is updated as

```python
if k0 == 0:
    u[n+1, 0] = f0[n+1] / a0
```

The same approach is used at the right boundary.

Notice that the initial temperature distribution is not overwritten. The first row of the solution array, `u[0, :]`, retains the original user-defined profile, while the new Dirichlet conditions are applied from the first subsequent time level.

#### Neumann and Robin Boundary Conditions

For Neumann and Robin conditions, $k_b\neq0$, and the boundary derivative must be approximated numerically.

At the left boundary, the outward-normal derivative is

$$
\left.\frac{\partial u}{\partial n}\right|_{x=0}=-\left.\frac{\partial u}{\partial x}\right|_{x=0}.
$$

Using a second-order central difference with a ghost node,

$$
\left.\frac{\partial u}{\partial n}\right|_{x=0}\approx\frac{u_{-1}^n-u_1^n}{2\Delta x}.
$$

Substituting this approximation into the general boundary condition gives

$$
a_0u_0^n+k_0\frac{u_{-1}^n-u_1^n}{2\Delta x}=f_0^n.
$$

Solving for the ghost-node temperature,

$$
u_{-1}^n=u_1^n+\frac{2\Delta x}{k_0}\left(f_0^n-a_0u_0^n\right).
$$

The FTCS discretization at the left boundary is

$$
u_0^{n+1}=
u_0^n+
r\left(
u_1^n-2u_0^n+u_{-1}^n
\right).
$$

Substituting the ghost-node expression yields

$$
u_0^{n+1}=
2ru_1^n+
\left(
1-2r-\frac{2r\Delta x a_0}{k_0}
\right)u_0^n+
\frac{2r\Delta x}{k_0}f_0^n.
$$

The same derivation applies to the right boundary using the outward-normal convention, giving

$$
u_N^{n+1}=
2ru_{N-1}^n+
\left(
1-2r-\frac{2r\Delta x a_N}{k_N}
\right)u_N^n+
\frac{2r\Delta x}{k_N}f_N^n.
$$

Both equations have the same algebraic structure and are implemented through a single function:

```python
def FTCS_boundary(u_boundary, u_adjacent, r, dx, a, k, f):

    return (
        2*r*u_adjacent
        + (1 - 2*r - 2*r*dx*(a/k))*u_boundary
        + 2*r*dx*(f/k)
    )
```

The left boundary is updated by

```python
u[n+1, 0] = FTCS_boundary(
    u[n, 0], u[n, 1],
    r, dx, a0, k0, f0[n]
)
```

while the right boundary uses

```python
u[n+1, -1] = FTCS_boundary(
    u[n, -1], u[n, -2],
    r, dx, aN, kN, fN[n]
)
```

For Neumann conditions, $a_b=0$, and the boundary-temperature contribution proportional to $a_b/k_b$ vanishes. For Robin conditions, both the boundary temperature and the forcing contribute to the update.

#### Boundary Forcing and Time Integration

An important distinction concerns the time level at which the boundary forcing is evaluated.

For Dirichlet conditions, the prescribed boundary temperature is imposed at the new time level, using $f_b^{n+1}$.

For Neumann and Robin conditions, FTCS uses the forcing at the current time level, $f_b^n$, because the explicit update is calculated entirely from quantities known at time $t_n$.

For example, suppose the boundary forcing is initially zero and changes to a constant value of one:

$$
f_b^0=0,\qquad f_b^n=1\quad(n\geq1).
$$

During the first FTCS update, from $n=0$ to $n=1$, a Neumann or Robin boundary uses $f_b^0=0$. The new forcing $f_b^1=1$ first contributes to the following update, from $n=1$ to $n=2$.

This distinction follows directly from the explicit time discretization and explains why the two boundary treatments use different indices in the implementation.

### Time Discretization of Boundary Conditions

The two numerical methods handle the boundary forcing differently.

For FTCS, Neumann and Robin boundary updates use the forcing at the current time level, $f_b^n$.

For Crank–Nicolson, the boundary forcing enters through the trapezoidal average of $f_b^n$ and $f_b^{n+1}$.

Consequently, sudden boundary changes can produce different first-step responses in the two methods.

Dirichlet conditions are imposed directly at the new time level in both methods.

### Mixed Boundary Conditions

The current implementation supports independent boundary-condition types at the left and right boundaries.

For example, the user can prescribe:

- Dirichlet on the left and Neumann on the right.
- Dirichlet on the left and Robin on the right.
- Neumann on the left and Robin on the right.

The program identifies the boundary-condition type from the coefficients:

| Coefficients | Boundary condition |
|---|---|
| $k_b=0,\ a_b\neq0$ | Dirichlet |
| $a_b=0,\ k_b\neq0$ | Neumann |
| $a_b\neq0,\ k_b\neq0$ | Robin |

The same numerical solver is used for all supported configurations.

## Graphical User Interface

The GUI is implemented in `Mini_1D_Heat_App.py` using Tkinter.

The application starts with a method-selection window, where the user chooses between FTCS and Crank–Nicolson.

After selecting a method, the parameter window opens.

The updated interface uses a **two-column layout**:

- **Left panel:** initial temperature distribution, boundary coefficients, initial and new boundary values, and numerical parameters.
- **Right panel:** mathematical boundary-condition formulation and explanations for Dirichlet, Neumann, and Robin conditions.

The user provides the following parameters:

| Parameter | Description |
|---|---|
| Initial Temperature | Mathematical expression defining $\phi(x)$ |
| Left and Right $a$ | Boundary temperature coefficients |
| Initial Left and Right Boundary Values | Boundary forcing values at $n=0$ |
| New Left and Right Boundary $f$ | Prescribed forcing expressions for $n\geq1$ |
| Left and Right $k$ | Boundary gradient coefficients |
| Spatial step `dx` | Spatial grid resolution |
| Time step `dt` | Temporal resolution |
| `L` | Domain length |
| `T` | Final simulation time |

Mathematical expressions such as `sin(pi*x)`, `cos(pi*x)`, and `sin(2*pi*t)` are evaluated using `numexpr`.

The GUI performs basic numerical input validation, while HeatLib performs additional checks on grid parameters, boundary coefficients, and FTCS stability.

After the simulation is completed, the temperature distribution is displayed through a Matplotlib animation.

The GUI remains separate from the numerical library, allowing the solvers to be used independently for numerical experiments and further analysis.

## Example Problem

Consider the initial temperature distribution

$$
u(x,0)=\cos(\pi x),
\qquad 0\leq x\leq1.
$$

Suppose Dirichlet boundary conditions are imposed, with both boundary temperatures changing to $0.5$ from $n=1$.

The initial boundary temperatures are

$$
u(0,0)=1,\qquad u(1,0)=-1.
$$

For $a_0=a_N=1$ and $k_0=k_N=0$, the corresponding forcing values are:

| Parameter | Value |
|---|---|
| Initial Temperature | `cos(pi*x)` |
| Left $a$ | `1` |
| Right $a$ | `1` |
| Initial Left Boundary Value | `1` |
| Initial Right Boundary Value | `-1` |
| New Left Boundary $f$ | `0.5` |
| New Right Boundary $f$ | `0.5` |
| Left $k$ | `0` |
| Right $k$ | `0` |

The initial cosine profile is preserved at $n=0$. From the first new time level onward, the boundary temperatures are prescribed as $0.5$.

The numerical solution then evolves according to the one-dimensional heat equation.

## Solution Visualization

The computed temperature field is stored in a two-dimensional NumPy array,

$$
u[n,i]\approx u(x_i,t_n).
$$

The numerical library returns the spatial coordinates, time coordinates, and temperature field:

```python
x, t, u = FTCS(...)
```

or

```python
x, t, u = CrankNicolson(...)
```

The solution is also visualized through a Matplotlib animation.

For small time steps, intermediate frames may be skipped to improve animation performance without changing the numerical solution.

Returning the numerical arrays allows additional post-processing, including comparison with analytical solutions, error analysis, and grid-convergence studies.

## Future Development

- **Numerical verification:** systematic comparison against analytical solutions and convergence studies for both numerical methods.
- **Variable thermal diffusivity:** extension to diffusivity coefficients that depend on space or temperature.
- **Additional visualization tools:** improved plotting controls, export of numerical results, and comparison of multiple simulations.
- **Extended physical models:** possible inclusion of source terms and advection–diffusion equations.

## References

Ferziger, J. H., & Perić, M. *Computational Methods for Fluid Dynamics*. Springer.

