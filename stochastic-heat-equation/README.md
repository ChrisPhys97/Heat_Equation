## Stochastic Heat Equation

<div style="text-align: justify;text-justify: inter-word; font-size: 18px;">
  
This 2D simulation originally started as a heat diffusion homework assignment during numerical methods courses and was progressively upgraded to its current form. Therefore, the project serves both educational and visualization purposes.

The upgrades include:
- Multiple localized heat sources,
- Source locations whose x-coordinates follow the probability density function $$f(x)=\pi \sin(\pi x)$$ while the y-coordinates follow a uniform distribution within \[0.1,0.9]\.
- Stochastic source activation modeled as a discrete-time approximation of a Poisson process,
- Advection-diffusion transport,
- Numba JIT and parallel acceleration for improved computational performance,
- Time-dependent visualization and animation with Plotly.

</div>

### Notebook and code description

<div style="text-align: justify;text-justify: inter-word; font-size: 18px;">

1. Source coordinates generation using the Inverse Transform Sampling method, followed by a brief theoretical background. The probability density function was selected since this distribution is not directly available as a predefined sampling distribution in either <i>numpy.random</i> or <i>scipy.stats</i>. Additionally, the function is analytically integrable and invertible, allowing direct implementation of the inversion method. For more complex probability density functions without analytical inversion, the "hit-or-miss" method represents an alternative sampling approach.

2. Construction of localized Gaussian source fields and implementation of stochastic source activation through a discrete-time Poisson process.

3. Information regarding the discretization procedure, stability conditions, and boundary conditions. The finite volume formulation includes:

- Surface integrals evaluated using the midpoint rule,
- Interpolation schemes:
  - Upwind Differencing Scheme (UDS) for the advective terms,
  - Central Differencing Scheme (CDS) for the diffusive terms,
- Volume integrals approximated by the product of the cell-average value and the control-volume size, namely

$$
\bar{u}\,\Delta V,
\qquad \text{with} \qquad \bar{u}=u_{P}.
$$
- Time discretization by Forward Explict Euler method.
- Boundary Conditions: Linear interpolation.
  
The resulting formulation yields the same numerical scheme as the Finite Difference Method (FDM), namely the Forward Time-Centered Space (FTCS) scheme.

4. Time-dependent animation and visualization using Plotly.
   
</div>
   
