<h2>Stochastic Heat Equation</h2>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

This 2D simulation originally started as a heat diffusion homework assignment during numerical methods courses and was progressively upgraded to its current form. Therefore, the project serves both educational and visualization purposes.

<p>The upgrades include:</p>

<ul>
<li>Multiple localized heat sources,</li>

<li>
Source locations whose <i>x</i>-coordinates follow the probability density function
</li>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?f(x)=\pi\sin(\pi&space;x)" />
</p>

<li>
while the <i>y</i>-coordinates follow a uniform distribution within
<img src="https://latex.codecogs.com/svg.image?[0.1,0.9]" />.
</li>

<li>Stochastic source activation modeled as a discrete-time approximation of a Poisson process,</li>

<li>Advection-diffusion transport,</li>

<li>Numba JIT and parallel acceleration for improved computational performance,</li>

<li>Time-dependent visualization and animation with Plotly.</li>

</ul>

</div>

<h3>Notebook and Code Description</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

<p>
1. Source coordinates generation using the Inverse Transform Sampling method, followed by a brief theoretical background. The probability density function was selected since this distribution is not directly available as a predefined sampling distribution in either <i>numpy.random</i> or <i>scipy.stats</i>. Additionally, the function is analytically integrable and invertible, allowing direct implementation of the inversion method. For more complex probability density functions without analytical inversion, the "hit-or-miss" method represents an alternative sampling approach.
</p>

<p>
2. Construction of localized Gaussian source fields and implementation of stochastic source activation through a discrete-time Poisson process.
</p>

<p>
3. Information regarding the discretization procedure, stability conditions, and boundary conditions. The finite volume formulation includes:
</p>

<ul>

<li>Surface integrals evaluated using the midpoint rule,</li>

<li>
Interpolation schemes:
<ul>
<li>Upwind Differencing Scheme (UDS) for the advective terms,</li>
<li>Central Differencing Scheme (CDS) for the diffusive terms,</li>
</ul>
</li>

<li>
Volume integrals approximated by the product of the cell-average value and the control-volume size, namely
</li>

</ul>

<p align="center">
<img src="https://latex.codecogs.com/png.image?\bar{u}{\Delta} V"/>
</p>

<ul>

<li>Time discretization by Forward Explicit Euler method,</li>

<li>Boundary Conditions: Linear interpolation.</li>

</ul>

<p>
The resulting formulation yields the same numerical scheme as the Finite Difference Method (FDM), namely the Forward Time-Centered Space (FTCS) scheme.
</p>

<p>
4. Time-dependent animation and visualization using Plotly.
</p>

</div>

<h3>Requirements</h3>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

<ul>
<li>NumPy</li>
<li>Plotly</li>
<li>Numba</li>
</ul>

</div>
