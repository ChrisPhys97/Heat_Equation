<h2>Stochastic Heat Equation</h2>

<p>This 2D simulation originally started as a heat diffusion homework assignment during numerical methods courses and was progressively upgraded to its current form. Therefore, the project serves both educational and visualization purposes.</p>

<p>The upgrades include:</p>

<ul>
<li>Multiple localized heat sources,</li>

<li>Source locations whose <i>x</i>-coordinates follow the probability density function $$f(x)=\pi\sin(\pi x)$$ </li>

<li>while the <i>y</i>-coordinates follow a uniform distribution within $$[0.1, 0.9]$$ </li>

<li>Stochastic source activation modeled as a discrete-time approximation of a Poisson process,</li>

<li>Advection-diffusion transport,</li>

<li>Numba JIT and parallel acceleration for improved computational performance,</li>

<li>Time-dependent visualization and animation with Plotly.</li>

</ul>

</div>

<h3>Notebook and Code Description</h3>

<p>
<b>1.</b> Source coordinates generation using the Inverse Transform Sampling method, followed by a brief theoretical background. The probability density function was selected since this distribution is not directly available as a predefined sampling distribution in either <i>numpy.random</i> or <i>scipy.stats</i>. Additionally, the function is analytically integrable and invertible, allowing direct implementation of the inversion method. For more complex probability density functions without analytical inversion, the "hit-or-miss" method represents an alternative sampling approach.
</p>

<p>
<b>2.</b> Construction of localized Gaussian source fields and implementation of stochastic source activation through a discrete-time Poisson process.
</p>

<p>
<b>3.</b> Information regarding the discretization procedure, stability conditions, and boundary conditions. The finite volume formulation includes:
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

<li>Volume integrals approximated by the product of the cell-average value and the control-volume size, namely $$\overline{u}\Delta{V}\quad$$ where$$\overline{u}=u_{P}$$</li>

</ul>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?{\bar{u}}{\Delta}V,\qquad{\bar{u}=u_P}"/>
</p>

<ul>

<li>Time discretization by Forward Explicit Euler method,</li>

<li>Boundary Conditions: Linear interpolation.</li>

</ul>

<p>
The resulting formulation yields the same numerical scheme as the Finite Difference Method (FDM), namely the Forward Time-Centered Space (FTCS) scheme.
</p>

<p>
<b>4.</b> Time-dependent animation and visualization using Plotly.
</p>

<h3>Requirements</h3>

<ul>
<li>NumPy</li>
<li>Plotly</li>
<li>Numba</li>
</ul>

