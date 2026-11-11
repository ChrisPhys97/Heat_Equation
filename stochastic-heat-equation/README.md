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

<li>Volume integrals approximated by the product of the cell-average value and the control-volume size, namely </li>

</ul>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?{\bar{u}}{\Delta}V,\qquad{\bar{u}=u_P}"/>
</p>

<ul>

<li>Time discretization by Forward Explicit Euler method,</li>

<li>Boundary Conditions: Linear interpolation.</li>

</ul>

<p>
The resulting formulation yields the same numerical scheme as the Finite Difference Method (FDM) for the interior cells, namely the Forward Time-Centered Space (FTCS) scheme for diffusion, combined with an upwind discretization for advection. By eliminating the ghost points algebraically through the imposed boundary conditions, four separate equations are obtained for the boundary sides and another four for the corresponding corners. This approach highlights the difficulty of deriving the correct boundary equations, as well as the need to implement additional equations beyond the general interior stencil. The procedure becomes even more challenging when dealing with three-dimensional problems. Although explicitly storing ghost points may provide a more convenient alternative for such numerical simulations, deriving the boundary equations can simplify the stability analysis by making the stencil coefficients directly available. A brief derivation of a sufficient stability condition is presented below.</p>

<p>The stability condition depends on what you want to guarantee. Requiring all stencil coefficients to be nonnegative gives a sufficient condition for monotonicity, while Gershgorin's Circle Theorem can provide a less restrictive sufficient condition for spectral stability.</p> 

<p>Following the first approach, the stability condition can be derived more directly. The explicit numerical scheme can be expressed in matrix form as

$$
\mathbf{u}^{n+1}=A\mathbf{u}^{n},
$$

where $\mathbf{u}^{n}$ is the vector of temperature values at time level $n$, and $A$ is the amplification matrix containing the coefficients of the discretized equations, including the boundary conditions. For the numerical scheme to be linearly stable, the spectral radius of $A$ must satisfy

$$
\rho(A)=\max_k|\lambda_k|\leq1,
$$

where $\lambda_k$ denotes the eigenvalues of $A$. Rather than calculating these eigenvalues explicitly, Gershgorin's Circle Theorem can be used to estimate their possible locations. According to the theorem, every eigenvalue of $A$ lies within at least one disk defined by

$$
|\lambda-a_{ss}|\leq R_s,
$$

where $a_{ss}$ is the diagonal element of row $s$ and $R_s$ is the sum of the absolute values of its off-diagonal elements,

$$
R_s=\sum_{j\neq s}|a_{sj}|.
$$

Therefore, requiring all Gershgorin disks to lie within the unit disk $|\lambda|\leq1$ provides a sufficient condition for the linear stability of the numerical scheme.</p>




<p>
<b>4.</b> Time-dependent animation and visualization using Plotly.
</p>

<h3>Requirements</h3>

<ul>
<li>NumPy</li>
<li>Plotly</li>
<li>Numba</li>
</ul>

