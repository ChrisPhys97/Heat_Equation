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

Therefore, requiring all Gershgorin disks to lie within the unit disk $|\lambda|\leq1$ provides a sufficient condition for the linear stability of the numerical scheme. Smith, G. D. (1985).
<i>Numerical Solution of Partial Differential Equations: Finite Difference Methods</i>.
3rd ed. Oxford University Press.</p>

<p>
The amplification matrix $A$ does not need to be constructed explicitly to apply Gershgorin's Circle Theorem, since every discretized equation corresponds to one row of the matrix. Therefore, the diagonal and off-diagonal elements can be identified directly from the coefficients of each discretized equation. For example, considering the bottom boundary, excluding the corners, the explicit diffusion scheme is written as

$$
u_{0,i}^{n+1}=(1-3r)u_{0,i}^{n}+r(u_{0,i-1}^{n}+u_{0,i+1}^{n}+u_{1,i}^{n}).
$$

The diagonal element is $a_{ss}=1-3r$, corresponding to the temperature of the cell itself, while the three neighboring cells contribute to the off-diagonal elements. Consequently, the Gershgorin radius is

$$
R_s=|r|+|r|+|r|=3r,
$$

since $r\geq0$. Applying Gershgorin's Circle Theorem gives

$$
|\lambda-(1-3r)|\leq3r.
$$

Since the amplification matrix for the diffusion-only problem is symmetric, its eigenvalues are real, and the corresponding Gershgorin disk can be examined as an interval on the real axis,

$$
1-6r\leq\lambda\leq1.
$$

Requiring this interval to lie within $[-1,1]$ gives the sufficient stability restriction $r\leq1/3$ for the bottom-boundary rows. However, this condition alone does not guarantee the stability of the complete numerical scheme, since all the discretized equations must be considered.

The same procedure can be repeated for the interior cells, the four boundary sides, and the four corners. The corresponding diagonal elements, Gershgorin radii, and lower eigenvalue bounds are summarized below.

| Cell type | Diagonal element $a_{ss}$ | Radius $R_s$ | Lower bound |
|---|---|---|---|
| Interior | $1-4r$ | $4r$ | $1-8r$ |
| Bottom, top, right | $1-3r$ | $3r$ | $1-6r$ |
| Left | $1-5r$ | $3r$ | $1-8r$ |
| Bottom-left, top-left | $1-4r$ | $2r$ | $1-6r$ |
| Bottom-right, top-right | $1-2r$ | $2r$ | $1-4r$ |

From these results, the most restrictive lower bound is obtained from the interior cells and the left Dirichlet boundary, while the maximum upper bound is equal to one. Therefore, all eigenvalues of the amplification matrix are bounded by

$$
1-8r\leq\lambda\leq1.
$$

To satisfy the spectral stability requirement, $\rho(A)=\max_k|\lambda_k|\leq1$, the lower bound must be greater than or equal to $-1$. Thus,

$$
1-8r\geq-1
$$

and finally,

$$
\boxed{r\leq\frac14}.
$$

This provides a sufficient spectral stability condition for the complete two-dimensional diffusion scheme without explicitly constructing the amplification matrix or calculating its eigenvalues. For example, on a grid containing $100\times100$ computational cells, the amplification matrix would have dimensions $10{,}000\times10{,}000$, although only the nine types of discretized equations need to be examined to apply the theorem.

Although Gershgorin's Circle Theorem provides a sufficient condition for spectral stability, it does not necessarily guarantee that the numerical solution remains nonnegative. An additional restriction can be obtained by requiring all coefficients of the explicit discretization to remain nonnegative, ensuring positivity preservation for nonnegative initial conditions and source terms. This approach is also discussed in the [TU Delft MUDE Textbook — Diffusion Equation](https://mude.citg.tudelft.nl/book/2026/numerical_methods_for_PDEs/diffusion_equation.html).

For the present discretization, the most restrictive coefficient occurs at the left Dirichlet boundary, yielding $r\leq1/5$. When advection is included, the stability condition must account for both diffusion and transport. Requiring all stencil coefficients to remain nonnegative gives

$$
\boxed{5r+2c\leq1},
$$

where $r=a\Delta t/\Delta x^2$ and $c=U\Delta t/\Delta x$. This provides a sufficient condition for stability and positivity preservation.

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

