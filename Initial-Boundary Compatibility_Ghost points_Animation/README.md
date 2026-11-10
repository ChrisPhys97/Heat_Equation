<h2>Initial-Boundary Compatibility, Ghost Points and Animation</h2>

**1. Initial-boundary compatibility:** 
<div style="text-align: justify;
text-justify: inter-word; font-size: 15px;">
Is refered to the agreement between the initial and boundary conditions at the corners of the space-time domain. At basic level, this requires their values to match at $(t=0)$. Higher-order compatibility conditions also involve their derivatives and are important for the regularity of the solution.</div> 

For example, consider a rod of length \(L=1)\, with the initial temperature profile: $$\quad u(x,0)=\sin(\pi x)$$

And suppose that the temperature at both ends remains zero: $$\quad u(0,t)=u(L,t)=0,\qquad t>0$$

In this case, the initial and boundary conditions are compatible because the sine profile gives zero temperature at both ends. In fact, this example satisfies the higher-order compatibility conditions as well, and the solution is smooth. 

An incompatible case arises when we retain the same initial temperature profile but impose different boundary conditions. For example, both ends of the rod are suddenly brought into contact with heat sources maintained at a temperature of 0.5. The boundary conditions can then be represented by a step function in time:

$$
u(0,t)=u(L,t)=
\begin{cases}
0, & t=0,\\
0.5, & t\gt 0.
\end{cases}
$$

<p>These conditions are incompatible with the initial profile because the temperature at the boundaries changes abruptly from 0 to 0.5. Consequently, the solution cannot remain continuous at the initial space-time corners. Nevertheless, a solution still exists for positive time, and its visualization shows the sudden "jump" in temperature at the boundaries. One possible approach is to introduce a time-dependent boundary temperature that increases gradually instead of changing instantaneously. This is explored in the second example of the notebook. In particular, we first introduce a linear ramp function:</p>

$$
u(0,t) = u(L,t) =
\begin{cases}
0.5\frac{t}{T_r} &  0 \lt t \leq T_r \\
0.5 &  t \gt T_r
\end{cases}
$$

where $T_{r}$ is the transition time. This function preserves continuity and satisfies the zeroth-order initial-boundary compatibility condition. However, it does not satisfy the first-order compatibility because its first derivative is discontinuous at $t=T_{r}$. Therefore, although the linear ramp removes the initial jump in temperature, it does not guarantee higher-order smoothness. Further treatments using smoother boundary functions are introduced and discussed in more detail in the notebook.<br>
</br>

**2. Ghost points:** When ghost points are employed to impose boundary conditions, two different approaches can be implemented. 
<ul>
<li>The first, used in the Stochastic Heat Equation project, is to express the temperature at the ghost points in terms of the neighboring interior values and the boundary conditions through interpolation. We then apply the discretized equations to the boundary cells, where some terms involve points outside the physical domain. These ghost-point terms are replaced using the expressions derived earlier, so the final equations contain only temperatures at points inside the domain.</li>
  
<br>

<li>The second approach is used in this project. Suppose that the physical domain contains $N$ cells. Instead of eliminating the ghost-point values from the equations, we calculate them using the boundary conditions and store them in additional positions of the solution array. This allows us to use the same discretized equation for all physical cells, including those next to the boundaries. When the equation is applied to a boundary cell, it involves indices outside the physical domain. However, these indices correspond to the ghost points, whose values have already been calculated and stored. Therefore, the same equation can be applied without any additional treatment of the boundary cells.</li>
</ul>

An illustration of the one-dimentional domain used in both notebook examples is shown in the next picture

                            Physical domain: x ∈ [0, L]

       Ghost     |              N physical cells             |   Ghost
        cell     |                                           |    cell
    ┌────────────┬────────────┬─────────────────┬────────────┬────────────┐
    │            │            │                 |            │            │
    │     ●      │     ●      │       ....      |     ●      │     ●      │
    │            │     1      │                 |     N      │            │
    └────────────┴────────────┴─────────────────┴────────────┴────────────┘
       Ghost_L      x = Δx/2                      x = L−Δx/2     Ghost_R
                 ↑                                           ↑
                x=0                                         x=L                                      

The **extended solution array**, which includes the two ghost-cell values, has $N+2$ entries:

$$
\mathbf{u}_{ext}=[u_{GL},u_1,u_2,\ldots,u_N,u_{GR}].
$$

However, the **physical numerical solution** contains only the $N$ temperatures calculated at the centers of the physical cells:

$$
\mathbf{u}_{phys}=[u_1,u_2,\ldots,u_N].
$$

The main advantage of this approach is that we can use the same discretized equation for all physical cells. The numerical loop runs over the indices $i=1,\ldots,N$, while the discretized equation also involves the neighboring values $u_{i-1}$ and $u_{i+1}$. When the loop reaches the first cell ($i=1$), the term $u_{i-1}$ refers to the left ghost cell ($u_{GL}$). Similarly, at the last cell ($i=N$), the term $u_{i+1}$ refers to the right ghost cell ($u_{GR}$).

Since these ghost-cell values have already been calculated from the boundary conditions and stored in the extended array, they are accessed directly through their indices. This allows us to apply the same equation throughout the physical domain without introducing separate discretized equations for the boundary cells.
<ul>
<li>The first and simplest example considers an insulated rod of length $L=1$, with the initial temperature profile

$$u(x,0)=\cos(\pi x)$$

Since both ends of the rod are insulated, the imposed boundary conditions are

$$\frac{\partial{u}}{\partial{x}}(0,t)=\frac{\partial{u}}{\partial{x}}(L,t)=0 \qquad (1)$$

This is an initially compatible problem, even at higher orders, since the spatial derivatives required by the boundary conditions agree with the initial temperature profile.

To impose these conditions, we use CDS for the first derivative at the boundary faces. At the right boundary, Eq. (1) gives

$$\frac{u^n_{GR}-u^n_N}{\Delta x}=0\quad\Rightarrow\quad u^n_{GR}=u^n_N \qquad (2)$$

Similarly, at the left boundary,

$$\frac{u^n_1-u^n_{GL}}{\Delta x}=0\quad\Rightarrow\quad u^n_{GL}=u^n_1 \qquad (3)$$

Since the algorithm uses two arrays that are swapped between time steps, the ghost-cell values must be updated before applying the discretized equation. In particular, Eqs. (2) and (3) are evaluated using the temperatures at time step $n$, ensuring that the correct boundary values are available when calculating the solution at time step $n+1$.</li>


<li>The second example considers a rod of length $L=1$, with the initial temperature profile

$$u(x,0)=\sin(\pi x)$$

Unlike the previous example, both ends of the rod are now brought into contact with heat sources, whose temperature gradually increases from 0 to 0.5. The imposed Dirichlet boundary conditions are

$$u(0,t)=u(L,t)=g(t), \qquad (1)$$

where $g(t)$ is initially defined as a linear ramp function:

$`g(t)=\begin{cases}\frac{0.5t}{T_r} & \text{if } 0\leq t<T_r \\ 0.5 & \text{if } t\geq T_r\end{cases}`$

Here, $T_r$ is the transition time required for the boundary temperature to reach 0.5.

This function satisfies the zeroth-order initial-boundary compatibility condition, since $g(0)=u(0,0)=u(L,0)=0$. However, it does not satisfy the first-order compatibility condition because its initial time derivative is nonzero, whereas the second spatial derivative of the initial temperature profile vanishes at both boundaries. Furthermore, the first derivative of the ramp function is discontinuous at $t=T_r$.

To impose the boundary conditions, we use linear interpolation between the ghost-cell centers and the neighboring physical-cell centers. Since the physical boundaries lie halfway between these centers, the temperature at the right boundary is approximated by

$$g(t_n)=\frac{u_N^n+u_{GR}^n}{2}\quad\Rightarrow\quad u_{GR}^n=2g(t_n)-u_N^n \qquad (2)$$

Similarly, at the left boundary,

$$g(t_n)=\frac{u_{GL}^n+u_1^n}{2}\quad\Rightarrow\quad u_{GL}^n=2g(t_n)-u_1^n \qquad (3)$$

As in the first example, the algorithm uses two arrays that are swapped between time steps. Therefore, the ghost-cell values are updated using Eqs. (2) and (3) before applying the discretized equation. Both expressions are evaluated at time step $n$, ensuring that the boundary conditions are correctly imposed when calculating the physical solution at time step $n+1$.

Although the linear ramp removes the initial discontinuity in temperature, it does not provide higher-order compatibility. For this reason, smoother boundary functions are also introduced in the notebook to investigate how the regularity of the boundary conditions affects the solution.</li>
</ul><br>
</br>

**3. Animation:**
The final aspect of this project concerns the visualization of the numerical solution. Since a cell-centered FVM is used, the temperatures are calculated at the centers of the physical cells rather than at the boundaries. Therefore, the first and last calculated values are located at $x=\Delta x/2$ and $x=L-\Delta x/2$, leaving a small gap between the numerical solution and the physical boundaries.

To visualize the temperature profile across the entire rod, we can introduce two additional plotting points at $x=0$ and $x=L$. However, their temperatures must first be determined from the imposed boundary conditions.

<li>For the general Neumann case, suppose that the temperature gradients at both boundaries are prescribed as</li>

$$
\frac{\partial u}{\partial x}(0,t)=q_L(t),
\qquad
\frac{\partial u}{\partial x}(L,t)=q_R(t).
$$

Using a central difference approximation at the boundary faces, the ghost-cell expressions become

$$
u_{GL}^n=u_1^n-q_L(t_n)\Delta x,
$$

$$
u_{GR}^n=u_N^n+q_R(t_n)\Delta x.
$$

Since the physical boundaries lie halfway between the ghost-cell centers and the neighboring physical-cell centers, their temperatures can be approximated through linear interpolation:

$$
u(0,t_n)\approx\frac{u_{GL}^n+u_1^n}{2}
=u_1^n-\frac{q_L(t_n)\Delta x}{2},
$$

$$
u(L,t_n)\approx\frac{u_N^n+u_{GR}^n}{2}
=u_N^n+\frac{q_R(t_n)\Delta x}{2}.
$$

For the first example, where the rod is insulated, $q_L=q_R=0$. Therefore, the boundary temperatures reduce to

$$
u(0,t_n)\approx u_1^n,
\qquad
u(L,t_n)\approx u_N^n.
$$

This means that the additional plotting points would simply repeat the temperatures of the first and last physical cells. For this reason, the animation in the first example is constructed directly from the physical cell-center values.

<li>For the second example, where time-dependent Dirichlet boundary conditions are imposed, the boundary temperatures are already prescribed by $g(t)$. Therefore, the plotting arrays are defined as</li>

$$
\mathbf{x}_{plot}=[0,x_1,x_2,\ldots,x_N,L],
$$

$$
\mathbf{u}_{plot}^{n}=[g(t_n),u_1^n,u_2^n,\ldots,u_N^n,g(t_n)].
$$

These additional boundary points are used only for visualization and do not affect the numerical solution or the ghost-cell treatment. The saved temperature profiles are then displayed at successive time steps using Matplotlib, allowing us to observe how the temperature evolves along the rod.
