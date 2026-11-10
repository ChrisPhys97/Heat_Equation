<h2>Initial-Boundary Compatibility, Ghost Points and Animation</h2>

**1. Initial-boundary compatibility:** 
Is refered to the agreement between the initial and boundary conditions at the corners of the space-time domain. At basic level, this requires their values to match at \(t=0)\. Higher-order compatibility conditions also involve their derivatives and are important for the regularity of the solution. 

For example, consider a rod of length \(L=1)\, with the initial temperature profile $$u(x,0)=\sin(\pi x)$$

And suppose that the temperature at both ends remains zero: $$u(0,t)=u(L,t)=0,\qquad t>0$$

In this case, the initial and boundary conditions are compatible because the sine profile gives zero temperature at both ends. In fact, this example satisfies the higher-order compatibility conditions as well, and the solution is smooth. On the other hand, an incompatible case arises when we retain the same initial temperature profile but impose different boundary conditions. For example, both ends of the rod are suddenly brought into contact with heat sources maintained at a temperature of 0.5. The boundary conditions can then be represented by a step function in time:

$$
u(0,t)=u(L,t)=
\begin{cases}
0, & t=0,\\
0.5, & t\gt 0.
\end{cases}
$$

These conditions are incompatible with the initial profile because the temperature at the boundaries changes abruptly from 0 to 0.5. Consequently, the solution cannot remain continuous at the initial space-time corners. Nevertheless, a solution still exists for positive time, and its visualization shows the sudden "jump" in temperature at the boundaries. One possible approach is to introduce a time-dependent boundary temperature that increases gradually instead of changing instantaneously. This is explored in the second example of the notebook. In particular, we first introduce a linear ramp function:

$$
u(0,t) = u(L,t) =
\begin{cases}
0.5\frac{t}{T_r} &  0 \lt t \leq T_r \\
0.5 &  t \gt T_r
\end{cases}
$$

where $T_{r}$ is the transition time. This function preserves continuity and satisfies the zeroth-order initial-boundary compatibility condition. However, it does not satisfy the first-order compatibility because its first derivative is discontinuous at $t=T_{r}$. Therefore, although the linear ramp removes the initial jump in temperature, it does not guarantee higher-order smoothness. Further treatments using smoother boundary functions are introduced and discussed in more detail in the notebook.

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

<li>The first and simplest example intoduces an insulated rod of lenght $L=1$ and initial temperature $$u(x,0)=\cos(\pi x)$$ 

The imposed boundary conditions are $$\frac{\partial{u(0,t)}}{\partial{x}} = \frac{\partial{u(L,t)}}{\partial{x}} = 0 \qquad (1)$$ 

An already compatible case in all orders because all functions and derivatives of both initial and boundary conditions are in agreement. To employ the boundary conditions, a CDS interpolation is implimented to Eq. (1). Thus,
$$\frac{u^n_{GR}-u^n_N}{\DELTA x}=0 \Rightarrow u^n_{GR}=u^n_N \qquad(2)$$
Similary at the left boundary,
$$\frac{u^n_1-u^n_{GL}}{\DELTA x}=0 \Rightarrow u^n_1=u^n_{GL} \qquad(3)$$
Becasue swapping arrays are used, expressions (2,3) is updated before the descritized equation since the interpolation implimented at time step $n$.</li>
