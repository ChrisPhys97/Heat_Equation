<h2>Initial-Boundary Compatibility, Ghost Points and Animation</h2>

<div style="text-align: justify;
text-justify: inter-word; font-size: 15px;">
As already mentioned, this project aims to demonstrate three aspects of numerical methods through two simple one-dimensional examples of the heat equation. The first aspect, "initial-boundary compatibility", refers to the agreement between the initial and boundary conditions at the corners of the space-time domain. At the most basic level, this requires their values to match at \(t=0)\. Higher-order compatibility conditions also involve their derivatives and are important for the regularity of the solution. It is worth mentioning that initial-boundary compatibility differs from the consistency of a numerical scheme, which concerns the local truncation error, and from the well-posedness of a problem, which concerns the existence, uniqueness, and continuous dependence of its solution.
</div>

<div style="text-align: justify;
text-justify: inter-word; font-size: 15px;">
For example, consider a rod of length \(L=1)\, with the initial temperature profile

$$
u(x,0)=\sin(\pi x).
$$

Suppose that the temperature at both ends remains zero:

$$
u(0,t)=u(L,t)=0,\qquad t>0.
$$

</div> 

In this case, the initial and boundary conditions are compatible because the sine profile gives zero temperature at both ends. In fact, this example satisfies the higher-order compatibility conditions as well, and the solution is smooth.

On the other hand, an incompatible case arises when we retain the same initial temperature profile but impose different boundary conditions. For example, suppose that both ends of the rod are suddenly brought into contact with heat sources maintained at a temperature of 0.5. The boundary conditions can then be represented by a step function in time:
</div>

$$
u(0,t)=u(L,t)=
\begin{cases}
0, & t=0,\\
0.5, & t>0.
\end{cases}
$$

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

These conditions are incompatible with the initial profile because the temperature at the boundaries changes abruptly from 0 to 0.5. Consequently, the solution cannot remain continuous at the initial space-time corners. Nevertheless, a solution still exists for positive time, and its visualization shows the sudden change in temperature at the boundaries.

One possible approach is to introduce a time-dependent boundary temperature that increases gradually instead of changing instantaneously. This is explored in the second example of the notebook. In particular, we first introduce a linear ramp function:
</div>

```math
u(0,t) = u(L,t) =
\begin{cases}
0.5\frac{t}{T_r} & \text{if } 0 \le t < T_r \\
0.5 & \text{if } t \ge T_r
\end{cases}
```

where $T_{r}$ is the transition time.

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">
  
This function preserves continuity and satisfies the zeroth-order initial-boundary compatibility condition. However, it does not satisfy the first-order compatibility condition for the chosen initial profile. Furthermore, its first derivative is discontinuous at $t=T_{r}$.

Therefore, although the linear ramp removes the initial jump in temperature, it does not guarantee higher-order smoothness. Further treatments using smoother boundary functions are introduced and discussed in more detail in the notebook.
</div>

