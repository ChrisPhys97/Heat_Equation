<h2>Initial-Boundary Compatibility, Ghost points and Animation</h2>

<div style="text-align: justify; text-justify: inter-word; font-size: 18px;">

As already mentioned this project aims to presents three aspects of numerical methods through just two simple one-dimentional examples of Heat Equation. The first aspect of initial-boundary compatibility is refered as the value matching of the intial and boundary condioions at $t=0$, which is not always the case. It is worth to mention that initial-boundary compatibility differs from the usual compatibity where numerical solution should tend to the analytical one and also differs for the aspect of well-posted problem. For example, a compatible case is when there is a rod of length $L=1$ having initial temperature porfile, 

$$u(x,0)=\sin(\pi x)$$ 

and the boundary conditions at both sides are 

$$u(0,t)=u(L,t)=0$$ for $t>0$. 

In this case $u$ at space-time domain is continous because sine at both boundaries give zero temperature. 

On the other hand, an incompatible case is when we want the same initial profile but we want to impose different boundary condition, for example, we want the rod to touch a heat source in the edges where temperature remains at 0.5, which means $u(0,t)=u(L,t)=0.5$ for $t>0$. This is translated as a step function of time and $u$ becomes discontinous in the space-time of boundaries,

$$
u(0 & L,t)=
\begin{cases}
0, & t=0,\\[4pt]
0.5, & t>0.
\end{cases}
$$
</div>
