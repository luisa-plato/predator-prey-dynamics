# predator-prey-dynamics

This is a finite element code to simulate the following predator-prey dynamics over a two-dimensional domain $\Omega$ and a finite time horizon $T>0$:
```math
\begin{align*}
\partial_t u - \nu \Delta u + \kappa \nabla \cdot (u \nabla w) 
		 &= (\alpha w - \beta) u  &&\text{ in } (0,T) \times \Omega,\label{eq:model_u}\\
		\partial_t w - \mu \Delta w
		&= (\gamma - \delta u) w &&\text{ in } (0,T) \times \Omega,\\
		\nabla u \cdot \boldsymbol{n} &= 0 &&\text{ on } [0,T] \times \Gamma,\\
		\nabla w \cdot \boldsymbol{n} &= 0 &&\text{ on }  [0,T] \times \Gamma .
\end{align*}
```


In a conda environment, where FEniCS is installed you can run the code via:

`python pred_pey.py 1 0.1 0.01` 

where the first argument is the prey-taxis coefficient $\kappa$, the second argument is the diffusion coefficient of the predator $\nu$
and the third argument is the diffusion coefficient of the prey $\mu$.

All other parameters of the system and the simulation need to be changed within the script.
