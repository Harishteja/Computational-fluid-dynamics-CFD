## 2D Convection–Diffusion

This problem involves combined **convection and diffusion** of heat in a 2D domain with a prescribed velocity field (see Fig. 4). The mesh data and velocity field are provided in the dataset folder.

**Domain properties:**

* Density: \$ρ = 1\$
* Thermal diffusivity: \$k/C\_p = 1/50\$
* Geometry ratios: \$h\_A/H = h\_C/H = 0.068\$

**Boundary conditions:**

* \$U\_A = 1\$, \$U\_B = 0\$, \$U\_C = 1\$, \$V\_D = 0\$
* \$T\_A = 20^\circ C\$ at boundary A
* At \$x = L\$ (except the outlet): \$T = 50^\circ C\$

**Governing equation:**

```math
\frac{\partial}{\partial x} (\rho U T) + \frac{\partial}{\partial y} (\rho V T) =
\frac{\partial}{\partial x} \left( \Gamma \frac{\partial T}{\partial x} \right) +
\frac{\partial}{\partial y} \left( \Gamma \frac{\partial T}{\partial y} \right) + S,
\quad \text{where } \Gamma = \frac{k}{C_p}.
```

### Discretization

* Diffusion terms → **Central Difference Scheme**
* Convection terms → **Hybrid Scheme**

  * Uses **Central Differencing** if \$|Pe| < 2\$
  * Uses **Upwind Differencing** if \$|Pe| > 2\$

Example for the west face:

```math
T_w = 
\begin{cases} 
\frac{T_W + T_P}{2}, & |Pe_w| < 2 \\ 
T_W, & Pe_w > 2 \\ 
T_P, & Pe_w < -2
\end{cases}
\quad \text{where } Pe_w = \frac{\rho U_w \delta x_w}{\Gamma}.
```

The discretized equations are solved using:

* **Gauss–Seidel Iteration**
* **2D TDMA (Tri-Diagonal Matrix Algorithm)** with directional sweeps

### Results

* The temperature field is obtained for the given domain (Fig. 5).
* Outlet temperature is computed from the solution.
* Additional studies include the effect of boundary conditions and iterative solvers on convergence.


