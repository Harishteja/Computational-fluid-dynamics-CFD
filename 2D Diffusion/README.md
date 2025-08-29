This work demonstrates the implementation of the Finite Volume Method (FVM) in MATLAB for solving a 2D diffusion problem. All MATLAB scripts are organized within the MATLAB codes folder, along with the corresponding reports. The sections below provide an overview of the governing equations, the discretization strategies employed, and a discussion of the obtained results.

---
### Domain and boundary conditions

The Governing diffusion equation for temperature T (i.e., heat conduction equation), which can be written as,

$$ 0 = \frac{\partial}{\partial x}\left( k \frac{\partial T}{\partial x} \right)+ \frac{\partial}{\partial y}\left( k \frac{\partial T}{\partial y} \right) + S $$
<p align="center">
  <img align=center width="560" height="359" alt="image" src="https://github.com/user-attachments/assets/6d3751af-5cbe-47c4-acf7-aedf993bfddb" />
</p>
<p align="center"><em>Figure 1: Computational domain</em></p>
The problem is defined with the following boundary conditions:

- **Boundary 1:** $T = 15^\circ \text{C}$
- **Boundary 2:** $T(y) = 5\left(1 - \frac{y}{H}\right) + 15 \sin\left(\frac{\pi y}{H}\right)$
- **Boundary 3:** $T = 10^\circ \text{C}$
- **Boundary 4 (Insulated):** $\dfrac{\partial T}{\partial x} = 0$

The thermal conductivity varies spatially as:
$k = 16 \left(\frac{y}{H} + 1\right)$

and a uniform heat source term is specified as:
$S = -1.5 \quad \text{per unit area}.$

---
### Discretization stencil

In this study, the **Finite Volume Method (FVM)** is employed for discretization, using a central difference scheme along with source term linearization. The governing differential equation is first integrated over each control volume. 

The computational domain is divided into a set of control volumes, where each cell is represented by a central control volume **P**, surrounded by its neighboring cells to the east (**E**), west (**W**), north (**N**), and south (**S**). The faces of the control volume are denoted as east (**e**), west (**w**), north (**n**), and south (**s**). 

The grid spacing is defined by: $\Delta x$, $\Delta y$, while the half-width distances from the cell center to its respective faces are represented as: $\delta x_e$ \; $\delta x_w$ \; $\delta y_n$ \; $\delta y_s$.
This discretization stencil provides the framework for evaluating fluxes across the control volume boundaries and forms the foundation of the finite volume formulation.

---

### Discretization
The governing differential equation is integrated over a 2D control volume :

$$
0 = \int_w^e \int_s^n \frac{\partial}{\partial x}\left( k \frac{\partial T}{\partial x} \right) \ dx \ dy + \int_s^n \int_w^e \frac{\partial}{\partial y}\left( k \frac{\partial T}
{\partial y} \right) \ dx \ dy + \int_w^e \int_s^n S \ dx \ dy
$$

Here, the source term \( S \) is assumed uniform over the control volume, so its integral reduces to \( S $\Delta x \Delta y \$).

$$
\big[ \ k_e \ \Delta y \left(\frac{\partial T}{\partial x}\right)_e - k_w \ \Delta y \left(\frac{\partial T}{\partial x}\right)_w \ \big] + \big[ \ k_n \ \Delta x \left(\frac{\partial T}{\partial y}\right)_n - k_s \ \Delta x \left(\frac{\partial T}{\partial y}\right)_s \ \big] + S \ \Delta x \ \Delta y = 0
$$

At each cell center \( P \), gradients at the faces (east, west, north, south) are evaluated using the **central difference scheme**:


$$
\left(\frac{\partial T}{\partial x}\right)_e = \frac{T_E - T_P}{\delta x_e}, \quad
\left(\frac{\partial T}{\partial x}\right)_w = \frac{T_P - T_W}{\delta x_w}, \quad
\left(\frac{\partial T}{\partial y}\right)_n = \frac{T_N - T_P}{\delta y_n}, \quad
\left(\frac{\partial T}{\partial y}\right)_s = \frac{T_P - T_S}{\delta y_s}
$$

Substituting these into the integrated equation gives a balance in the standard finite volume form:

$$
a_P T_P = a_E T_E + a_W T_W + a_N T_N + a_S T_S + S \Delta x \Delta y
$$

---

### Coefficients
The coefficients are expressed as:

$$
a_W = \frac{k_w \Delta y}{\delta x_w}; \quad
a_E = \frac{k_e \Delta y}{\delta x_e}; \quad
a_S = \frac{k_s \Delta x}{\delta y_s}; \quad
a_N = \frac{k_n \Delta x}{\delta y_n}
$$

$$
a_P = a_W + a_E + a_S + a_N
$$

























