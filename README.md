# Problem Statement
The objective is to solve the two-dimensional steady-state heat conduction equation for a rectangular plate with no internal heat generation.
The plate is discretized using a Uniform Grid and Non Uniform Grid, and appropriate boundary conditions are applied. 
The left and bottom sides of the plate are subjected to Neumann boundary conditions (specified heat flux), while the remaining boundaries are treated with Dirichlet boundary conditions.


The aim is to implement Node Based and Cell Centred approaches correctly using the finite difference method and obtain the steady-state temperature distribution within the 2D domain.
In short the problem must be solved for the following four combinations:
1. Node-based formulation with a uniform grid
2. Node-based formulation with a non-uniform grid
3. Cell-centered formulation with a uniform grid
4. Cell-centered formulation with a non-uniform grid


# Governing Equation
The steady-state heat equation is defined as $\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$.

## Boundary Conditions

### Dirichlet (Fixed Temperature)
Temperature is directly specified on boundaries (top and right):$T = T_{\text{specified}}$.

### Neumann (Zero Heat Flux)
Temperature gradient is specified (left and bottom):$\frac{\partial T}{\partial x} = 0$.

Finite difference form:$\frac{T_1 - T_0}{\Delta x} = 0 \;\Rightarrow\; T_0 = T_1$




