# 1D-Finite-Volume-Solver-for-Elliptic-PDEs
This repository contains a MATLAB implementation of a 1D finite difference (cell-centered finite volume) method for solving linear elliptic partial differential equations. It supports user-defined grids, boundary conditions, and exact solutions, making it versatile for numerical experimentation.

# Description 
The script ELLIPTIC1d.m solves the elliptic equation: -(ku_x)_x = f on (a,b) with options for Dirichlet and Neumann boundary conditions. The default domain is (a,b) = (0,1), and the user can modify the conductivity k(x) and the source term f(x) in the script. 

# Key Features 
- Custom Grid: Supports uniform and non-uniform grids.
- Boundary Conditions: Dirichlet or Neumann conditions specified via input arguments.
- Exact Solution Comparison: Computes error norms if the exact solution is provided.
- Visualization: Plots numerical and exact solutions (if available).
- Matrix Assembly: Implements transmissibility calculations for diffusion terms.

# Usage 
Syntax: 
```matlab
[xc, nsol] = ELLIPTIC1d(nxdx, bcond, ifexact, ifplot, ifdemo)
```

Inputs: 
- nxdx: Number of grid cells for a uniform grid (integer), Non-uniform grid node positions (vector).
- bcond: Vector [flag_a, val_a, flag_b, val_b] specifying boundary conditions. flag_a and flag_b: 0 for Dirichlet, non-zero for Neumann. val_a and val_b: Corresponding boundary values.
- ifexact: 0 if No exact solution provided. 1 if Exact solution known for comparison.
- ifplot: 1 to plot solutions (default), 0 for no plot.
- ifdemo: 1 for verbose matrix info, 0 for standard operation.

Outputs:
- xc: Cell center locations.
- nsol: Numerical solution at cell centers.

Example: 
Non-uniform Grid, No Exact Solution: 
```matlab
ELLIPTIC1d([0, 1/3, 1/2, 2/3, 0.9, 1], [0, 0, 0, 1], 0, 1, 0);
```
Uniform Grid, Exact Solution:
```matlab
ELLIPTIC1d(5, [0, 0, 0, 0], 1, 1, 0);
```
Mixed Boundary Conditions 
```matlab
ELLIPTIC1d(5, [0, 1, 0, 0], 1, 1, 0);
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.


