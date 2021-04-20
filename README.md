# 2-D Compressible Euler Equations -- Hydrodynamic Shock-Tube Problem (C++ Implementation)
## Author: Alexandros Papados ##

---------------------------------------------------------------------------------------------------------------------------------
This repository is dedicated to provide users of interests with the ability to solve 2-D hydrodynamic shock-tube problems using 
the finite volume method with flux limiting in C++. The 2-D problem of interest stems from (Liska, 2003). The paper may be found below.

This code makes use of RK-2 to integration in time, the HLLE flux solver, and the monotonized central flux limiter. As of now, the initial state is 
hardcoded, but may be easily changed to solve various shock-tube problems. The solution to this problem is shown in the figure below.

<p align="center"><img src=2D_Euler.png width="650" height="550" /></p>
                             
<p align="center">Initial State of density and pressure (top) Final State of density and pressure (bottom)</p>

---------------------------------------------------------------------------------------------------------------------------------
## In-progress ##

* Implement various flux solvers and limiters
* Generalized initial condition function
* MPI implementation
* Write-up discussing schemes and problem

---------------------------------------------------------------------------------------------------------------------------------
## References ##
Paper: https://www.researchgate.net/publication/228781823_Comparison_of_Several_Difference_Schemes_on_1D_and_2D_Test_Problems_for_the_Euler_Equations

Code: https://github.com/wme7/Euler


