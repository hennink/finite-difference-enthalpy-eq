# Finite Difference Methods for the Non-linear Enthalpy Equation

This code accompanies the paper

    "A Pressure-based Solver for Low-Mach Number Flow using a Discontinuous Galerkin Method", A. Hennink, M. Tiberga, D. Lathouwers

that was submitted to the Journal of Computational Fluids on 2019-11.

The Fortran code is a finite-difference solver for the non-linear ODE

    d/dt (rho * h)  =  - lambda * h + Q ,

    h  =  cp * T - h0
       =  cp * (T - T0)

where 

* `lambda`,         is a constant,
* `h0`, `T0`, `cp`      are constants,
* `t`               is the independent variable (the time),
* `T`               is the temperature (the unknown),
* `h = h(t)`        is the specific enthalpy,
* `rho = rho(T)`    is the density, which is an arbitrary function of T,
* `Q = Q(t)`        is the source term.

It uses various techniques to deal with the non-linearity in the ODE. 

