# Finite Difference Methods for the Non-linear Enthalpy Equation

This code accompanies the paper

```citation
"A Pressure-based Solver for Low-Mach Number Flow
using a Discontinuous Galerkin Method",
A. Hennink, M. Tiberga, D. Lathouwers
```

that was submitted to the Journal of Computational Fluids on 2019-11.

The Fortran code is a finite-difference solver for the non-linear ODE

```math
d/dt (rho · h)  =  - lambda · h + Q ,

h  =  cp · T - h0
   =  cp · (T - T0)
```

where

* `lambda`,         is a constant;

* `h0`, `T0`, `cp`      are constants;

* `t`               is the independent variable (the time);

* `T`               is the temperature (the unknown);

* `h = h(t)`        is the specific enthalpy;

* `rho = rho(T)`    is the density, which is an arbitrary function of T;

* `Q = Q(t)`        is the source term.

It uses various techniques to deal with the non-linearity in the ODE.

## Usage

All code is in a single file, so compilation is straightforward.
Simply `gfortran fd.f90 -o fd`,
or `ifort fd.f90 -o fd`,
or the equivalent.

Invoke the program with command line options
to set the fluid properties and the solution strategy.
For example,

```bash
    ./fd \
        fluid_props=affine lambda=0.1 cp=1.0 rho0=0.5 rho1=2.0 \
        T0=1.0 \
        derivhr_strategy=deriv_rh order_BDF=2 order_EX=3 nsteps=2048
```

(The rational behind putting all options in the command line is that
it is easier to loop over values from a script,
compared to if an input file were used.)

The output is of the form

```text
    <relative-error> <uniq-rho-h> <nonpos-ddh-rhoh> <nonpos-implicit-weight> <h-too-low> <h-too-high>
```

where `<relative-error>` is `(T_numerical - T_exact) / T_exact` at the end of the finite difference time-stepping,
and the remaining terms indicate whether a particular condition was encountered during the time-stepping,
with `0`/`1` indicating `false`/`true`:

* `<nonuniq-rho-h>`:
there is no one-to-one correspondence between
the temperature *T* and the volumetric enthalpy *(rho · h*);

* `<nonpos-ddh-rhoh>`:
for some value of the specific enthalpy *h*,
the derivative **d(*rho · h*) / d(*h*)** is negative;

* `<nonpos-implicit-weight>`:
There was a time step during which
the implicit coefficient of the finite difference scheme became non-positive;

* `<h-too-low>`:
There was some time step during which
***h*\* <= -*cp*\* / *beta*\***,
where the \* indicates a predictor value,
and *beta* is the thermal expansibility;

* `<h-too-high>`:
There was some time step during which
***h*\* >= +*cp*\* / *beta*\***.
