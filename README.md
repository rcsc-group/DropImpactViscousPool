# Drop impact onto viscous pools

Code to simulation the impact of a droplet onto a deep pool of another immiscible fluid using the [Basilisk](http://basilisk.fr/) flow solver.

## Installation

-Uses a standard [installation](http://basilisk.fr/src/INSTALL) of Basilisk with the extra visualisation [packages](http://basilisk.fr/src/gl/INSTALL) required to produce video outputs.

-The three-phase method of [Chizari](http://basilisk.fr/sandbox/chizari/threephase/) is included here to allow independent variation of the droplet, pool and surrounding air fluids including interfacial tensions.

-[gfsview](http://gfs.sourceforge.net/wiki/index.php/Main_Page) can also be installed in order to view the outputed simulation slices.

## Setup

Currently the code is setup to simulate the impact of a [Fluorinert FC-770](https://www.3m.co.uk/3M/en_GB/p/d/b40006507/) droplet onto a pool of variable properties defined by the user surrounded by air at standard conditions. The droplet properties (density and viscosity) can be changed by editing the values on lines 71 and 72 of the main file. (The air properties could also be changed by editing lines 67 and 68, for example to represent a change in air pressure).

```
// Fixed FC-770 Droplet Properties, would need to change if using a different droplet
#define REAL_DROP_DENSITY 1793 // Real Droplet Density in kg/m^3 (fixed)
#define REAL_DROP_VISCOSITY 1.4e-3 // Real Droplet Dynamic Viscosity in Pa.s (fixed)
```

The rest of the simulations parameters are provided as inputs to the main function when the code is instantiated. These are (and have to be in this order

1. Reynolds number of the impacting droplet based on the droplet diameter. Defined as $\textrm{Re} = \rho_dDV_0/\mu_d$.
2. Weber number of the impacting droplet based on the droplet diameter. Defined as $\textrm{We} = \rho_dDV_0^2/\sigma_{da}$.
3. Froude number of the impacting droplet based on the droplet diameter. Defined as $\textrm{Fr} = V_0/\sqrt{gD}$.
4. The real pool density in $\textrm{kgm}^{-3}$ (i.e. for water this would be 1000).
5. The real pool kinematic viscosity in $\textrm{cSt}$ (centistokes) as this is usually how silicone oils (used in the work here) are described. The conversion to dynamic viscosity to find the viscosity ratio to the droplet is calculated automatically.
6. The surface tension coefficient between the droplet and the air in $\textrm{mNm}$ (i.e for water this would be 72).
7. The surface tension coefficient between the pool and the air in $\textrm{mNm}$.
8. The interfacial tension coefficient between the droplet and the pool in in $\textrm{mNm}$.
9. The maximum resolution level of the simulations (i.e. the maximum resolution corresponds to $2^{\textrm{MAXLEVEL}}$ grid points per domain size).

## Execution

If running on multiple cores the number of cores needs to first be defined by 
```
export OMP_NUM_THREADS=6
```
Then the code needs to be compiled (including links to the graphical libraries if output videos are desired
```
qcc -O2 -w -fopenmp 3phasedroponpoolexample.c -o 3phasedroponpoolexample -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm
```
and then instantiated with the parameters defined above 
```
./3phasedroponpoolexample Re We Fr Real_Pool_Density Real_Pool_Viscosity Real_Drop_Air_ST Real_Pool_Air_ST Real_Drop_Pool_ST MaxLevel
```
for example
```
./3phasedroponpoolexample 6660 2020 25.9 934 100 14 20 5 10
```
An example shell script which automates this process is included.
