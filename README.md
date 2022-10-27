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

1. Reynolds number of the impacting droplet based on the droplet diameter. Defined as ` $\textrm{Re} = \rho_dDV_0/\mu_d$ `
