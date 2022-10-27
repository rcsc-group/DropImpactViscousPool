# Drop impact onto viscous pools

Code to simulation the impact of a droplet onto a deep pool of another immiscible fluid using the [Basilisk](http://basilisk.fr/) flow solver.

## Installation

-Uses a standard [installation](http://basilisk.fr/src/INSTALL) of Basilisk with the extra visualisation [packages](http://basilisk.fr/src/gl/INSTALL) required to produce video outputs.

-The three-phase method of [Chizari](http://basilisk.fr/sandbox/chizari/threephase/) is included here to allow independent variation of the droplet, pool and surrounding air fluids including interfacial tensions.

-[gfsview](http://gfs.sourceforge.net/wiki/index.php/Main_Page) can also be installed in order to view the outputed simulation slices.

## Setup

Currently the code is setup to simulate the impact of a [Fluorinert FC-770](https://www.3m.co.uk/3M/en_GB/p/d/b40006507/) droplet onto a pool of variable properties defined by the user surrounded by air at standard conditions. The droplet properties (density and viscosity) can be changed by editing the values on lines 71 and 72 of the main file

```
// Fixed FC-770 Droplet Properties, would need to change if using a different droplet
#define REAL_DROP_DENSITY 1793 // Real Droplet Density in kg/m^3 (fixed)
#define REAL_DROP_VISCOSITY 1.4e-3 // Real Droplet Dynamic Viscosity in Pa.s (fixed)
```
