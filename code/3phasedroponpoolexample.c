/* 3phasedroponpool.c

	Ben Fudge, 08/08/2022
	bfudge@live.co.uk

	3 Phase drop on pool impact

	Compile with 

		export OMP_NUM_THREADS=X

	where X is the number of threads you want to use

		qcc -O2 -w -fopenmp 3phasedroponpoolexample.c -o 3phasedroponpoolexample -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

	And executed with

		./3phasedroponpoolexample Re We Fr Real_Pool_Density Real_Pool_Viscosity Real_Drop_Air_ST Real_Pool_Air_ST Real_Drop_Pool_ST MaxLevel
	
	where
		-Real_Pool_Density is the pool density in kg/m^3
		-Real_Pool_Viscosity is the pool kinematic viscosity in cSt
		-Real_Drop_Air_ST is the surface tension between the droplet and air in mNm (i.e. water is 72) and the same for the other two pairs
		-MaxLevel is the maximum resolution level   

	e.g. 
	
		./3phasedroponpoolexample 5000 1500 20 1000 2 16 20 5 11

	The parameters have to be in that order, if an incorrect number of parameters are provided them the simulation will abort.
*/

// Sets harmonic mean for density which provides better performance for high viscosity ratios
#define mu(f1, f2, f3) (1./(clamp(f1,0,1)*1./mu1 + clamp(f2,0,1)*1./mu2 + clamp(f3,0,1)*1./mu3))

// Standard Basilisk files inclusion
#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "tag.h"

// Extra files required for three-phase and axisymmetric surface area calculation
#include "three-phaseCHIZARI.h"
#include "interface_area_axi.h"

#include "tension.h"

// Needed to sort out directories for results
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// Grid refinement
#define LEVEL 6
int MAXLEVEL;

// Non dimensional droplet size, velocity and density as well as initial distance from the surface
#define D0 1.0
#define V0 1.0
#define RHO1 1.0
#define S0 0.6

// Domain Size
#define SIZE 2.5

// Fixed Air Properties
#define REAL_AIR_DENSITY 1.2754 // Air Density in kg/m^3 (fixed)
#define REAL_AIR_VISC 18.1e-6 // Air Dynamic Viscosity in Pa.s (fixed)

// Fixed FC-770 Droplet Properties, would need to change if using a different droplet
#define REAL_DROP_DENSITY 1793 // Real Droplet Density in kg/m^3 (fixed)
#define REAL_DROP_VISCOSITY 1.4e-3 // Real Droplet Dynamic Viscosity in Pa.s (fixed)

// Varying Pool Properties to be provided on execution
double REAL_POOL_DENSITY; // Real Pool Density in kg/m^3 (given)
double REAL_POOL_VISCOSITY; // Real Pool Kinematic Viscosity in cSt (convenient for Silicone Oils) (given)

// Density Ratios
double RHO_AIR2DROP; // Air to Droplet Density Ratio (calculated)
double RHO_POOL2DROP; // Pool to Droplet Density Ratio (calculated)

// Viscosity Ratios
double MU_AIR2DROP; // Air to Droplet Dynamic Viscosity Ratio (calculated) 
double MU_POOL2DROP; // Pool to Droplet Dynamic Viscosity Ratio (calculated)

// Surface Tensions to be provided on execution
double REAL_DROP_AIR_ST; // Droplet-Air Surface Tension in mN/m (given)
double REAL_POOL_AIR_ST; // Pool-Air Surface Tension in mN/m (given)
double REAL_DROP_POOL_ST; // Droplet-Pool Surface Tension in mN/m (given)

// Dimensionless Groups to be provided on execution
double RE;
double WE;
double FR;

// For Surface Tensions
double surf_drop_air;
double sigmaSF_sigma_FA;
double surf_drop_pool;
double sigmaSA_sigma_FA;
double surf_pool_air;
double SIGMA1;
double SIGMA2;
double SIGMA3;

// For Adaptivity
double f1tol;
double f2tol;
double f3tol;
double ff1tol;
double ff2tol;
double ff3tol;
double vorttol;

// For Filtering
int filtersize;
double filtertol;

// For saving results
char resultsfolder[400];

char restring[100];
char westring[100];
char frstring[100];
char viscratiostring[100];
char maxlevelstring[100];

char interfacesfolder[400];
char snapshotsfolder[400];

char videoname[400];
char runstats[400];
char energylog[400];
char paramslog[400];
char maxlogs[400];


// Bottom of Pool = LEFT of domain
// No slip, no permeability
u.n[left]=dirichlet(0); // Impermeable
u.t[left]=dirichlet(0); // No Slip
p[left]=neumann(0);
pf[left]=neumann(0);

// Side of Pool = TOP of domain
// No slip and no permeability
u.n[top]=dirichlet(0); // Impermeable
u.t[top]=neumann(0); // Slip 
p[top]=neumann(0);
pf[top]=neumann(0);

// Above the Pool = RIGHT of domain
// Outflow
u.n[right]=neumann(0);
p[right]=dirichlet(0);
pf[right]=dirichlet(0);

// Default for right edge is symmetry

int main(int argc, char * argv[])
{
	// Makes sure that correct number of inputs are provided
	assert(argc==10);

	RE=atof(argv[1]);
	WE=atof(argv[2]);
	FR=atof(argv[3]);
	REAL_POOL_DENSITY=atof(argv[4]);
	REAL_POOL_VISCOSITY=atof(argv[5]);
	REAL_DROP_AIR_ST=atof(argv[6]);
	REAL_POOL_AIR_ST=atof(argv[7]);
	REAL_DROP_POOL_ST=atof(argv[8]);
	MAXLEVEL=atoi(argv[9]);

	// Calculate density and viscosity Ratios
	RHO_AIR2DROP=REAL_AIR_DENSITY/REAL_DROP_DENSITY;
	RHO_POOL2DROP=REAL_POOL_DENSITY/REAL_DROP_DENSITY;

	MU_AIR2DROP=REAL_AIR_VISC/REAL_DROP_VISCOSITY;

	double REAL_POOL_VISCOSITY_Pas=REAL_POOL_VISCOSITY*REAL_POOL_DENSITY/1e6; // Converts pool visc in cSt to Pas
	MU_POOL2DROP=REAL_POOL_VISCOSITY_Pas/REAL_DROP_VISCOSITY;

	// Creates folder and file names for saving results
	sprintf(resultsfolder, "./Results_3_Phase_Pool/");

	sprintf(restring, "Re_%1.03f_", RE);
	sprintf(westring, "We_%1.03f_", WE);
	sprintf(frstring, "Fr_%2.1f_", FR);
	sprintf(viscratiostring, "Visc_Ratio_%1.03f_", MU_POOL2DROP);
	sprintf(maxlevelstring,"MAX_LEVEL_%i", MAXLEVEL);

	strcat(resultsfolder, restring);
	strcat(resultsfolder, westring);
	strcat(resultsfolder, frstring);
	strcat(resultsfolder, viscratiostring);
	strcat(resultsfolder, maxlevelstring);

	strcpy(interfacesfolder, resultsfolder);
	strcpy(snapshotsfolder, resultsfolder);
	strcpy(videoname, resultsfolder);
	strcpy(runstats, resultsfolder);
	strcpy(energylog, resultsfolder);
	strcpy(paramslog, resultsfolder);
	strcpy(maxlogs, resultsfolder);

	strcat(interfacesfolder, "/Interfaces/");
	strcat(snapshotsfolder, "/Snapshots/");

	strcat(runstats, "/Simulation_Stats_");
	strcat(runstats, restring);
	strcat(runstats, westring);
	strcat(runstats, frstring);
	strcat(runstats, viscratiostring);
	strcat(runstats, maxlevelstring);

	strcat(videoname, "/Video_");
	strcat(videoname, restring);
	strcat(videoname, westring);
	strcat(videoname, frstring);
	strcat(videoname, viscratiostring);
	strcat(videoname, maxlevelstring);
	strcat(videoname, ".mp4");

	strcat(energylog, "/Energy_Log_");
	strcat(energylog, restring);
	strcat(energylog, westring);
	strcat(energylog, frstring);
	strcat(energylog, viscratiostring);
	strcat(energylog, maxlevelstring);
	
	strcat(maxlogs, "/Max_Values_Log_");
	strcat(maxlogs, restring);
	strcat(maxlogs, westring);
	strcat(maxlogs, frstring);
	strcat(maxlogs, viscratiostring);
	strcat(maxlogs, maxlevelstring);		
		
	strcat(paramslog, "/Simulation_Parameters");

	// Checks if folders already exist and if not automatically makes folders for results
	struct stat st1 = {0};

		if (stat("./Results_3_Phase_Pool", &st1) == -1) {
			mkdir("./Results_3_Phase_Pool", 0700);
		}

	struct stat st2 = {0};

		if (stat(resultsfolder, &st2) == -1) {
			mkdir(resultsfolder, 0700);
		}

	struct stat st3 = {0};

		if (stat(interfacesfolder, &st3) == -1) {
			mkdir(interfacesfolder, 0700);
		}

	struct stat st4 = {0};

		if (stat(snapshotsfolder, &st4) == -1) {
			mkdir(snapshotsfolder, 0700);
		}

	// Delete text files if they already exist
	if ( access( videoname, F_OK ) != -1 ) {
		remove(videoname);
	}

	if ( access( runstats, F_OK ) != -1 ) {
		remove(runstats);
	}

	if ( access( energylog, F_OK) != -1) {
		remove(energylog);
	}

	if ( access( maxlogs, F_OK) != -1) {
		remove(maxlogs);
	}

	// Simulation size and resolution
	size(SIZE*D0);
	origin(-0.5*SIZE,0);
	init_grid(1<<LEVEL);

	// Densities using earlier calculated ratios, phase 1 is the droplet, 2 the air and 3 the pool
	rho1 = RHO1;
	rho2 = RHO1 * RHO_AIR2DROP;
	rho3 = RHO1 * RHO_POOL2DROP;

	// Defines Viscosity to get correct Re for droplet and uses earlier calculated ratios
	// to get the viscosities for the other phases 
	mu1 = 1.0/(double)RE;
	mu2 = mu1 * MU_AIR2DROP;
	mu3 = mu1 * MU_POOL2DROP;

	// Defines surface tensions to get correct We for droplet and uses ratios
	// to get the tensions for the other phases 
	surf_drop_air = 1.0/WE;

	sigmaSF_sigma_FA = REAL_DROP_POOL_ST/REAL_DROP_AIR_ST;
	surf_drop_pool = surf_drop_air*sigmaSF_sigma_FA;

	sigmaSA_sigma_FA = REAL_POOL_AIR_ST/REAL_DROP_AIR_ST;
	surf_pool_air = surf_drop_air*sigmaSA_sigma_FA;

	// As there is a VOF for each phase it means that at each interface there are actually 
	// two VOFs present. Consequently we define the tension assigned to each phase such that
	// when two VOFs are present it results in the correct interfacial tension. For example
	// if phase 1 (the droplet) and phase 2 (the air) are present then SIGMA1+SIGMA2 sums to
	// surf_drop_air the droplet-air interfacial tension as required.
	SIGMA1 = 0.5*(surf_drop_air + surf_drop_pool - surf_pool_air);
	SIGMA2 = 0.5*(surf_drop_air + surf_pool_air - surf_drop_pool);
	SIGMA3 = 0.5*(surf_drop_pool + surf_pool_air - surf_drop_air);

	f1.sigma = SIGMA1;
	f2.sigma = SIGMA2;
	f3.sigma = SIGMA3;

	run();
}

double levelwidth; 

event init (t=0)
{
	// Width of bands for setting up the initial grid around the pool-air interface
	if (MAXLEVEL >= 10)
	{
		levelwidth=0.005;
	}
	else
	{
		levelwidth=0.015;
	}

	// Setting the pool to initially be in the very centre of the domain is a problem
	// as this will always be a cell boundary which is undesirable. Consequently we 
	// offset the pool interface to be 1/3 of a cell off the centreline to avoid this
	// problem. (1/3 is used as it shouldn't become a cell boundary if they change).
	double offset=SIZE*pow(2,-MAXLEVEL)/3;

	// Have initial grid highly refined around droplet and the pool-air interface.
	refine(sq(x-S0) + sq(y) < sq(0.505*D0) && sq(x-S0) + sq(y) > sq(0.495*D0) && level < MAXLEVEL);
	refine(fabs(x) < 32.0*levelwidth && level < MAXLEVEL-5);
	refine(fabs(x) < 16.0*levelwidth && level < MAXLEVEL-4);
	refine(fabs(x) < 8.0*levelwidth && level < MAXLEVEL-3);
	refine(fabs(x) < 6.0*levelwidth && level < MAXLEVEL-2);
	refine(fabs(x) < 4.0*levelwidth && level < MAXLEVEL-1);
	refine(fabs(x) < 1.0*levelwidth && level < MAXLEVEL);
	
	// Initialises droplet with centre S0 from pool level.
	fraction(f1, -sq(x-S0+offset)-sq(y)+sq(0.5*D0)); // Drop
	fraction(f2, difference(x+offset, -sq(x-S0+offset)-sq(y)+sq(0.5*D0))); // Air
	fraction(f3, difference(-x-offset, -sq(x-S0+offset)-sq(y)+sq(0.5*D0))); // Pool

	// Droplet has initial velocity towards pool
	foreach()
	{
		u.x[] = -V0*f1[];
	}

	// Save a snaphot of the initial setup
	char snapshotname[400];
	char snapshottimestamp[100];
	sprintf(snapshottimestamp, "Initial.gfs");

	strcpy(snapshotname, snapshotsfolder);
	strcat(snapshotname, snapshottimestamp);

	output_gfs(file = snapshotname, translate = true, t=t);

	// Saves the interfaces of the three phases at the initial setup
	char interfacename[400];

	char interfacetimestamp[100];
	sprintf(interfacetimestamp, "Interface_Initial");

	strcpy(interfacename, interfacesfolder);
	strcat(interfacename, interfacetimestamp);

	char interface2name[400];
	strcpy(interface2name, interfacename);
	strcat(interface2name,"_2");

	char interface3name[400];
	strcpy(interface3name, interfacename);
	strcat(interface3name,"_3");

	strcat(interfacename, "_1");

	FILE * fpp = fopen(interfacename, "w");
	output_facets(f1,fpp);
	fclose(fpp);

	FILE * fpp2 = fopen(interface2name, "w");
	output_facets(f2,fpp2);
	fclose(fpp2);

	FILE * fpp3 = fopen(interface3name, "w");
	output_facets(f3,fpp3);
	fclose(fpp3);

	// Defines the tolerances for the adaptivity
	f1tol=7.5e-3;
	f2tol=7.5e-3;
	f3tol=7.5e-3;
	ff1tol=7.5e-3;
	ff2tol=7.5e-3;
	ff3tol=7.5e-3;

	// Defines the parameters for the droplet (noise) removal
	filtersize=4;
	filtertol=1e-3;

	// Prints all simulation parameters to a file, useful if you need to check 
	// after the simulation is done.
	FILE * parlog = fopen(paramslog, "w");
	fprintf(parlog,"Level %i\nMax Level %i\nDomain Size %1.1f\nDroplet Dimensionless Diameter %1.1f\nDroplet Initial Dimensionless Velocity %1.1f\n",
	LEVEL,MAXLEVEL,SIZE,D0,V0); 
	fprintf(parlog,"Droplet Dimensionless Density %1.1f\nDroplet Dimensionless Distance Above Pool %1.1f\nAir to Droplet Density Ratio %1.6f\n",
	RHO1,S0,RHO_AIR2DROP);
	fprintf(parlog,"Pool to Droplet Density Ratio %1.6f\nAir to Droplet Viscosity Ratio %1.6f\nPool to Droplet Viscosity Ratio %1.6f\n",
	RHO_POOL2DROP,MU_AIR2DROP,MU_POOL2DROP);
	fprintf(parlog,"Real Droplet-Air Surface Tension (mN/m) %1.2f\nReal Droplet-Pool Surface Tension (mN/m) %1.2f\n",REAL_DROP_AIR_ST,REAL_DROP_POOL_ST);
	fprintf(parlog,"Real Pool-Air Surface Tension (mN/m) %1.2f\nRe %4.1f\nWe %3.1f\nFr %1.1f\nf1 (Droplet) Density %1.6f\nf2 (Air) Density %1.6f\n",
	REAL_POOL_AIR_ST,RE,WE,FR,rho1,rho2);
	fprintf(parlog,"f3 (Pool) Density %1.6f\nf1 (Droplet) Viscosity %1.6f\nf2 (Air) Viscosity %1.6f\nf3 (Pool) Viscosity %1.6f\n",rho3,mu1,mu2,mu3);
	fprintf(parlog,"Drop-Air Dimensionless Surface Tension %1.6f\nDrop-Pool/Drop Air Surface Tension Ratio %1.6f\n",surf_drop_air,sigmaSF_sigma_FA);
	fprintf(parlog,"Drop-Pool Dimensionless Surface Tension %1.6f\nAir-Pool/Drop Air Surface Tension Ratio %1.6f\n",surf_drop_pool,sigmaSA_sigma_FA);
	fprintf(parlog,"Air-Pool Dimensionless Surface Tension %1.6f\nf1 (Droplet) Surface Tension %1.6f\nf2 (Air) Surface Tension %1.6f\n",
	surf_pool_air,SIGMA1,SIGMA2);
	fprintf(parlog,"f3 (Pool) Surface Tension %1.6f\n",SIGMA3);
	fprintf(parlog,"Adaptivity Level Widths %1.4f\nf1 (Droplet) Adaptivity Threshold %1.4f\nf2 (Air) Adaptivity Threshold %1.4f\n",
	levelwidth,f1tol,f2tol);
	fprintf(parlog,"f3 (Pool) Adaptivity Threshold %1.4f\nf1 (Droplet) Boundary Adaptivity Threshold %1.4f\nf2 (Air) Boundary Adaptivity Threshold %1.4f\n",
	f3tol,ff1tol,ff2tol);
	fprintf(parlog,"f3 (Pool) Boundary Adaptivity Threshold %1.4f\nVorticity Adaptivity Threshold %1.4f\nFilter Size %i\nFilter Threshold %1.4f\n",
	ff3tol,vorttol,filtersize,filtertol);
	fclose(parlog);
}

event acceleration (i++)
{
	// Defines gravity to get correct Fr
	// Also in this case downwards is in the 
	// negative x direction
	face vector av=a;
	foreach_face(x)
	{
		av.x[] -= 1.0/(FR*FR);
	}
}


// Defines the fields for the kinetic, gravitiational potential, surface
// and viscous dissipated energies as well as any derivatives needed to 
// calculate them and their initial values to incriment them.
scalar energydissipationrate[];
scalar energydissipationratedroplet[];
scalar energydissipationrateair[];
scalar energydissipationratepool[];

scalar dudx[];
scalar dvdy[];
scalar dvdx[];
scalar dudy[];

scalar kineticenergypervol[];
scalar kineticenergypervoldroplet[];
scalar kineticenergypervolair[];
scalar kineticenergypervolpool[];

scalar gravenergypervol[];
scalar gravenergypervoldroplet[];
scalar gravenergypervolair[];
scalar gravenergypervolpool[];

scalar viewingfield[];
scalar mymu[];
scalar myrho[];

scalar vel2[];
scalar mylevel[];

scalar pdrop[];
scalar ppool[];
scalar viscstressdrop[];
scalar viscstresspool[];

double totaldissipatedenergy = 0;
double totaldissipateddropletenergy = 0;
double totaldissipatedpoolenergy = 0;
double totaldissipatedairenergy = 0;

double currentdissipatedenergy = 0;
double currentdissipateddropletenergy = 0;
double currentdissipatedpoolenergy = 0;
double currentdissipatedairenergy = 0;

double totalkineticenergy = 0;
double totalkineticdropletenergy = 0;
double totalkineticpoolenergy = 0;
double totalkineticairenergy = 0;

double totalgravenergy = 0;
double totalgravdropletenergy = 0;
double totalgravpoolenergy = 0;
double totalgravairenergy = 0;

double f1SurfaceArea;
double f2SurfaceArea;
double f3SurfaceArea;
double totalSurfaceArea;

double f1SurfaceEnergy;
double f2SurfaceEnergy;
double f3SurfaceEnergy;
double totalSurfaceEnergy;

double totalenergy;

double PressureDropMaxima;
double PressurePoolMaxima;
double ViscStressPoolMaxima;
double ViscStressDropMaxima;


event energiesandsurfaceareas (i++)
{
	foreach()
	{
		// Central difference derivates of both velocity components (u and v)
		// in both directions (x and y which correspond to r and z)
		dudx[]=((u.x[1,0]-u.x[-1,0])/(2.0*Delta));
		dvdy[]=((u.y[0,1]-u.y[0,-1])/(2.0*Delta));
		dvdx[]=((u.y[1,0]-u.y[-1,0])/(2.0*Delta));
		dudy[]=((u.x[0,1]-u.x[0,-1])/(2.0*Delta));

		// Useful fields for viewing things/checking snapshots
		viewingfield[]=0.000*f1[]+0.5*f2[]+1.000*f3[];
		mymu[]=mu(f1[], f2[], f3[]);
		myrho[]=rho(f1[], f2[], f3[]);

		vel2[]=u.x[]*u.x[]+u.y[]*u.y[];
		mylevel[]=level;

		// Pressures in the droplet and pool used to find their maxima
		pdrop[]=p[]*f1[];
		ppool[]=p[]*f3[];

		// Viscous stresses in the droplet and pool
		viscstressdrop[]=sqrt(f1[]*(dudx[]*dudx[]+dudy[]*dudy[]+dvdx[]*dvdx[]+dvdy[]*dvdy[]))*mymu[];
		viscstresspool[]=sqrt(f3[]*(dudx[]*dudx[]+dudy[]*dudy[]+dvdx[]*dvdx[]+dvdy[]*dvdy[]))*mymu[];

		// Viscous energy dissipation rate per unit volume, including for each phase
		energydissipationrate[]=2.0*mu(f1[], f2[], f3[])*(dudx[]*dudx[]+dvdy[]*dvdy[]+(u.y[]*u.y[]/(cm[]*cm[])))
        +mu(f1[], f2[], f3[])*(dvdx[]+dudy[])*(dvdx[]+dudy[]);

		energydissipationratedroplet[]=energydissipationrate[]*f1[];
		energydissipationrateair[]=energydissipationrate[]*f2[];
		energydissipationratepool[]=energydissipationrate[]*f3[];

		// Kinetic energy per unit volume, including for each phase
		kineticenergypervol[]=0.5*rho(f1[], f2[], f3[])*(u.x[]*u.x[] + u.y[]*u.y[]);

		kineticenergypervoldroplet[]=kineticenergypervol[]*f1[];
		kineticenergypervolair[]=kineticenergypervol[]*f2[];
		kineticenergypervolpool[]=kineticenergypervol[]*f3[];

		// Gravitiational potential energy per unit volume, including for each phase, takes the bottom
		// of the domain as the reference zero height
		gravenergypervol[]=rho(f1[], f2[], f3[])*(x+0.5*SIZE*D0)/(FR*FR);

		gravenergypervoldroplet[]=gravenergypervol[]*f1[];
		gravenergypervolair[]=gravenergypervol[]*f2[];
		gravenergypervolpool[]=gravenergypervol[]*f3[];	
	}

	// Find the maxima of the pressure and viscous stress in the droplet and pool
	PressureDropMaxima=statsf(pdrop).max;
	PressurePoolMaxima=statsf(ppool).max;
	ViscStressDropMaxima=statsf(viscstressdrop).max;
	ViscStressPoolMaxima=statsf(viscstresspool).max;

	// Finds the current energy dissipation rate by multiplying current value by cell
	// volume and summing over all cells. 2pi is needed for axisymmetry
    currentdissipatedenergy = statsf(energydissipationrate).sum*2.0*pi;
	currentdissipateddropletenergy = statsf(energydissipationratedroplet).sum*2.0*pi;
	currentdissipatedpoolenergy = statsf(energydissipationratepool).sum*2.0*pi;
	currentdissipatedairenergy = statsf(energydissipationrateair).sum*2.0*pi;

	// Incriments the total disspiated energy by multiplying the current rate by the 
	// timestep duration and adding to the currently accumulated value
    totaldissipatedenergy += currentdissipatedenergy*dt;
	totaldissipateddropletenergy += currentdissipateddropletenergy*dt;
	totaldissipatedpoolenergy += currentdissipatedpoolenergy*dt;
	totaldissipatedairenergy += currentdissipatedairenergy*dt;

	// Evaluates the current kinetic energy (again 2pi due to axisymmetry)
	totalkineticenergy = statsf(kineticenergypervol).sum*2.0*pi;
	totalkineticdropletenergy = statsf(kineticenergypervoldroplet).sum*2.0*pi;
	totalkineticpoolenergy = statsf(kineticenergypervolpool).sum*2.0*pi;
	totalkineticairenergy = statsf(kineticenergypervolair).sum*2.0*pi;

	// Evaluates the current graviatational potential energy (again 2pi due to axisymmetry)
	totalgravenergy = statsf(gravenergypervol).sum*2.0*pi;
	totalgravdropletenergy = statsf(gravenergypervoldroplet).sum*2.0*pi;
	totalgravpoolenergy = statsf(gravenergypervolpool).sum*2.0*pi;
	totalgravairenergy = statsf(gravenergypervolair).sum*2.0*pi;

	// Need to use a modified version of the interfacial area calculation for axisymmetric
	// setup, standard version only finds the arc lengths
	f1SurfaceArea=interface_area_axi(f1);
    f2SurfaceArea=interface_area_axi(f2);
    f3SurfaceArea=interface_area_axi(f3);
	totalSurfaceArea=f1SurfaceArea+f2SurfaceArea+f3SurfaceArea;

	// Calcutes the current surface energy of each phase from surface area multiplied
	// by the surface energies given by the tensions
	f1SurfaceEnergy=f1SurfaceArea*SIGMA1;
	f2SurfaceEnergy=f2SurfaceArea*SIGMA2;
	f3SurfaceEnergy=f3SurfaceArea*SIGMA3;
	totalSurfaceEnergy=f1SurfaceEnergy+f2SurfaceEnergy+f3SurfaceEnergy;

	// Calculates the (theoretically constant) total energy in the simulation
	totalenergy=totaldissipatedenergy+totalkineticenergy+totalgravenergy+totalSurfaceEnergy;
}

// Notes down the maximum pressure and viscous stress in both the droplet and pool in a log file
event maximanoting(t += 0.001)
{
	FILE * maxlog = fopen(maxlogs, "a");

	fprintf(maxlog, "%i %1.5f %1.5f %1.5f %1.5f %1.5f\n", i, t, PressureDropMaxima,PressurePoolMaxima, ViscStressDropMaxima, ViscStressPoolMaxima);
	fclose(maxlog);
}

// Notes down the various measured energies in a log file
event energylogging (t += 0.001)
{
 	FILE * enlog = fopen(energylog, "a");
	fprintf(enlog, "%i %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n"
	,i, t, currentdissipatedenergy, totaldissipatedenergy, currentdissipatedpoolenergy, totaldissipatedpoolenergy, 
	currentdissipateddropletenergy, totaldissipateddropletenergy, currentdissipatedairenergy, totaldissipatedairenergy, 
	totalkineticenergy, totalkineticdropletenergy, totalkineticpoolenergy, totalkineticairenergy,
	totalgravenergy, totalgravdropletenergy, totalgravpoolenergy, totalgravairenergy,
	totalSurfaceEnergy, f1SurfaceEnergy, f2SurfaceEnergy, f3SurfaceEnergy, totalenergy);
	fclose(enlog);
}

// Adaptation is based on the errors on the tracer fields and the velocity defined by earlier parameters
event controlledadapt (i++)
{
	scalar ff1[], ff2[], ff3[];
		
	foreach()
	{
		ff1[]=f1[];
		ff2[]=f2[];
		ff3[]=f3[];
	}
	boundary({ff1, ff2, ff3});

	adapt_wavelet( {f1, f2, f3, ff1, ff2, ff3, u.x, u.y }, (double[])
	{f1tol, f2tol, f3tol, ff1tol, ff2tol, ff3tol, 1e-2, 1e-2}, MAXLEVEL, LEVEL-1);
}


// Remove small droplets (Noise) defined by earlier parameters
event FilterDroplets (i++)
{
 		remove_droplets(f1, filtersize, filtertol);
 		remove_droplets(f1, filtersize, filtertol, true);
 		remove_droplets(f2, filtersize, filtertol);
 		remove_droplets(f2, filtersize, filtertol, true);
 		remove_droplets(f3, filtersize, filtertol);
 		remove_droplets(f3, filtersize, filtertol, true);
}

// Save video of simulation, has the interfaces and a differently coloured field for 
// each phase on the left hand side and the interfaces, mesh and mesh refinement level
// on the right hand side. Also displays the simulation time in the top left corner.
event viewing (t += 0.001)
{
	view(width=1900, height=1050, fov=22.5, ty= 0, quat = { 0, 0, -0.707, 0.707 });
	
    clear();
	draw_vof("f1", lw=2);
	draw_vof("f2", lw=2);
	draw_vof("f3", lw=2);
	squares("viewingfield", map = cool_warm, min = -0.2, max = 1.2);
	mirror({0,1}) {
		draw_vof("f1", lw=2);	
		draw_vof("f2", lw=2);
		draw_vof("f3", lw=2);
		cells(lw=0.25);
		squares("mylevel", map = cool_warm, min = 5, max = 11);
	} 

	char timestring[100];
	sprintf(timestring, "t=%2.03f",t);
	draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);

	save(videoname);
}


// Useful information about simulation
// Timestep, Current Sim Time, Current Timestep, Number of Gridpoints
// Wall Clock Time elapsed and CPU time elapsed
event SimStats (i += 10) 
{
	timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
	FILE * runlog = fopen (runstats, "a");
	fprintf (runlog, "%i %g %g %ld %g %g\n", i, t, dt, grid->n, s.real, s.cpu);
	fclose(runlog);
	// i, time, no of cells, real time elapsed, cpu time
}


// Saves Gerris view snapshots which are quite useful
event snapshot (t += 0.1)
{
	char snapshotname[200];
	char snapshottimestamp[100];
	sprintf(snapshottimestamp, "Snapshot-%1.02f.gfs", t);

	strcpy(snapshotname, snapshotsfolder);
	strcat(snapshotname, snapshottimestamp);

	output_gfs(file = snapshotname, translate = true, t=t);
}

// Saves the interfaces of the three phases 
event interfacelogging (t += 0.001)
{
	char interfacename[200];

	char interfacetimestamp[100];
	sprintf(interfacetimestamp, "Interface_%1.03f",t);

	strcpy(interfacename, interfacesfolder);
	strcat(interfacename, interfacetimestamp);

	char interface2name[200];
	strcpy(interface2name, interfacename);
	strcat(interface2name,"_2");

	char interface3name[200];
	strcpy(interface3name, interfacename);
	strcat(interface3name,"_3");

	strcat(interfacename, "_1");
	
	FILE * fpp = fopen(interfacename, "w");
	output_facets(f1,fpp);
	fclose(fpp);

	FILE * fpp2 = fopen(interface2name, "w");
	output_facets(f2,fpp2);
	fclose(fpp2);

	FILE * fpp3 = fopen(interface3name, "w");
	output_facets(f3,fpp3);
	fclose(fpp3);	
}


// Defines end time and writes final logfile entry
event end (t=0.25)
{
	// End Timing
	timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
	FILE * runlog = fopen (runstats, "a");
	fprintf (runlog, "%i %g %g %ld %g %g\n", i, t, dt, grid->n, s.real, s.cpu);
	fclose(runlog);
	// i, time, no of cells, real time elapsed, cpu time
}
