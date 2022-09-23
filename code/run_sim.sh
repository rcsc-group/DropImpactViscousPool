#!/bin/bash

export OMP_NUM_THREADS=6

qcc -O2 -w -fopenmp 3phasedroponpoolexample.c -o 3phasedroponpoolexample -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

#	RE=atof(argv[1]);
#	WE=atof(argv[2]);
#	FR=atof(argv[3]);
#	REAL_POOL_DENSITY=atof(argv[4]);
#	REAL_POOL_VISCOSITY=atof(argv[5]);
#	REAL_DROP_AIR_ST=atof(argv[6]);
#	REAL_POOL_AIR_ST=atof(argv[7]);
#	REAL_DROP_POOL_ST=atof(argv[8]);
#	MAXLEVEL=atoi(argv[9]);


# 2cSt K=165,000
./3phasedroponpoolexample 6660 2020 25.9 934 100 14 20 5 10
