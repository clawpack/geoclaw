# Land-spill example: inclined_plane

This example simulates the overland flow with one point source, which continues 
leaking fluid into the computational domain.

The domain is an inclined plane with an angle of 2.5 degree. The working fluid
is silicone oil (kinematic viscosity 1.13e-3 m^2 / sec). And the volumetric rate
of the point source is 1.48e-6 m^3 / sec.

Note, the grid resolution of this simulation is much higher than other exmaples, 
so it may take a long time to finish the run on a personal computer or a laptop. 
On a personal computer with Intel i7-5930K (6 physical cores @ 3.5GHz), it takes 
about 30 hours to finish.

## Build and run the example

First, set up the necessary environment variables. For example

```shell
$ export FC=gfortran 
$ export CLAW=$HOME/clawpack
$ export PYTHONPATH=$HOME/clawpack:$PYTHONPATH
```

Next, a simple `$ make all` will create the topography file, compile the 
program, setup input parameters, and then run the simulation. Environment
variable `FFLAGS` can be used for optimization, and `OMP_NUM_THREADS` is used 
to specify number of threads used if OpenMP is enabled. For example, to enable
OpenMP and machine-specific optimization, and to use 6 threads: 

```
$ FFLAGS="-O2 -fopenmp -march=native" OMP_NUM_THREADS=6 make all
```

After simulation, execute `$ python plotresults.py` to plot the flow at each 
output time frame. The generated plots will be in the folder `_plots`.

In the folder `lister_1992`, there are `.csv` files for the experimental data
obtained from Lister's paper (1992). The experimental data are plotted in the
plots created by `plotresults.py`.

Instead of a simple `$ make all`, users can do a fine control of the process 
with the following steps:

#### 1. Create topography:

Either
```
$ make topo
```
or
```
$ python ./produce_topo.py
```

The python script is written and tested with Python 3.6. Other versions are not
tested.

#### 2. Build the executable of the solver

Use
```
$ OMP_NUM_THREADS=<# of threads> FFLAGS=<compiler flags> make .exe
```
to compile and build the executable with customized flags and additional 
libraries.
For example, to enable OpenMP and NetCDF supports, use
```
$ FFLAGS="-O2 -fopenmp -DNETCDF -lnetcdff -I${NETCDF_INCLUDE_DIR}" make .exe
```
Both compiler flags and linker flags go to `FFLAGS`.

(Depending on the system, the flag for Fortran binding of NetCDF may be `-lnetcdf`
instead of `-lnetcdff`.)

#### 3. Set up input parameters

The input parameters for the simulation are in `setrun.py`. 
After setting up or modifying `setrun.py`, use
```
$ make .data
```
or 
```
$ python ./setrun.py geoclaw
```
to create ASCII files that the solver understands.

If `make` complains that `.data` is up to date, then just remove the hidden file
`.data`.

#### 4. Run the simulation

Run
```
$ make .output
```
to launch the simulation.
If `make` complains that `.output` is up to date, then just remove the hidden file
`.output`.

`$ make .output` uses a python script provided by Clawpack to launch the executable
and set some CMD arguments for users. If launching the simulations with simple
`$ ./xgeoclaw`, the simulation will run, but output files will not be located 
in proper locations.

## Reference

Lister, J. R. (1992). Viscous flows down an inclined plane from point and line 
sources. Journal of Fluid Mechanics, 242, 631-653.
