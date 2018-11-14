# Land-spill example: inclined_plane

This folder contains an example of overland flow with one point source that
continues leaking fluid into the computational domain.

The domain is an inclined plane with an angle of 2.5 degree. The working fluid
is silicon oil (kinematic viscosity 1.13e-3 m^2 / sec). And the volumetric rate
of the point source is 1.48e-6 m^3 / sec.

## Build and run the example

First, set up the environment variables necessary for using Clawpack. 
For example

```shell
$ export FC=gfortran 
$ export CLAW=$HOME/Sync/repos/clawpack-git
$ export PYTHONPATH=$HOME/Sync/repos/clawpack-git:$PYTHONPATH
```

Next, a simple `$ make all` will create the topography file, compile the 
program, setup input parameters, and then run the simulation.

After the simulation, open the file `plot.claw` with 
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
to visualize the results.

Alternatively, the Python script `plotresults.py` can be used to plot the flow
at each output time step. The generated plots will be in the folder `_plots`.

In the folder `lister_1992`, there are `.csv` files for the experimental data
obtained from Lister's paper (1992). If one uses `plotresults.py` to plot the 
flow, the experimental data will also be plotted.

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
$ FFLAGS=<compiler flags> make .exe
```
to compile and build the executable with customized flags and additional 
libraries.
For example, to enable OpenMP and NetCDF supports, use
```
$ FFLAGS="-O2 -fopenmp -DNETCDF -lnetcdff -I${NETCDF_INCLUDE_DIR}" make .exe
```
Both compiler flags and linker flags go to `FFLAGS`.

(Depends on your system, the flag for Fortran binding of NetCDF may be `-lnetcdf`
instead of `-lnetcdff`.)

#### 3. Set up input parameters

The input parameters for the simulation is in `setrun.py`. 
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

Lister, J. R. (1992). Viscous flows down an inclined plane from point and line sources. Journal of Fluid Mechanics, 242, 631-653.
