# Land-spill example: utah_dem

This example simulates the overland flow with one point source, and the point
source continues leaking fluid into the computational domain. The working
fluid is Maya crude.

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

After simulation, execute `$ python plotframes.py` to plot the flow at each 
output time frame. The generated plots will be in the folder `_plots`.

If GIS software will be used to visualize the result, users can use the script
`writenc.py` to create a multidimension NetCDF file. The output NetCDF file has
been tested with [QGIS](https://www.qgis.org/) and [ArcGIS Pro](https://pro.arcgis.com/).

Instead of a simple `$ make all`, users can do a fine control of the process 
with the following steps:

#### 1. Download topography:

Either
```
$ make topo
```
or
```
$ python ./download_topo.py
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
