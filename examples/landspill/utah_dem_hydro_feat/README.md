# Land-spill example: utah_dem_hydro_feat

This folder contains an example of overland flow with one point source that
continues leaking fluid into the computational domain. The difference between
this case and the `utah_dem` example is that now we have hydrologic features
in the domain. Whenever working fluid flows into a hydrologic feature, we will
assume the hydrologic feature captures all working fluid. From the computational
perspective, the solver removes working fluid if a cell is marked as hydrologic
feature cell.

Currently, only ESRI ASCII raster file format is supported. Hydrologic features
are usually in vector formats. It's users' responsibility to convert the vector
files to raster files before simulations.

## Build and run the example

First, set up the environment variables necessary for using Clawpack. 
For example

```shell
$ export FC=gfortran 
$ export CLAW=$HOME/Sync/repos/clawpack-git
$ export PYTHONPATH=$HOME/Sync/repos/clawpack-git:$PYTHONPATH
```

Next, a simple `$ make all` will download the topography file, compile the 
program, setup input parameters, and then run the simulation.

After the simulation, open the file `plot.claw` with 
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
to visualize the results.

Alternatively, the Python script `plotframes.py` can be used to plot the flow
at each output time step. The generated plots will be in the folder `_plots`.

If GIS software will be used to visualize the result, users can use the script
`writenc.py` to create a multidimension NetCDF file. The output NetCDF file has
been tested with [QGIS](https://www.qgis.org/) 
and [ArcGIS Pro](https://pro.arcgis.com/).

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

#### 2. Download hydrologic features

```
$ make hydro_feature
```
or
```
$ python ./download_hydro_feature.py
```

#### 3. Build the executable of the solver

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

#### 4. Set up input parameters

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

#### 5. Run the simulation

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
