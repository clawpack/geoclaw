#!/usr/bin/env python

"""

Classes representing parameters for GeoClaw runs

:Classes:

 - GeoClawData
 - RefinementData
 - TopographyData
 - FGoutData
 - FGmaxData
 - DTopoData
 - QinitData
 - SurgeData
 - MultilayerData
 - FrictionData
 - BoussData
 - GridData1D
 - BoussData1D

:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees
 - LAT2METER factor to convert degrees in latitude to meters
"""

import os
import numpy as np
import warnings

import clawpack.clawutil.data

# Radius of earth in meters.
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi
LAT2METER = Rearth * DEG2RAD

class GeoClawData(clawpack.clawutil.data.ClawData):
    r"""
    Object containing the basic .

    Note that this data object will write out multiple files.
    """
    def __init__(self):
        super(GeoClawData,self).__init__()

        # GeoClaw physics parameters
        self.add_attribute('gravity',9.8)
        self.add_attribute('rho', 1025.0)  # Density of water kg/m^3
        self.add_attribute('rho_air',1.15) # Density of air kg/m^3
        self.add_attribute('ambient_pressure', 101.3e3) # Nominal atmos pressure
        self.add_attribute('earth_radius',Rearth)
        self.add_attribute('coordinate_system',1)
        self.add_attribute('sphere_source',1)  # New starting in v5.10.0
        self.add_attribute('coriolis_forcing',True)
        self.add_attribute('theta_0',45.0)
        self.add_attribute('friction_forcing',True)
        self.add_attribute('manning_coefficient',[0.025])
        self.add_attribute('manning_break',[])

        # GeoClaw algorithm parameters
        self.add_attribute('dry_tolerance',1e-3)
        self.add_attribute('friction_depth',1.0e6)
        self.add_attribute('sea_level',0.0)
        self.add_attribute('speed_limit',50.)


    def write(self,data_source='setrun.py', out_file='geoclaw.data'):

        self.open_data_file(out_file, data_source)

        self.data_write('gravity',
                               description="(gravitational acceleration m/s^2)")
        self.data_write('rho', description="(Density of water kg/m^3)")
        self.data_write('rho_air',description="(Density of air kg/m^3)")
        self.data_write('ambient_pressure',
                                description="(Nominal atmospheric pressure Pa)")
        self.data_write('earth_radius', description="(Radius of the earth m)")
        self.data_write('coordinate_system',
                        description="(1=meters, 2=lon-lat)")
        self.data_write('sphere_source',
                        description="(0=none, 1=only in mass eqn, 2=all)")
        self.data_write('sea_level')
        self.data_write()

        # Forcing terms
        self.data_write('coriolis_forcing')
        if self.coordinate_system == 1 and self.coriolis_forcing:
            self.data_write('theta_0')
        self.data_write('friction_forcing')
        if self.friction_forcing:
            if type(self.manning_coefficient) in [int,float]:
                self.manning_coefficient = [self.manning_coefficient]
            num_manning = len(self.manning_coefficient)
            if len(self.manning_break) != num_manning - 1:
                raise IOError("***manning_break array has wrong length")
            self.data_write(value=num_manning,alt_name='num_manning')
            self.data_write('manning_coefficient')
            self.data_write('manning_break')
            self.data_write('friction_depth')

        self.data_write()

        self.data_write('dry_tolerance')
        self.data_write('speed_limit')

        self.close_data_file()



class RefinementData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(RefinementData,self).__init__()

        # Refinement controls
        self.add_attribute('wave_tolerance',1.0e-1)
        self.add_attribute('speed_tolerance',[1.0e12]*6)
        self.add_attribute('deep_depth',None)      # deprecated
        self.add_attribute('max_level_deep',None)  # deprecated
        self.add_attribute('variable_dt_refinement_ratios',False)


    def write(self,data_source='setrun.py', out_file='refinement.data'):
        # Refinement controls
        self.open_data_file(out_file, data_source)
        self.data_write('wave_tolerance')

        # check if user set deprecated parameters:
        if self.deep_depth is not None:
            w = '\n  *** WARNING: deep_depth parameter ignored as of v5.8.0'
            warnings.warn(w, UserWarning)
        if self.max_level_deep is not None:
            w = '\n  *** WARNING: max_level_deep parameter ignored as of v5.8.0'
            warnings.warn(w, UserWarning)

        if not isinstance(self.speed_tolerance,list):
            self.speed_tolerance = [self.speed_tolerance]
        self.data_write('speed_tolerance')
        self.data_write()
        self.data_write('variable_dt_refinement_ratios',
                        description="(Set dt refinement ratios automatically)")
        self.close_data_file()



class TopographyData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(TopographyData,self).__init__()

        # Topography data
        self.add_attribute('topo_missing', 99999.0)
        self.add_attribute('test_topography', 0)
        self.add_attribute('override_order', False)
        self.add_attribute('topofiles', [])

        # Jump discontinuity
        self.add_attribute('topo_location',-50e3)
        self.add_attribute('topo_left',-4000.0)
        self.add_attribute('topo_right',-200.0)
        self.add_attribute('topo_angle',0.0)

        # Simple oceanic shelf
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.008)


    def write(self,data_source='setrun.py', out_file='topo.data'):

        self.open_data_file(out_file, data_source)
        self.data_write(name='topo_missing',
                        description='replace no_data_value in topofile')
        self.data_write(name='test_topography',description='(Type topography specification)')
        if self.test_topography == 0:
            ntopofiles = len(self.topofiles)
            self.data_write(value=ntopofiles,alt_name='ntopofiles')
            self.data_write(name="override_order", description="(Override order topo files are used)")
            for tfile in self.topofiles:

                if len(tfile) == 6:
                    w = '\n  *** WARNING: topofile specs changed in v5.8.0 -- ' + \
                          'Flag level info now ignored'
                    warnings.warn(w, UserWarning)
                    tfile = [tfile[0], tfile[-1]] # drop minlevel,maxlevel,t1,t2
                elif len(tfile) == 2:
                    pass  # now expect only topo_type, filename
                else:
                    raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
                self._out_file.write("\n'%s' \n " % fname)
                self._out_file.write("%3i   # topo_type\n" % tfile[0])
        elif self.test_topography == 1:
            self.data_write(name='topo_location',description='(Bathymetry jump location)')
            self.data_write(name='topo_left',description='(Depth to left of bathy_location)')
            self.data_write(name='topo_right',description='(Depth to right of bathy_location)')
        elif self.test_topography == 2 or self.test_topography == 3:
            self.data_write(name='x0',description='(Location of basin end)')
            self.data_write(name='x1',description='(Location of shelf slope end)')
            self.data_write(name='x2',description='(Location of beach slope)')
            self.data_write(name='basin_depth',description='(Depth of basin)')
            self.data_write(name='shelf_depth',description='(Depth of shelf)')
            self.data_write(name='beach_slope',description='(Slope of beach)')
        else:
            raise NotImplementedError("Test topography type %s has not been"
                                        " implemented." % self.test_topography)

        self.close_data_file()


class FixedGridData(clawpack.clawutil.data.ClawData):

    """
    Deprecated, starting in 5.9.0 use FGoutData instead.
    """

    def __init__(self):

        super(FixedGridData,self).__init__()

        # Fixed Grids
        self.add_attribute('fixedgrids',[])


    def write(self,data_source='setrun.py', out_file='fixed_grids.data'):
        # Fixed grid settings
        msg = 'rundata.fixed_grid_data is deprecated starting in v5.9.0,' \
            + ' use rundata.fgout_data instead'
        #warnings.warn(msg)
        if len(self.fixedgrids) > 0:
            raise AttributeError(msg)


class FGoutData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(FGoutData,self).__init__()

        # File name for fgout points and parameters:
        self.add_attribute('fgout_grids',[])


    def write(self,data_source='setrun.py', out_file='fgout_grids.data'):
        self.open_data_file(out_file, data_source)
        num_fgout_grids = len(self.fgout_grids)
        self.data_write(value=num_fgout_grids,alt_name='num_fgout_grids')
        self.data_write()

        fgno_unset = 0  # to use if fg.fgno not set by user
        fgno_list = []  # to check for uniqueness of fgno's

        for fg in self.fgout_grids:
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from

            if fg.fgno is None:
                # not set by user in setrun
                fgno_unset += 1
                fg.fgno = fgno_unset

            if fg.fgno in fgno_list:
                msg = 'Trying to set fgout grid number to fgno = %i' % fg.fgno \
                      + '\n             but this fgno was already used' \
                      + '\n             Set unique fgno for each fgout grid'
                raise ValueError(msg)

            fgno_list.append(fg.fgno)
            fg.write_to_fgout_data(self._out_file)
        self.close_data_file()

class FGmaxData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(FGmaxData,self).__init__()

        # File name for fgmax points and parameters:
        self.add_attribute('fgmax_files',[])
        self.add_attribute('num_fgmax_val',1)
        self.add_attribute('fgmax_grids',[])


    def write(self,data_source='setrun.py', out_file='fgmax_grids.data'):
        if len(self.fgmax_files) > 0:
            msg = '*** fgmax_files has been deprecated, ' \
                  + 'use fgmax_grids instead.'
            raise ValueError(msg)

        # new style:
        self.open_data_file(out_file, data_source)
        num_fgmax_val = self.num_fgmax_val
        if num_fgmax_val not in [1,2,5]:
            raise NotImplementedError(
                   "Expecting num_fgmax_val in [1,2,5], got %s" % num_fgmax_val)
        self.data_write(value=num_fgmax_val, alt_name='num_fgmax_val')
        num_fgmax_grids = len(self.fgmax_grids)
        self.data_write(value=num_fgmax_grids, alt_name='num_fgmax_grids')
        self.data_write()

        fgno_unset = 0  # to use if fg.fgno not set by user
        fgno_list = []  # to check for uniqueness of fgno's

        for fg in self.fgmax_grids:
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            if fg.xy_fname is not None:
                fg.xy_fname = os.path.abspath(os.path.join(\
                              os.path.dirname(out_file),fg.xy_fname))

            if fg.fgno is None:
                # not set by user in setrun
                fgno_unset += 1
                fg.fgno = fgno_unset

            if fg.fgno in fgno_list:
                msg = 'Trying to set fgmax grid number to fgno = %i' % fg.fgno \
                      + '\n             but this fgno was already used' \
                      + '\n             Set unique fgno for each fgmax grid'
                raise ValueError(msg)

            fgno_list.append(fg.fgno)
            fg.write_to_fgmax_data(self._out_file)
        self.close_data_file()


    def read(self, path="fgmax_grids.data", force=False):
        r"""Read a FGMax data file."""

        super(FGmaxData, self).read(path, force=force)

        # Look for basic parameters
        fig_numbers = []
        with open(os.path.abspath(path), 'r') as data_file:
            # Forward to first parameter
            for line in data_file:
                # Regular parameter setting
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    if varname == "num_fgmax_val":
                        self.num_fgmax_val = int(value)
                    elif varname == "num_fgmax_grids":
                        num_fgmax_grids = int(value)

                # Contains a fixed grid number
                elif "# fgno" in line:
                    value, tail = line.split("#")
                    fig_numbers.append(int(value))

        if len(fig_numbers) != num_fgmax_grids:
            raise ValueError("Number of FGMaxGrid numbers found does not ",
                             "equal the number of grids recorded.")

        # Read each fgmax grid
        import clawpack.geoclaw.fgmax_tools
        for (i, grid_num) in enumerate(fig_numbers):
            new_fgmax_grid = clawpack.geoclaw.fgmax_tools.FGmaxGrid()
            new_fgmax_grid.read_fgmax_grids_data(grid_num, data_file=path)
            self.fgmax_grids.append(new_fgmax_grid)



class DTopoData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(DTopoData,self).__init__()

        # Moving topograhpy
        self.add_attribute('dtopofiles',[])
        self.add_attribute('dt_max_dtopo', 1.e99)

    def write(self, data_source='setrun.py', out_file='dtopo.data'):

        # Moving topography settings
        self.open_data_file(out_file, data_source)
        mdtopofiles = len(self.dtopofiles)
        self.data_write(value=mdtopofiles,alt_name='mdtopofiles')
        self.data_write()
        for tfile in self.dtopofiles:

            if len(tfile) == 4:
                w = '\n  *** WARNING: dtopofile specs changed in v5.8.0 -- ' + \
                      'Flag level info now ignored'
                warnings.warn(w, UserWarning)
                tfile = [tfile[0], tfile[-1]]  # drop minlevel,maxlevel
            elif len(tfile) == 2:
                pass   # now expect only dtopo_type, filename
            else:
                raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
            self._out_file.write("\n'%s' \n" % fname)
            self._out_file.write("%3i   # dtopo_type\n" % tfile[0])
        self.data_write()
        self.data_write(value=self.dt_max_dtopo,alt_name='dt_max_dtopo')
        self.close_data_file()


    def read(self, path="dtopo.data", force=False):
        r"""Read a dtopography data file."""

        print(self.dtopofiles)

        with open(os.path.abspath(path), 'r') as data_file:

            file_name = None

            # Forward to first parameter
            for line in data_file:

                # Regular parameter setting
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    if varname == "mdtopofiles":
                        num_dtopo_files = int(value)
                    elif varname == "dt_max_dtopo":
                        self.dt_max_dtopo = float(value)

                # Assume this is the second line of a record, complete and add
                # to dtopofiles list
                elif file_name is not None:
                    base_values = [int(value) for value in line.split()]
                    base_values.append(file_name)
                    self.dtopofiles.append(base_values)
                    file_name = None

                # Non-empty line, assume start of dtopo_file record
                elif line[0] == "'":
                    file_name = line.strip()[1:-1]

        # Check to make sure we have all the dtopofiles
        if len(self.dtopofiles) != num_dtopo_files:
            raise IOError("The number of dtopo files specified does not equal ",
                          "the number found.")



class ForceDry(clawpack.clawutil.data.ClawData):

    def __init__(self):
        r"""
        A single force_dry array and associated data
        """

        super(ForceDry,self).__init__()
        self.add_attribute('tend',None)
        self.add_attribute('fname','')


class QinitData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(QinitData,self).__init__()

        # Qinit data
        self.add_attribute('qinit_type',0)
        self.add_attribute('qinitfiles',[])
        self.add_attribute('variable_eta_init',False)
        self.add_attribute('force_dry_list',[])
        self.add_attribute('num_force_dry',0)

    def write(self,data_source='setrun.py', out_file='qinit.data'):

        # Initial perturbation
        self.open_data_file(out_file, data_source)
        self.data_write('qinit_type')

        # Perturbation requested
        if self.qinit_type == 0:
            pass
        else:
            # Check to see if each qinit file is present and then write the data
            for tfile in self.qinitfiles:

                if len(tfile) == 3:
                    w = '\n  *** WARNING: qinit specs changed in v5.8.0 -- ' + \
                          'Flag level info now ignored'
                    warnings.warn(w, UserWarning)
                    tfile = [tfile[-1]]  # drop minlevel,maxlevel
                elif len(tfile) == 1:
                    pass  # now expect only filename
                else:
                    raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
                self._out_file.write("\n'%s' \n" % fname)
        # else:
        #     raise ValueError("Invalid qinit_type parameter %s." % self.qinit_type)


        self.data_write('variable_eta_init')

        self.num_force_dry = len(self.force_dry_list)
        self.data_write('num_force_dry')

        for force_dry in self.force_dry_list:

            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),\
                    force_dry.fname))
            self._out_file.write("\n'%s' \n" % fname)
            self._out_file.write("%.3f \n" % force_dry.tend)


        self.close_data_file()


# Storm data
class SurgeData(clawpack.clawutil.data.ClawData):
    r"""Data object describing storm surge related parameters"""

    # Provide some mapping between model names and integers
    storm_spec_dict_mapping = {"data": -1,
                               None: 0,
                               'holland80': 1,
                               'holland08': 8,
                               'holland10': 2,
                               'cle': 3,
                               'slosh': 4,
                               'rankine': 5,
                               'modified-rankine': 6,
                               'DeMaria': 7,
                               'willoughby': 9,
                              }
    storm_spec_not_implemented = ['CLE']

    def __init__(self):
        super(SurgeData, self).__init__()

        # Source term controls
        self.add_attribute('wind_forcing', False)
        self.add_attribute('drag_law', 1)
        self.add_attribute('pressure_forcing', False)
        self.add_attribute('rotation_override', 0)

        # Algorithm parameters - Indexing is python based
        self.add_attribute("wind_index", 4)
        self.add_attribute("pressure_index", 6)
        self.add_attribute("display_landfall_time", False)

        # AMR parameters
        self.add_attribute('wind_refine', [20.0,40.0,60.0])
        self.add_attribute('R_refine', [60.0e3,40e3,20e3])

        # Storm parameters
        self.add_attribute('storm_type', None)  # Backwards compatibility
        self.add_attribute('storm_specification_type', 0) # Type of parameterized storm
        self.add_attribute("storm_file", None) # File containing data

    def write(self, out_file='surge.data', data_source="setrun.py"):
        """Write out the data file to the path given"""

        # print "Creating data file %s" % out_file
        self.open_data_file(out_file,data_source)

        self.data_write('wind_forcing', description='(Wind source term used)')
        self.data_write('drag_law', description='(Type of drag law to use)')
        self.data_write('pressure_forcing',
                        description="(Pressure source term used)")
        if isinstance(self.rotation_override, str):
            if self.rotation_override.lower() == "normal":
                self.rotation_override = 0
            elif "n" in self.rotation_override.lower():
                self.rotation_override = 1
            elif "s" in self.rotation_override.lower():
                self.rotation_override = 2
            else:
                raise ValueError("Unknown rotation_override specification.")
        else:
            self.rotation_override = int(self.rotation_override)
        self.data_write('rotation_override',
                        description="(Override storm rotation)")
        self.data_write()

        self.data_write("wind_index", value=self.wind_index + 1,
                        description="(Index into aux array - fortran indexing)")
        self.data_write("pressure_index", value=self.pressure_index +  1,
                        description="(Index into aux array - fortran indexing)")
        self.data_write("display_landfall_time",
                        description='(Display time relative to landfall)')
        self.data_write()

        if isinstance(self.wind_refine, bool):
            if not self.wind_refine:
                self.data_write('wind_refine', value=False,
                                description='(Refinement ratios)')
        elif isinstance(self.wind_refine, type(None)):
            self.data_write('wind_refine', value=False,
                            description='(Refinement ratios)')
        else:
            self.data_write('wind_refine',description='(Refinement ratios)')
        if isinstance(self.R_refine, bool):
            if not self.R_refine:
                self.data_write('R_refine', value=False,
                                description='(Refinement ratios)')
        elif isinstance(self.R_refine, type(None)):
            self.data_write('R_refine', value=False,
                            description='(Refinement ratios)')
        else:
            self.data_write('R_refine', description='(Refinement ratios)')
        self.data_write()

        # Storm specification
        # Handle deprecated member value
        if self.storm_type is not None:
            self.storm_specification_type = self.storm_type
        # Turn value into integer descriptor
        if isinstance(self.storm_specification_type, int):
            spec_type = self.storm_specification_type
        elif isinstance(self.storm_specification_type, str):
            if self.storm_specification_type.lower() in self.storm_spec_dict_mapping.keys():
                spec_type = self.storm_spec_dict_mapping[self.storm_specification_type.lower()]
            else:
                raise TypeError(f"Unknown storm specification type" +
                                f" '{self.storm_specification_type}' provided.")
        else:
            raise TypeError(f"Unknown storm specification type" +
                            f" '{self.storm_specification_type}' provided.")
        # Check to see if spec type is in supported formats
        if spec_type in self.storm_spec_not_implemented:
            raise NotImplementedError(f"'{spec_type}' has not been implemented.")
        # Write out values
        self.data_write(name="storm_specification_type", value=spec_type,
                        description="(Storm specification)")
        self.data_write(name="storm_file", description='(Path to storm data)')

        self.close_data_file()


class FrictionData(clawpack.clawutil.data.ClawData):
    r"""Data class representing complex variable friction"""

    def __init__(self):
        r""""""

        super(FrictionData, self).__init__()

        # Variable friction support
        self.add_attribute('variable_friction', False)

        # Index where the variable friction field is stored (Python indexed)
        self.add_attribute('friction_index', 3)

        # Region support
        self.add_attribute('friction_regions', [])

        # File support
        self.add_attribute('friction_files', [])


    def read(self, path="friction.data", force=False):
        r"""Read friction data file"""

        with open(os.path.abspath(path), 'r') as data_file:
            # Header
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()

            # Generic data
            self.variable_friction = bool(data_file.readline().split("=:")[0])
            self.friction_index = int(data_file.readline().split("=:")[0])
            data_file.readline()
            num_regions = int(data_file.readline().split("=:")[0])
            data_file.readline()
            # Regions
            self.friction_regions = []
            for n in range(num_regions):
                lower = self._convert_line(data_file.readline())
                upper = self._convert_line(data_file.readline())
                depths = self._convert_line(data_file.readline())
                coeff = self._convert_line(data_file.readline())
                self.friction_regions.append([lower, upper, depths, coeff])
                data_file.readline()
            self.friction_files = [] # Is not supported

    def _convert_line(self, line):
        values = []
        for value in line.split("=:")[0].split(" "):
            if len(value) > 1:
                values.append(float(value))
        return values


    def write(self, out_file='friction.data', data_source='setrun.py'):

        self.open_data_file(out_file, data_source)

        self.data_write('variable_friction',
                        description="(method for setting variable friction)")
        self.data_write('friction_index', value=self.friction_index + 1,
                        description=("(Index into aux array ",
                                     "- fortran indexing)"))
        self.data_write()
        if self.variable_friction:
            # Region based friction
            self.data_write(value=len(self.friction_regions),
                            alt_name='num_friction_regions',
                            description="(Friction Regions)")
            self.data_write()
            for region in self.friction_regions:
                self.data_write(value=region[0], alt_name="lower")
                self.data_write(value=region[1], alt_name="upper")
                self.data_write(value=region[2], alt_name="depths")
                self.data_write(value=region[3],
                                alt_name="manning_coefficients")
                self.data_write()

            # File based friction
            self.data_write(value=len(self.friction_files),
                            alt_name='num_friction_files')
            for friction_file in self.friction_files:
                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),
                                                     friction_file))
                self._out_file.write("'%s' %s\n " % fname)

        self.close_data_file()


class MultilayerData(clawpack.clawutil.data.ClawData):
    r"""
    Multilayer SWE data object
    """

    def __init__(self):
        super(MultilayerData, self).__init__()

        # Physics parameters
        self.add_attribute('num_layers', 1)
        self.add_attribute('rho', [1025.0, 1028.0])
        self.add_attribute('eta', [0.0, -200.0])
        self.add_attribute('wave_tolerance', [1.e-1, 1.e-1])

        # Algorithm parameters
        self.add_attribute('eigen_method', 4)
        self.add_attribute('inundation_method', 2)
        self.add_attribute('check_richardson', True)
        self.add_attribute('richardson_tolerance', 0.95)
        self.add_attribute('layer_index', 8)

        # Need to adjust refinement module for this, dry_limit is in geodata
        self.add_attribute('wave_tolerance', [1e-1, 2e-1])
        self.add_attribute('dry_limit', False)

    def write(self, out_file='multilayer.data', datasource="setrun.py"):

        self.open_data_file(out_file, datasource)

        self.data_write('num_layers', description='(Number of layers)')
        self.data_write('eta',
                        description='(Initial top surface of each layer)')
        self.data_write('wave_tolerance',
                        description=('(Tolerance of surface perturbation per',
                                     ' layer, used for refinement criteria)'))
        self.data_write('layer_index', value=self.layer_index + 1,
                        description=("(Index into aux array -",
                                     " fortran indexing)"))
        self.data_write(None)
        self.data_write('check_richardson',
                        description="(Check Richardson number)")
        self.data_write('richardson_tolerance',
                        description='(Tolerance for Richardson number)')
        self.data_write('eigen_method',
                        description='(Method for calculating eigenspace)')
        self.data_write('inundation_method',
                        description=('(Method for calculating inundation ',
                                     'eigenspace)'))
        self.close_data_file()



class BoussData(clawpack.clawutil.data.ClawData):
    r"""
     data object for Boussinesq info in 2D geoclaw

    """
    def __init__(self):
        super(BoussData,self).__init__()

        self.add_attribute('bouss_equations',2)
        self.add_attribute('bouss_min_level', 1)
        self.add_attribute('bouss_max_level', 10)
        self.add_attribute('bouss_min_depth', 10.)
        self.add_attribute('bouss_solver', 3)
        self.add_attribute('bouss_tstart', 0.)

    def write(self,out_file='bouss.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)
        self.data_write('bouss_equations', description='0=SWE, 1=MS, 2=SGN')
        self.data_write('bouss_min_level',
                        description='coarsest level to apply bouss')
        self.data_write('bouss_max_level',
                        description='finest level to apply bouss')
        self.data_write('bouss_min_depth',
                        description='depth to switch to SWE')
        self.data_write('bouss_solver', description='1=GMRES, 2=Pardiso, 3=PETSc')
        self.data_write('bouss_tstart', description='time to switch from SWE')

        self.close_data_file()


# ==================================
# data objects for 1d_classic code
# ==================================


#  Gauge data object removed, version from amrclaw works in 1d
#class GaugeData1D(clawpack.clawutil.data.ClawData):


class GridData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for grid info

    """
    def __init__(self):
        super(GridData1D,self).__init__()

        self.add_attribute('grid_type',0)
        self.add_attribute('fname_celledges',None)
        self.add_attribute('monitor_fgmax',False)
        self.add_attribute('monitor_runup',False)
        self.add_attribute('monitor_total_zeta',False)

    def write(self,out_file='grid.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('grid_type')
        if self.grid_type == 2:
            if self.fname_celledges is None:
                self.fname_celledges = 'celledges.txt'
                print('*** grid_type ==2 and fname_celledges not specified,')
                print('*** using celledges.txt')
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),
                                    self.fname_celledges))
            self._out_file.write("\n'%s'   =: fname_celledges\n " % fname)

        self._out_file.write("\n%s   =: monitor_fgmax" \
                             % str(self.monitor_fgmax)[0])
        self._out_file.write("\n%s   =: monitor_runup" \
                             % str(self.monitor_runup)[0])
        self._out_file.write("\n%s   =: monitor_total_zeta" \
                             % str(self.monitor_total_zeta)[0])
        self.close_data_file()

    def read(self, path, force=False):
        with open(os.path.abspath(path), 'r') as data_file:
            for line in data_file:
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]
                    if varname == 'grid_type':
                        self.grid_type = int(value)
                    elif varname == 'fname_celledges':
                        self.fname_celledges = value.strip()


class BoussData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for Boussinesq info

    """
    def __init__(self):
        super(BoussData1D,self).__init__()

        self.add_attribute('bouss_equations',2)
        self.add_attribute('bouss_min_depth',20.)

    def write(self,out_file='bouss.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('bouss_equations')
        self.data_write('bouss_min_depth')

        self.close_data_file()

