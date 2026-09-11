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
from pathlib import Path
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
        self.add_attribute('speed_tolerance', None)
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

        if isinstance(self.speed_tolerance, float):
            self.speed_tolerance = [self.speed_tolerance]
        self.data_write('speed_tolerance')
        self.data_write()
        self.data_write('variable_dt_refinement_ratios',
                        description="(Set dt refinement ratios automatically)")
        self.close_data_file()



def _reject_remote_path(path, kind, fetch_hint):
    """Raise if *path* is a URL rather than a local file.

    GeoClaw's Fortran reader opens a filesystem path, and the ``.data`` writers
    run every path through ``os.path.abspath``, which turns a URL into a bogus
    local path (see ``netcdf_utils.is_remote_url``).  The result used to be a
    FileNotFoundError naming a path the user never typed, which gives no hint
    that the real problem is "remote sources must be fetched first".

    *kind* is the noun for the message ("topography"/"dtopography") and
    *fetch_hint* the recipe to show.
    """
    from clawpack.geoclaw.netcdf_utils import is_remote_url

    if is_remote_url(path):
        raise ValueError(
            f"{kind} path is a URL, which GeoClaw's Fortran reader cannot "
            f"open:\n    {path}\n"
            f"Read it in Python first and write a local file, then reference "
            f"that:\n{fetch_hint}")


def _write_preprocessing_block(f, t):
    """Write the 8 preprocessing-attribute lines for one topo/dtopo file.

    The lines are, in order: crop_extent, coarsen, buffer, align, x_shift,
    y_shift, z_shift, negate_z.

    Shared by TopographyData.write() (topo.data) and DTopoData.write()
    (dtopo.data); Fortran reads the same 8 lines in read_topo_settings and
    read_dtopo_settings.  *t* is a Topography or DTopography object.

    Float values use repr (shortest round-trip representation) so coordinates
    reach Fortran at full precision; %g would truncate to 6 significant
    digits.
    """
    # crop_extent: write all-zero sentinel or 4 space-separated floats
    # All-zero is safe: a valid extent requires x1<x2 and y1<y2.
    if t.crop_extent is None:
        f.write("0. 0. 0. 0.   # crop_extent [x1 x2 y1 y2]\n")
    else:
        vals = ' '.join(repr(float(v)) for v in t.crop_extent)
        f.write(f"{vals}   # crop_extent [x1 x2 y1 y2]\n")

    f.write(f"{int(t.coarsen):d}   # coarsen\n")
    # buffer is a count of grid points; Fortran reads an integer.
    f.write(f"{int(t.buffer):d}   # buffer\n")

    # align: write all-zero sentinel or 2 space-separated floats
    if t.align is None:
        f.write("0. 0.   # align [x y]\n")
    else:
        vals = ' '.join(repr(float(v)) for v in t.align)
        f.write(f"{vals}   # align [x y]\n")

    f.write(f"{float(t.x_shift)!r}   # x_shift\n")
    f.write(f"{float(t.y_shift)!r}   # y_shift\n")
    f.write(f"{float(t.z_shift)!r}   # z_shift\n")
    f.write(f"{'T' if t.negate_z else 'F'}   # negate_z\n")


class TopographyData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(TopographyData,self).__init__()

        # Topography data
        self.add_attribute('topo_missing', 99999.0)
        self.add_attribute('test_topography', 0)
        self.add_attribute('override_order', False)
        self.add_attribute('topofiles', [])

        # Run coordinate system (mirrors geodata.coordinate_system: 1 =
        # Cartesian, 2 = lon-lat).  Registered but never data_write-n, so it is
        # NOT emitted to topo.data -- it only gates the Python-side antimeridian-
        # wrap early check at write time (the Fortran read applies the
        # authoritative gate).  None => per-file heuristic, no cross-check.
        self.add_attribute('coordinate_system', None)

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


    def _normalize_topofiles(self):
        """Return topofiles as a list of Topography objects.

        Converts legacy list/tuple and dict entries with a DeprecationWarning.
        Topography instances pass through unchanged.
        """
        from clawpack.geoclaw.topotools import Topography

        result = []
        for entry in self.topofiles:
            if isinstance(entry, Topography):
                result.append(entry)
                continue

            if isinstance(entry, (list, tuple)):
                # Handle old 6-element format [topo_type,minlev,maxlev,t1,t2,path]
                if len(entry) == 6:
                    warnings.warn(
                        "6-element topofile entries are deprecated since v5.8.0 "
                        "(level/time info is ignored). Construct Topography "
                        "objects directly instead.",
                        DeprecationWarning, stacklevel=3,
                    )
                    entry = [entry[0], entry[-1]]

                # [topo_type, path, TopoMetadata] — from topo_entries(); not deprecated
                if len(entry) == 3:
                    from clawpack.geoclaw.netcdf_utils import TopoMetadata
                    if isinstance(entry[2], TopoMetadata):
                        topo = Topography()
                        topo.topo_type = int(entry[0])
                        topo.path = str(entry[1])
                        topo._netcdf_meta = entry[2]
                        result.append(topo)
                        continue

                # [topo_type, path] — legacy deprecated format
                if len(entry) >= 2:
                    warnings.warn(
                        "List/tuple topofile entries are deprecated. Construct "
                        "Topography objects directly and append them to "
                        "rundata.topo_data.topofiles instead.",
                        DeprecationWarning, stacklevel=3,
                    )
                    topo = Topography()
                    topo.topo_type = int(entry[0])
                    topo.path = str(entry[1])
                    result.append(topo)
                    continue

            if isinstance(entry, dict):
                warnings.warn(
                    "Dict topofile entries are deprecated. Construct Topography "
                    "objects directly and set attributes on them instead. "
                    "Note: dict key 'extent' maps to Topography.crop_extent.",
                    DeprecationWarning, stacklevel=3,
                )
                topo = Topography()
                topo.path = entry.get('topo_path', None)
                raw_type = entry.get('topo_type', None)
                if raw_type is not None:
                    topo.topo_type = int(raw_type)
                # The legacy dict key 'extent' is an alias for 'crop_extent'
                # (the requested crop; see Topography "Region terminology").
                # Note: if a spec supplies BOTH 'extent' and 'crop_extent', the
                # 'crop_extent' key wins -- it is applied by the loop below,
                # which runs after this alias assignment.
                if 'extent' in entry:
                    topo.crop_extent = entry['extent']
                for attr in ('crop_extent', 'coarsen', 'buffer', 'align',
                             'x_shift', 'y_shift', 'z_shift', 'negate_z'):
                    if attr in entry:
                        setattr(topo, attr, entry[attr])
                result.append(topo)
                continue

            raise ValueError(
                f"Unrecognized topofile entry type {type(entry).__name__!r}: "
                f"{entry!r}.  Expected a Topography object, a [topo_type, path] "
                f"list, or a dict."
            )
        return result

    def _compute_priority_order(self, topos):
        """Return *topos* sorted coarsest-first (finest = highest priority last).

        Coarsest resolution (largest dx*dy) is placed at index 0 and written
        first in topo.data; the finest file is written last.  Fortran assigns
        rank 1 (highest priority) to the *last* file listed, so last-written =
        finest = highest priority in overlap resolution.  This matches the
        traditional GeoClaw convention of listing topography files from
        coarsest to finest.

        Stable sort: equal-area files preserve their relative input order, so
        among equal-resolution files the last one listed wins.

        When override_order=True the input list is returned unchanged.  The
        caller is responsible for placing the highest-priority (finest) file
        last.  Files are written in the order provided; the last file has the
        highest priority in GeoClaw.

        Header data is loaded via read_header() if coordinates are not yet
        available; the Z array is never loaded here.
        """
        if self.override_order:
            return list(topos)

        def _cell_area(topo):
            if topo.path is None:
                return float('inf')
            if topo._x is None:
                try:
                    topo.read_header()
                except Exception:
                    return float('inf')
            try:
                dx, dy = topo.delta
                return float(dx) * float(dy)
            except Exception:
                return float('inf')

        # Descending cell area => coarsest first, finest last.  reverse=True
        # keeps the sort stable for equal-area files (they retain input order).
        return sorted(topos, key=_cell_area, reverse=True)

    def _resolve_topo_records(self, topos, out_file):
        """Expand *topos* into the entries that will be written to topo.data.

        Returns a list of ``(fname, topo_type, topo, meta)``.  ``meta`` is a
        ``TopoMetadata`` for type-4 files (written as a CF descriptor block)
        and None otherwise.

        Most files produce exactly one record.  A ``topo_type=4`` file whose
        ``crop_extent`` runs past the file's longitude extent produces one
        record per side of the antimeridian, each with its own
        ``lon_wrap_offset`` -- this is what makes a cross-seam crop work from
        an ordinary ``topofiles.append(topo)`` instead of requiring the caller
        to build descriptors by hand.
        """
        import dataclasses
        from clawpack.geoclaw import netcdf_utils as _ncutils
        from clawpack.geoclaw import topotools

        records = []
        for topo in topos:
            _reject_remote_path(
                topo.path, "Topography",
                "    topo = topotools.fetch_remote_topo(\n"
                "        url, crop_extent=[...], coarsen=..., buffer=...)\n"
                "    topo.write('topo_cropped.tt3', topo_type=3)\n"
                "    rundata.topo_data.topofiles.append(topo)\n"
                "Only the requested hyperslab is read, so this does not "
                "download the whole file.")

            # Resolve path relative to out_file's directory, same as before
            fname = os.path.abspath(
                os.path.join(os.path.dirname(out_file), topo.path))

            # topo_type may still be None when a Topography was built by hand
            # and never read; the :3d format would raise a TypeError naming
            # neither the file nor the attribute.
            topo_type = topo.topo_type
            if topo_type is None:
                topo_type = topotools.determine_topo_type(topo.path,
                                                          default=None)
                if topo_type is None:
                    raise ValueError(
                        f"topo_type is not set for {topo.path} and cannot be "
                        f"inferred from its extension. Set topo.topo_type "
                        f"explicitly, or pass topo_type= to Topography().")
                topo.topo_type = topo_type

            # Type 5 (GeoTIFF) is readable in Python but has no case(5) in the
            # Fortran reader, which aborts with "Unrecognized topo_type".
            if abs(topo_type) == 5:
                warnings.warn(
                    f"{topo.path} is written to {out_file} as topo_type=5 "
                    f"(GeoTIFF). GeoClaw's Fortran reader does not support "
                    f"type 5 and will abort; convert to topo_type 3 or 4 "
                    f"(Topography.write) before running.",
                    UserWarning, stacklevel=3)

            # A descending crop_extent is not a rectangle in any frame.  The
            # *continuous* spelling of a cross-seam crop ([-190, -120]) is
            # handled below; the wrapped spelling ([170, -170]) cannot be, and
            # written directly it would emit "crop_bounds = 170.0 -170.0",
            # which Fortran resolves to mx=0, my=0 -- an empty topo, silently.
            if (topo.crop_extent is not None
                    and topo.crop_extent[0] >= topo.crop_extent[1]
                    and getattr(topo, '_netcdf_meta', None) is None):
                raise ValueError(
                    f"crop_extent {list(topo.crop_extent)} for {topo.path} is "
                    f"descending in longitude. To cross the antimeridian, "
                    f"write the crop in continuous coordinates instead -- "
                    f"[{topo.crop_extent[0] - 360.0}, {topo.crop_extent[1]}] "
                    f"for this one -- which is split across the seam "
                    f"automatically. A descending pair has no such reading and "
                    f"would produce an empty grid in Fortran with no error.")

            if abs(topo_type) != 4:
                records.append((fname, topo_type, topo, None))
                continue

            # --- type 4: build the CF descriptor block -------------------
            if getattr(topo, '_netcdf_meta', None) is not None:
                # Pre-computed by topo_entries(); already carries the right
                # lon_wrap_offset and file-coordinate crop_bounds.
                records.append((fname, topo_type, topo, topo._netcdf_meta))
                continue

            crop = (tuple(float(v) for v in topo.crop_extent)
                    if topo.crop_extent is not None else None)

            # Opened without crop_bounds so the longitude extent can be read
            # before deciding whether the crop needs wrapping; setting them at
            # construction would validate (and reject) a cross-seam crop first.
            with _ncutils.TopoInspector(fname) as insp:
                if insp.var_name is None:
                    insp.var_name = insp._find_topo_var_name()

                needs_wrap = False
                if crop is not None:
                    x_name = insp._find_x_name()
                    lon = insp.ds[x_name].values
                    file_lon_min = float(lon.min())
                    file_lon_max = float(lon.max())
                    tol = 1e-9
                    needs_wrap = (crop[0] < file_lon_min - tol
                                  or crop[1] > file_lon_max + tol)

                if not needs_wrap:
                    # Unchanged path: validates crop_bounds (including
                    # latitude) and skips the expensive fill scan, exactly as
                    # before.
                    insp.crop_bounds = crop
                    file_meta = insp.inspect(insp.var_name)
                    # Gate on the run's coordinate system.  This path never
                    # wraps (lon_wrap_offset stays 0.0), but a geographic file
                    # under a Cartesian run -- or the reverse -- is a setup
                    # error worth catching here rather than in Fortran.
                    _xc = insp.ds[file_meta.x_name]
                    _ncutils.resolve_wrap(
                        self.coordinate_system,
                        _ncutils.classify_lon_axis(
                            _xc.attrs.get('units'),
                            _xc.attrs.get('standard_name'),
                            float(_xc.min()), float(_xc.max())))
                    src_units = insp._check_topo_units()
                    scale = _ncutils._units_scale(
                        src_units, _ncutils.GEOCLAW_NETCDF_UNITS['topo'])
                    records.append((fname, topo_type, topo,
                                    _ncutils.TopoMetadata(
                                        **dataclasses.asdict(file_meta),
                                        var_name=insp.var_name,
                                        source_units=src_units,
                                        scale_factor=scale,
                                        fill_action='abort',
                                        lon_wrap_offset=0.0)))
                    continue

                # Cross-seam (or wholly off-seam) crop: one entry per side,
                # each with its own lon_wrap_offset.  fill_scan=False because
                # topo_entries inspects with crop_bounds unset, which would
                # otherwise scan the whole file -- ruinous for a global DEM
                # and prone to rejecting NaN far outside the crop.
                y_name = insp._find_y_name()
                lat = insp.ds[y_name].values
                file_lat_min = float(lat.min())
                file_lat_max = float(lat.max())
                if (crop[2] < file_lat_min - 1e-9
                        or crop[3] > file_lat_max + 1e-9):
                    # Longitude wraps; latitude never does.
                    raise ValueError(
                        f"crop_extent {list(topo.crop_extent)} for {topo.path} "
                        f"exceeds the file's latitude extent "
                        f"[{file_lat_min}, {file_lat_max}]. Longitude is "
                        f"wrapped across the antimeridian, but latitude cannot "
                        f"be; narrow the requested latitude range.")

                insp.crop_bounds = crop
                entries = insp.topo_entries(
                    fill_scan=False,
                    coordinate_system=self.coordinate_system)

            for _entry_type, _entry_path, entry_meta in entries:
                records.append((fname, topo_type, topo, entry_meta))

        return records

    def write(self, data_source='setrun.py', out_file='topo.data'):

        self.open_data_file(out_file, data_source)
        self.data_write(name='topo_missing',
                        description='replace no_data_value in topofile')
        self.data_write(name='test_topography',
                        description='(Type topography specification)')
        if self.test_topography == 0:
            topos = self._normalize_topofiles()
            topos = self._compute_priority_order(topos)

            # Warn if the topo files carry mismatched vertical datums: GeoClaw
            # applies no vertical-datum transformation, so mixing references is
            # a likely error.  Only datums that are actually available (set by
            # the user or read from NetCDF) are compared; None is ignored.
            datums = sorted({str(topo.datum) for topo in topos
                             if getattr(topo, 'datum', None) is not None})
            if len(datums) > 1:
                warnings.warn(
                    "Topography files have mismatched vertical datums (%s). "
                    "GeoClaw does not transform between datums; verify the "
                    "files share a common vertical reference."
                    % ", ".join(datums),
                    UserWarning, stacklevel=2,
                )

            # Resolve first, write second.  A type-4 crop that runs off the
            # file's longitude extent is covered by *two* descriptor entries
            # (one per side of the seam), so the entry count is not known until
            # every file has been resolved -- and ntopofiles is written before
            # the blocks.
            records = self._resolve_topo_records(topos, out_file)

            self.data_write(value=len(records), alt_name='ntopofiles')
            f = self._out_file
            for fname, topo_type, topo, meta in records:
                f.write(f"\n'{fname}'   # topo_path\n")
                f.write(f"{topo_type:3d}   # topo_type\n")
                # The originating Topography is reused for every entry it
                # expanded into, so buffer/coarsen/align/shifts reach Fortran
                # for each one.  (crop_extent is written as the user gave it;
                # for type 4 the descriptor's crop_bounds takes priority, and
                # for a wrapped pair it is the per-entry crop_bounds that
                # differ.)
                _write_preprocessing_block(f, topo)
                if meta is not None:
                    from clawpack.geoclaw import netcdf_utils as _ncutils
                    _ncutils.DescriptorWriter.write_topo_descriptor(f, meta)

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

    def _normalize_dtopofiles(self):
        """Return dtopofiles as a list of DTopography objects.

        Converts legacy list/tuple entries with a warning.  DTopography
        instances pass through unchanged.
        """
        from clawpack.geoclaw.dtopotools import DTopography

        result = []
        for entry in self.dtopofiles:
            if isinstance(entry, DTopography):
                result.append(entry)
                continue

            if isinstance(entry, (list, tuple)):
                if len(entry) == 4:
                    w = '\n  *** WARNING: dtopofile specs changed in ' + \
                        'v5.8.0 -- Flag level info now ignored'
                    warnings.warn(w, UserWarning)
                    entry = [entry[0], entry[-1]]  # drop minlevel,maxlevel
                if len(entry) == 2:
                    d = DTopography()
                    d.dtopo_type = int(entry[0])
                    d.path = str(entry[1])
                    result.append(d)
                    continue
                raise ValueError('Unexpected len(tfile) = %i' % len(entry))

            raise ValueError(
                f"Unrecognized dtopofile entry type {type(entry).__name__!r}: "
                f"{entry!r}.  Expected a DTopography object or a "
                f"[dtopo_type, path] list."
            )
        return result

    def write(self, data_source='setrun.py', out_file='dtopo.data'):

        # Moving topography settings
        self.open_data_file(out_file, data_source)
        dtopos = self._normalize_dtopofiles()
        self.data_write(value=len(dtopos), alt_name='mdtopofiles')
        self.data_write()
        for d in dtopos:
            # Unsupported preprocessing attributes fail loudly at write time
            # (mirrors DTopography.read() and the Fortran guard in
            # read_dtopo_settings): only z_shift and negate_z are implemented
            # for dtopography.
            unsupported = [name for name, is_set in (
                ("crop_extent", d.crop_extent is not None),
                ("coarsen", d.coarsen != 1),
                ("buffer", d.buffer != 0),
                ("align", d.align is not None),
            ) if is_set]
            if unsupported:
                raise NotImplementedError(
                    "Preprocessing attributes %s are not implemented for "
                    "dtopography (file %s). Only x_shift, y_shift, z_shift and "
                    "negate_z are supported." % (", ".join(unsupported), d.path))

            _reject_remote_path(
                d.path, "DTopography",
                "    dtopo = dtopotools.DTopography()\n"
                "    dtopo.read(url, dtopo_type=4)\n"
                "    dtopo.write('dtopo_local.tt3', dtopo_type=3)\n"
                "    rundata.dtopo_data.dtopofiles.append(dtopo)")

            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(
                os.path.join(os.path.dirname(out_file), d.path))
            self._out_file.write("\n'%s'   # dtopo_path\n" % fname)
            self._out_file.write("%3i   # dtopo_type\n" % d.dtopo_type)
            _write_preprocessing_block(self._out_file, d)

            # For NetCDF (type 4): write the descriptor block that Fortran's
            # read_dtopo_netcdf_descriptor parses (key=value lines terminated
            # by a blank line), carrying the time axis as (t0, dt) in
            # simulation seconds.
            if abs(d.dtopo_type) == 4:
                from clawpack.geoclaw import netcdf_utils as _ncutils
                from clawpack.geoclaw.topotools import extract_datum
                with _ncutils.DTopoInspector(fname) as _insp:
                    # A recognized non-meter deformation unit yields a
                    # scale_factor Fortran applies on read.
                    _meta = _insp.inspect_dtopo()
                    # Record an optional vertical datum (informational) for the
                    # consistency check below, if not already set on the object.
                    if d.datum is None:
                        d.datum = extract_datum(
                            *[_insp.ds[v].attrs for v in _insp.ds.data_vars],
                            _insp.ds.attrs)
                _ncutils.DescriptorWriter.write_dtopo_descriptor(
                    self._out_file, _meta)

        # Warn if the dtopo files carry mismatched vertical datums: GeoClaw
        # applies no vertical-datum transformation, so mixing references is a
        # likely error.  Only datums that are actually available are compared.
        datums = sorted({str(d.datum) for d in dtopos
                         if getattr(d, 'datum', None) is not None})
        if len(datums) > 1:
            warnings.warn(
                "dtopo files have mismatched vertical datums (%s). GeoClaw "
                "does not transform between datums; verify the files share a "
                "common vertical reference." % ", ".join(datums),
                UserWarning, stacklevel=2,
            )

        self.data_write()
        self.data_write(value=self.dt_max_dtopo,alt_name='dt_max_dtopo')
        self.close_data_file()


    def read(self, path="dtopo.data", force=False):
        r"""Read a dtopography data file written by write().

        Populates *dtopofiles* with DTopography objects carrying the path,
        dtopo_type, and preprocessing attributes of each 10-line file block.
        """
        from clawpack.geoclaw.dtopotools import DTopography

        def _data(line):
            # Strip trailing comments and whitespace
            return line.split('#')[0].strip()

        with open(os.path.abspath(path), 'r') as data_file:
            lines = data_file.readlines()

        num_dtopo_files = None
        self.dtopofiles = []
        i = 0
        while i < len(lines):
            line = lines[i]
            if "=:" in line:
                value, tail = line.split("=:")
                varname = tail.split()[0]
                if varname == "mdtopofiles":
                    num_dtopo_files = int(value)
                elif varname == "dt_max_dtopo":
                    self.dt_max_dtopo = float(value)
                i += 1
            elif _data(line).startswith("'"):
                # 10-line file block: path, dtopo_type, crop_extent, coarsen,
                # buffer, align, x_shift, y_shift, z_shift, negate_z
                d = DTopography()
                d.path = _data(line).strip("'")
                d.dtopo_type = int(_data(lines[i + 1]).split()[0])
                crop = [float(v) for v in _data(lines[i + 2]).split()]
                d.crop_extent = None if all(v == 0. for v in crop) else crop
                d.coarsen = int(_data(lines[i + 3]))
                d.buffer = int(_data(lines[i + 4]))  # grid-point count
                align = [float(v) for v in _data(lines[i + 5]).split()]
                d.align = None if all(v == 0. for v in align) else align
                d.x_shift = float(_data(lines[i + 6]))
                d.y_shift = float(_data(lines[i + 7]))
                d.z_shift = float(_data(lines[i + 8]))
                d.negate_z = _data(lines[i + 9]).upper().startswith("T")
                self.dtopofiles.append(d)
                i += 10
            else:
                i += 1

        # Check to make sure we have all the dtopofiles
        if len(self.dtopofiles) != num_dtopo_files:
            raise IOError("The number of dtopo files specified does not equal "
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


# =============================================================================
#  Meteorological forcing selection registry
#
#  Single source of truth mapping the met-forcing *family* + *subtype* to the
#  legacy signed-integer ``storm_specification_type`` encoding still carried on
#  the ``surge.data`` wire.  See met_forcing_refactor.md Sections 7 and 9.
#
#  family:  "parametric" - analytic model with a storm center/track;
#           "gridded"    - file-backed wind/pressure fields (no center);
#           "none"       - forcing off.
#
#  Each parametric subtype maps to that model's historical integer code; the
#  single gridded subtype maps to -1 and "none" to 0.  The gridded ASCII/NetCDF
#  distinction is NOT on this wire - it lives as ``file_format`` in the .storm
#  descriptor (see surge/gridded.py and gridded_met_forcing_module.f90).

# canonical subtype -> (family, legacy integer code)
forcing_subtype_registry = {"holland80":        ("parametric",  1),
                            "holland2008":      ("parametric",  8),
                            "holland2010":      ("parametric",  2),
                            "cle":              ("parametric",  3),
                            "slosh":            ("parametric",  4),
                            "rankine":          ("parametric",  5),
                            "modified_rankine": ("parametric",  6),
                            "demaria":          ("parametric",  7),
                            "willoughby":       ("parametric",  9),
                            "gridded":          ("gridded",    -1),
                            "none":             ("none",        0),
                           }

# Historical / alternate spellings -> canonical subtype.  Lookups are
# case-folded, so only spellings that differ by more than case need listing
# here.  These keep existing setrun.py scripts working unchanged.
forcing_subtype_aliases = {"holland08":        "holland2008",
                           "holland10":        "holland2010",
                           "modified-rankine": "modified_rankine",
                           "data":             "gridded",
                          }

# Subtypes recognized by the API but not yet implemented in the Fortran forcing
# code.  Currently empty: all nine parametric models and the gridded path have
# working Fortran implementations.  (The previous ``storm_spec_not_implemented``
# check compared an int against the string 'CLE' and so never fired; CLE is in
# fact implemented in parametric_met_forcing_module.f90, so nothing is blocked here.)
forcing_not_implemented = set()

# Reverse map: legacy integer code -> canonical subtype (each code is unique).
_forcing_code_to_subtype = {code: subtype for subtype, (family, code)
                            in forcing_subtype_registry.items()}


def resolve_forcing_subtype(value):
    r"""Resolve a user-supplied forcing spec to ``(family, subtype, code)``.

    ``value`` may be a canonical subtype string, a historical alias, a legacy
    integer code, or ``None`` (forcing off).  Returns the canonical family
    string, canonical subtype string, and legacy integer code.
    """
    if value is None:
        subtype = "none"
    elif isinstance(value, str):
        subtype = value.lower()
        subtype = forcing_subtype_aliases.get(subtype, subtype)
        if subtype not in forcing_subtype_registry:
            raise ValueError(f"Unknown forcing subtype '{value}'.")
    elif isinstance(value, (int, np.integer)):
        if int(value) not in _forcing_code_to_subtype:
            raise ValueError(f"Unknown forcing specification code '{value}'.")
        subtype = _forcing_code_to_subtype[int(value)]
    else:
        raise TypeError(f"Unknown forcing specification '{value}'.")

    if subtype in forcing_not_implemented:
        raise NotImplementedError(f"Forcing subtype '{subtype}' is not yet "
                                  "implemented.")

    family, code = forcing_subtype_registry[subtype]
    return family, subtype, code


# Storm data
class SurgeData(clawpack.clawutil.data.ClawData):
    r"""Meteorological (storm) forcing parameters written to ``surge.data``.

    Set as ``rundata.surge_data`` in a ``setrun.py`` and written to
    ``surge.data``, which the Fortran ``met_forcing_module`` reads to activate
    and configure wind/pressure forcing.  See :ref:`setrun_surge` for the
    per-attribute reference.

    :Forcing selection:
     The forcing family and subtype are the preferred, explicit way to select a
     model:

     - ``storm_family`` -- ``"parametric"`` (an analytic model with a storm
       center/track), ``"gridded"`` (file-backed wind/pressure fields, e.g.
       OWI/ASCII or NetCDF), or ``"none"`` (forcing off).
     - ``storm_subtype`` -- for a parametric family, a model name such as
       ``"holland80"``, ``"holland2010"``, ``"cle"``, ``"slosh"``,
       ``"rankine"``, ``"modified_rankine"``, ``"demaria"``, or
       ``"willoughby"``; for a gridded family, ``"gridded"``.

     The legacy ``storm_specification_type`` (a model-name string or the signed
     integer code, e.g. ``"holland80"``/``1`` or ``"data"``/``-1``) remains
     fully supported: when ``storm_family``/``storm_subtype`` are unset it is
     resolved through :data:`forcing_subtype_registry`.  All three map to the
     same forcing on the ``surge.data`` wire.

    :Key attributes:
     - ``wind_forcing`` / ``pressure_forcing`` (bool) -- enable the wind and
       pressure source terms.
     - ``drag_law`` (int) -- wind-drag law (0 none, 1 Garratt, 2 Powell).
     - ``wind_index`` / ``pressure_index`` (int) -- 0-based ``aux`` component
       indices for the forcing fields (Fortran indexing is +1).
     - ``storm_time_scale`` (float) -- multiplicative scale on the storm time
       axis (>1 slower, <1 faster).
     - ``t_ramp_on`` / ``t_ramp_off`` (float) -- seconds over which the forcing
       ramps on after ``t0`` and off before ``tfinal`` (0 disables).
     - ``rotation_override`` (int/str) -- override the hemisphere-based storm
       rotation sense.
     - ``wind_refine`` / ``R_refine`` (list/bool) -- AMR refinement thresholds
       on wind speed and (parametric only) distance to the storm center.
     - ``storm_file`` (str) -- path to the storm track file or gridded
       descriptor.

    See :ref:`storm_module` for the Python storm/track object model and
    :ref:`surgedata` for storm data sources.
    """

    # Legacy name/alias -> integer mapping, derived from the canonical
    # forcing_subtype_registry above.  Retained (with the historical spellings)
    # as a backwards-compatibility view; new code should prefer
    # storm_family/storm_subtype and resolve_forcing_subtype().
    storm_spec_dict_mapping = {None: 0,
                               **{sub: code for sub, (fam, code)
                                  in forcing_subtype_registry.items()},
                               **{alias: forcing_subtype_registry[canon][1]
                                  for alias, canon
                                  in forcing_subtype_aliases.items()},
                              }
    storm_spec_not_implemented = forcing_not_implemented

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

        # Model-storm time scaling and temporal onset/cutoff ramp.
        # storm_time_scale: >1 slower, <1 faster (model storms; data storms
        #   carry their own scale in the .storm descriptor).
        # t_ramp_on / t_ramp_off: seconds over which the wind/pressure forcing
        #   ramps on after t0 and off before tfinal (0 = no ramp).
        self.add_attribute("storm_time_scale", 1.0)
        self.add_attribute("t_ramp_on", 0.0)
        self.add_attribute("t_ramp_off", 0.0)

        # AMR parameters
        self.add_attribute('wind_refine', [20.0,40.0,60.0])
        self.add_attribute('R_refine', [60.0e3, 40e3, 20e3])

        # Storm / met-forcing selection.
        #
        # Preferred API: storm_family ("parametric" | "gridded" | "none") plus
        # storm_subtype (a model name for parametric, "gridded" for gridded).
        # These default to None; when unset, write() falls back to the legacy
        # storm_specification_type (string or integer), which is retained as a
        # plain stored attribute for backwards compatibility.  See the
        # forcing_subtype_registry above and met_forcing_refactor.md Section 7.
        self.add_attribute('storm_family', None)
        self.add_attribute('storm_subtype', None)
        self.add_attribute('storm_type', None)  # Backwards compatibility
        self.add_attribute('storm_specification_type', 0) # Legacy selection
        self.add_attribute("storm_file", None) # File containing data

    def read(self, path: Path=Path("surge.data"), force: bool=False):
        """Read surge data file"""

        with Path(path).open() as data_file:
            # Header
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()
            data_file.readline()

            self.wind_forcing = bool(data_file.readline())
            self.drag_law = int(data_file.readline().split("=:")[0])
            self.pressure_forcing = bool(data_file.readline().split("=:")[0])
            self.rotation_override = data_file.readline().split("=:")[0]
            data_file.readline()

            self.wind_index = int(data_file.readline().split("=:")[0]) - 1
            self.pressure_index = int(data_file.readline().split("=:")[0]) - 1
            self.display_landfall_time = bool(data_file.readline().split("=:")[0])
            self.storm_time_scale = float(data_file.readline().split("=:")[0])
            self.t_ramp_on = float(data_file.readline().split("=:")[0])
            self.t_ramp_off = float(data_file.readline().split("=:")[0])
            data_file.readline()
            data_file.readline() # TODO: Extra empty line, should fix

            # AMR parameters
            self.wind_refine = self._parse_value(data_file.readline().split("=:")[0])
            self.R_refine = self._parse_value(data_file.readline().split("=:")[0])
            # data_file.readline()

            # Storm / met-forcing selection (family + subtype tokens).
            self.storm_family = data_file.readline().split("=:")[0].strip()[1:-1]
            self.storm_subtype = data_file.readline().split("=:")[0].strip()[1:-1]
            # Populate the legacy integer selector for backwards compatibility.
            _, _, self.storm_specification_type = \
                resolve_forcing_subtype(self.storm_subtype)
            line = data_file.readline().split("=:")[0]
            if line.strip()[0] == "'":
                self.storm_file = line.strip()[1:-1]
            else:
                raise IOError("Error reading storm file name.")


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
        self.data_write("storm_time_scale",
                        description='(Model storm time scale: >1 slower, <1 faster)')
        self.data_write("t_ramp_on",
                        description='(Temporal onset ramp width in seconds)')
        self.data_write("t_ramp_off",
                        description='(Temporal cutoff ramp width in seconds)')
        self.data_write()

        # AMR storm refinement criteria
        # Handle older style of refinement for turning it off, F -> None
        if isinstance(self.wind_refine, bool):
            if not self.wind_refine:
                self.wind_refine = None
        if isinstance(self.R_refine, bool):
            if not self.R_refine:
                self.R_refine = None
        self.data_write()
        
        if isinstance(self.wind_refine, float):
            self.wind_refine = [self.wind_refine]
        self.data_write('wind_refine', description='(Wind speed refinement criteria)')
        if isinstance(self.R_refine, float):
            self.R_refine = [self.R_refine]
        self.data_write('R_refine', description='(Wind speed refinement criteria)')

        # Storm / met-forcing selection.
        #
        # Resolve to a canonical (family, subtype) pair.  The preferred inputs
        # are storm_family/storm_subtype; when unset we fall back to the legacy
        # storm_specification_type (or the deprecated storm_type alias), which
        # may be a name, an alias, or an integer code.
        if self.storm_family is not None or self.storm_subtype is not None:
            # storm_subtype is authoritative; family is validated for
            # consistency if the caller also supplied it.
            if self.storm_subtype is None:
                raise ValueError("storm_family set without storm_subtype.")
            family, subtype, _ = resolve_forcing_subtype(self.storm_subtype)
            if (self.storm_family is not None
                    and self.storm_family.lower() != family):
                raise ValueError(
                    f"storm_family '{self.storm_family}' is inconsistent with "
                    f"storm_subtype '{self.storm_subtype}' (family {family}).")
        else:
            legacy = self.storm_type if self.storm_type is not None \
                else self.storm_specification_type
            family, subtype, _ = resolve_forcing_subtype(legacy)

        # Write the explicit family/subtype tokens (quoted single words, read
        # list-directed by the Fortran side; see met_forcing_module.f90).
        self.data_write(name="storm_family", value=family,
                        description="(Forcing family: parametric | gridded | none)")
        self.data_write(name="storm_subtype", value=subtype,
                        description="(Forcing subtype / parametric model)")
        self.data_write(name="storm_file", description='(Path to storm data)')

        self.close_data_file()


#: Canonical alias for :class:`SurgeData` under the meteorological-forcing name.
#: ``clawpack.geoclaw.met`` re-exports this as ``MetData``.  The on-disk file is
#: still written as ``surge.data`` (unchanged wire contract with the Fortran).
MetData = SurgeData


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
                lower = self._parse_value(data_file.readline())
                upper = self._parse_value(data_file.readline())
                depths = self._parse_value(data_file.readline())
                coeff = self._parse_value(data_file.readline())
                self.friction_regions.append([lower, upper, depths, coeff])
                data_file.readline()
            self.friction_files = [] # Is not supported


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
