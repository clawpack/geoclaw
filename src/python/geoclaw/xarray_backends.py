r"""
xarray backends module: $CLAW/geoclaw/src/python/geoclaw/xarray_backends.py

Xarray backends for GeoClaw fixed grids and fgmax grids.

These only work for point_style = 2 (uniform regular) and have a dependency on
xarray and rioxarray. Xarray provides the core datastructure and rioxarray
provides an interface to rasterio, used to assign geospatial projection
information.

- https://docs.xarray.dev/en/stable/index.html
- https://corteva.github.io/rioxarray/stable/readme.html
- https://rasterio.readthedocs.io/en/latest/index.html

The expectation is that you run commands that use these backends from the same
directory you run "make output".

Includes:

- class FGMaxBackend: Xarray backend for fgmax grids.
- class FGOutBackend: Xarray backend for fgout grids.

Usage:

.. code-block:: python
    import rioxarray  # activate the rio accessor
    import xarray as xr
    from clawpack.geoclaw.xarray_backends import FGOutBackend, FGMaxBackend

    # An example of a fgout grid.

    filename = '_output/fgout0001.b0001
    # provide the .bxxx file if binary format is used or the
    # .qxxx file if ascii format is used.
    # the format, fg number, and frame number are inferred from the filename.

    ds = xr.open_dataset(
        filename,
        engine=FGOutBackend,
        backend_kwargs={
            'qmap': 'geoclaw',
            'epsg':epsg_code
                })
    # ds is now an xarray object. It can be interacted with directly or written to netcdf using
    ds.write_netcdf('filename.nc')

    # A 'qmap' backend_kwargs is required as it indicates the qmap used for
    # selecting elements of q to include in fgout. See
    # https://www.clawpack.org/dev/fgout.html#specifying-q-out-vars

    # Optionally, provide an epsg code to assign the associated coordinate system to the file.
    # default behavior assigns no coordinate system.

    # An example of a fgmax grid.
    filename = "_output/fgmax0001.txt"
    ds = xr.open_dataset(filename, engine=FGMaxBackend, backend_kwargs={'epsg':epsg_code})


Dimensions:

Files opened with FGOutBackend will have dimensions (time, y, x).
Files opened with FGMaxBackend will have dimensions (y, x).

Variable naming:

For fixed grid files, the dataset will have the elements of q specified
by qmap and q_out_vars. This may be as many as the full list described below or
as few as specified by the user.

For fixed grid geoclaw files, the full set of variables is:
- h
- hu
- hv
- eta
- B

For fixed grid geoclaw-bouss files, the full set of variables is:
- h
- hu
- hv
- huc
- hvc
- eta
- B

Fixed grid dclaw files, the full set of variables is:
- h
- hu
- hv
- hm
- pb
- hchi
- bdif
- eta
- B

Depending on the number of variables specified in the setrun.py fgmax files will
have a portion of the following variables:

If rundata.fgmax_data.num_fgmax_val == 1

- arrival_time, Wave arrival time (based on eta>sea_level + fg.arrival_tol)
- h_max, Maximum water depth
- eta_max, Maximum water surface elevation
- h_max_time, Time of maximum water depth
- B, Basal topography at the first time fgmax first monitored maximum amr level
- level, Maximum amr level

If rundata.fgmax_data.num_fgmax_val == 2:

- s_max, Maximum velocity
- s_max_time, Time of maximum velocity

If rundata.fgmax_data.num_fgmax_val == 5:

- hs_max, Maximum momentum
- hs_max_time, Time of maximum momentum
- hss_max, Maximum momentum flux
- hss_max_time, Time of maximum momentum flux
- h_min, Minimum depth
- h_min_time, Time of minimum depth


See the following links for additional information about xarray Backends.

- https://docs.xarray.dev/en/stable/generated/xarray.backends.BackendEntrypoint.html#xarray.backends.BackendEntrypoint
- https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html
"""

import os

import numpy as np

try:
    import rioxarray  # activate the rio accessor
    import xarray as xr
except ImportError:
    raise ImportError(
        "rioxarray and xarray are required to use the FGOutBackend and FGMaxBackend"
    )

from clawpack.geoclaw import fgmax_tools, fgout_tools
from xarray.backends import BackendEntrypoint

_qunits = {
    "h": "meters",
    "hu": "meters squared per second",
    "hv": "meters squared per second",
    "hm": "meters",
    "pb": "newton per meter squared",
    "hchi": "meters",
    "bdif": "meters",
    "eta": "meters",
    "B": "meters",
    "huc": "varies",
    "hvc": "varies",
}

time_unit = "seconds"
reference_time = "model start"
space_unit = "meters"
nodata = np.nan


def _prepare_var(data, mask, dims, units, nodata, long_name):
    "Reorient array into the format xarray expects and generate the mapping."
    if mask is not None:
        data[mask] = nodata
    data = data.T
    data = np.flipud(data)
    mapping = (
        dims,
        data,
        {"units": units, "_FillValue": nodata, "long_name": long_name},
    )
    return mapping


class FGOutBackend(BackendEntrypoint):
    "Xarray Backend for Clawpack fixed grid format."

    def open_dataset(
        self,
        filename,  # path to fgout file.
        qmap="geoclaw",  # qmap value for FGoutGrid ('geoclaw', 'dclaw', or 'geoclaw-bouss')
        epsg=None,  # epsg code
        dry_tolerance=0.001,
        # dry tolerance used for for masking elements of q that are not
        # eta or B. Default behavior is to mask all elements of q except for
        # eta and B where h>0.001.
        # used only if h or eta-B is available based on q_out_vars.
        # if dry_tolerance = None, no masking is applied to any variable.
        drop_variables=None,  # name of any elements of q to drop.
    ):

        if drop_variables is None:
            drop_variables = []

        full_path = os.path.abspath(filename)
        filename = os.path.basename(full_path)
        outdir = os.path.basename(os.path.dirname(full_path))

        # filename has the format fgoutXXXX.qYYYY (ascii)
        # or fgoutXXXX.bYYYY (binary)
        # where XXXX is the fixed grid number and YYYY is the frame
        # number.
        type_code = filename.split(".")[-1][0]
        fgno = int(filename.split(".")[0][-4:])
        frameno = int(filename.split(".")[-1][-4:])
        if type_code == "q":  # TODO, is this correct?
            output_format = "ascii"
        elif type_code == "b":
            output_format = "binary32"  # format of fgout grid output
        else:
            raise ValueError("Invalid FGout output format. Must be ascii or binary.")

        fgout_grid = fgout_tools.FGoutGrid(
            fgno=fgno, outdir=outdir, output_format=output_format, qmap=qmap
        )

        if fgout_grid.point_style != 2:
            raise ValueError("FGOutBackend only works with fg.point_style=2")

        try:
            fgout_grid.read_fgout_grids_data()
        except AssertionError:
            fgout_grid.read_fgout_grids_data_pre511()

        fgout = fgout_grid.read_frame(frameno)

        time = fgout.t
        # both come in ascending. flip to give expected order.
        x = fgout.x
        y = np.flipud(fgout.y)
        nj = len(x)
        ni = len(y)

        # mask based on dry tolerance
        if dry_tolerance is not None:
            try:
                mask = fgout.h < dry_tolerance
                # internally fgout_tools will try fgou.eta-fgout.B if h is not present.
            except AttributeError:
                print("FGOutBackend: No h, eta, or B. No mask applied.")
                mask = np.ones((nj, ni), dtype=bool)

        # create data_vars dictionary
        data_vars = {}
        for i, i_var in enumerate(fgout_grid.q_out_vars):

            # construct variable

            # Find the varname in fgout.qmap associated with
            # q_out_vars[i]
            varname = None
            for name, index in fgout_grid.qmap.items():
                if i_var == index:
                    varname = name

            # construct xarray dataset if varname not in drop vars.
            if varname not in drop_variables:

                Q = fgout.q[i, :, :]

                # mask all but eta based on h presence.
                if (varname not in ("B", "eta")) and dry_tolerance is not None:
                    Q[mask] = nodata

                # to keep xarray happy, need to transpose and flip ud.
                Q = Q.T
                Q = np.flipud(Q)
                Q = Q.reshape((1, ni, nj))  # reshape to add a time dimension.

                data_array_attrs = {"units": _qunits[varname], "_FillValue": nodata}

                data_vars[varname] = (
                    [
                        "time",
                        "y",
                        "x",
                    ],
                    Q,
                    data_array_attrs,
                )

        ds_attrs = {"description": "Clawpack model output"}

        ds = xr.Dataset(
            data_vars=data_vars,
            coords=dict(
                x=(["x"], x, {"units": space_unit}),
                y=(["y"], y, {"units": space_unit}),
                time=("time", [time], {"units": "seconds"}),
                reference_time=reference_time,
            ),
            attrs=ds_attrs,
        )

        if epsg is not None:
            ds.rio.write_crs(
                epsg,
                inplace=True,
            ).rio.set_spatial_dims(
                x_dim="x",
                y_dim="y",
                inplace=True,
            ).rio.write_coordinate_system(inplace=True)
            # https://corteva.github.io/rioxarray/stable/getting_started/crs_management.html#Spatial-dimensions
            # https://gis.stackexchange.com/questions/470207/how-to-write-crs-info-to-netcdf-in-a-way-qgis-can-read-python-xarray

        return ds

    open_dataset_parameters = ["filename", "drop_variables"]

    description = "Use Clawpack fixed grid output files in Xarray"
    url = "https://www.clawpack.org/fgout.html"


class FGMaxBackend(BackendEntrypoint):
    "Xarray Backend for Clawpack fgmax grid format."

    def open_dataset(
        self,
        filename,
        epsg=None,
        drop_variables=None,
        clip=False,  # if True, clip the entire array to the extent where arrival_time is defined.
    ):

        if drop_variables is None:
            drop_variables = []

        # expectation is for standard clawpack organization,
        # e.g., output and 'fgmax_grids.data' in _output
        fgno = int(os.path.basename(filename).split(".")[0][-4:])
        outdir = os.path.dirname(filename)
        data_file = os.path.join(outdir, "fgmax_grids.data")

        fg = fgmax_tools.FGmaxGrid()
        fg.read_fgmax_grids_data(fgno=fgno, data_file=data_file)
        if fg.point_style != 2:
            raise ValueError("FGMaxBackend only works with fg.point_style=2")

        fg.read_output(outdir=outdir)

        # Construct the x and y coordinates
        # Both come in ascending, therefore flip y so that it is ordered as expected.
        x = fg.x
        y = np.flipud(fg.y)

        # Construct the data_vars array. To organize the
        # data in the way expected by netcdf standards, need
        # to both transpose and flipud the array.

        data_vars = {}

        data_vars["arrival_time"] = _prepare_var(
            fg.arrival_time.data,
            fg.arrival_time.mask,
            [
                "y",
                "x",
            ],
            "seconds",
            nodata,
            "Wave arrival time",
        )

        data_vars["h_max"] = _prepare_var(
            fg.h.data,
            fg.h.mask,
            [
                "y",
                "x",
            ],
            "meters",
            nodata,
            "Maximum water depth",
        )

        data_vars["eta_max"] = _prepare_var(
            fg.h.data + fg.B.data,
            fg.h.mask,
            [
                "y",
                "x",
            ],
            "meters",
            nodata,
            "Maximum water surface elevation",
        )

        data_vars["h_max_time"] = _prepare_var(
            fg.h_time.data,
            fg.h_time.mask,
            [
                "y",
                "x",
            ],
            "seconds",
            nodata,
            "Time of maximum water depth",
        )

        data_vars["B"] = _prepare_var(
            fg.B,
            None,
            [
                "y",
                "x",
            ],
            "meters",
            nodata,
            "Basal topography at the first time fgmax first monitored maximum amr level",
        )

        data_vars["level"] = _prepare_var(
            fg.level,
            None,
            [
                "y",
                "x",
            ],
            "(no units)",
            -1,
            "Maximum amr level",
        )

        if hasattr(fg, "s"):
            if hasattr(fg.s, "data"):
                data_vars["s_max"] = _prepare_var(
                    fg.s.data,
                    fg.s.mask,
                    [
                        "y",
                        "x",
                    ],
                    "meters per second",
                    nodata,
                    "Maximum velocity",
                )
                data_vars["s_max_time"] = _prepare_var(
                    fg.s_time.data,
                    fg.s_time.mask,
                    [
                        "y",
                        "x",
                    ],
                    "seconds",
                    nodata,
                    "Time of maximum velocity",
                )
        if hasattr(fg, "hs"):
            if hasattr(fg.hs, "data"):
                data_vars["hs_max"] = _prepare_var(
                    fg.hs.data,
                    fg.hs.mask,
                    [
                        "y",
                        "x",
                    ],
                    "meters squared per second",
                    nodata,
                    "Maximum momentum",
                )
                data_vars["hs_max_time"] = _prepare_var(
                    fg.hs_time.data,
                    fg.hs_time.mask,
                    [
                        "y",
                        "x",
                    ],
                    "seconds",
                    nodata,
                    "Time of maximum momentum",
                )

                data_vars["hss_max"] = _prepare_var(
                    fg.hss.data,
                    fg.hss.mask,
                    [
                        "y",
                        "x",
                    ],
                    "meters cubed per second squared",
                    nodata,
                    "Maximum momentum flux",
                )
                data_vars["hss_max_time"] = _prepare_var(
                    fg.hss_time.data,
                    fg.hss_time.mask,
                    [
                        "y",
                        "x",
                    ],
                    "seconds",
                    nodata,
                    "Time of maximum momentum flux",
                )

                data_vars["h_min"] = _prepare_var(
                    fg.hmin.data,
                    fg.hmin.mask,
                    [
                        "y",
                        "x",
                    ],
                    "meters",
                    nodata,
                    "Minimum depth",
                )
                data_vars["h_min_time"] = _prepare_var(
                    fg.hmin_time.data,
                    fg.hmin_time.mask,
                    [
                        "y",
                        "x",
                    ],
                    "seconds",
                    nodata,
                    "Time of minimum depth",
                )

        # drop requested variables
        for var in drop_variables:
            if var in data_vars.keys():
                del data_vars[var]

        # Construct the values from
        ds_attrs = {"description": "D-Claw model output"}

        ds = xr.Dataset(
            data_vars=data_vars,
            coords=dict(
                x=(["x"], x, {"units": "meters"}),
                y=(["y"], y, {"units": "meters"}),
            ),
            attrs=ds_attrs,
        )

        if epsg is not None:

            ds.rio.write_crs(
                epsg,
                inplace=True,
            ).rio.set_spatial_dims(
                x_dim="x",
                y_dim="y",
                inplace=True,
            ).rio.write_coordinate_system(inplace=True)
            # https://corteva.github.io/rioxarray/stable/getting_started/crs_management.html#Spatial-dimensions
            # https://gis.stackexchange.com/questions/470207/how-to-write-crs-info-to-netcdf-in-a-way-qgis-can-read-python-xarray

        # clip
        if clip:
            clip_data = fg.arrival_time.data >= 0
            clip_data = clip_data.T
            clip_data = np.flipud(clip_data)
            ds = _clip(ds, clip_data)

        return ds

    open_dataset_parameters = ["filename", "drop_variables", "clip"]

    description = "Use Clawpack fix grid monitoring files in Xarray"
    url = "https://www.clawpack.org/fgmax.html#fgmax"


def _clip(ds, clip_data):
    # clip DS to the extent where clip_data is true.
    # useful for when there are large areas where nothing happened.

    nrow, ncol = clip_data.shape

    valid_col = np.sum(clip_data, axis=0) > 0
    valid_row = np.sum(clip_data, axis=1) > 0

    sel_rows = np.nonzero(valid_row)[0]
    sel_cols = np.nonzero(valid_col)[0]

    first_row = sel_rows[0]
    last_row = sel_rows[-1]
    first_col = sel_cols[0]
    last_col = sel_cols[-1]

    ds = ds.isel(y=np.arange(first_row, last_row), x=np.arange(first_col, last_col))
    return ds
