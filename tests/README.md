This directory is for unit tests for the Python portions of GeoClaw.  You can
run all of these by running `pytest tests` from `$CLAW/geoclaw` or individually
by 
```
pytest test_[TEST_NAME].py
```
You can also run all of the python tests by using the marker `python` with
```
pytest -m python
```
again in `$CLAW/geoclaw`.

## Subdirectories

### `netcdf/`

Unit tests for the NetCDF inspector and descriptor utilities in
`src/python/geoclaw/netcdf_utils.py`.  The classes under test are:

- `NetCDFInspector` — base class; coordinate discovery and convention detection
- `TopoInspector` — bathymetry-specific validation (units, NaN-in-crop checking)
- `MetInspector` — met-forcing validation (variable consistency, unit contracts,
  time-offset computation)
- `CFNormalizer` — standalone utility to fix CF metadata in-place
- `DescriptorWriter` — serialises inspector metadata to Fortran-readable
  descriptor files

Run just this suite with:
```
pytest tests/netcdf/
```

See more details on running and updating tests or debugging issues in the
[Clawpack documentation](https://www.clawpack.org/testing.html).
