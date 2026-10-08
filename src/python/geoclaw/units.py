#!/usr/bin/env python
# encoding: utf-8
r"""
Defines common units used throughout GeoClaw along with functions for
converting between them.
"""

from __future__ import print_function
from __future__ import absolute_import

import sys
import collections
import dataclasses

from clawpack.geoclaw.data import LAT2METER

# Define baseline units so that we can convert to these and then back to
# requested.
standard_units = {'time': 's',
                  'spherical_length': 'lat-long',
                  'length': 'm',
                  'speed': 'm/s',
                  'radius': 'm',
                  'pressure': 'Pa',
                  'temperature': 'C'}

# Unit conversion definitions - handles conversions to standard units above
# The dictionary contains a list of values, the first is the conversion factor
# and the second a human readable version of the unit.
conversion_func = {}
units = {}

# Length
conversion_func['m'] = [lambda L: L,
                        lambda L: L]
conversion_func['cm'] = [lambda L: L * 1e-2,
                         lambda L: L / 1e-2]
conversion_func['km'] = [lambda L: L * 1e3,
                         lambda L: L / 1e3]
conversion_func['miles'] = [lambda L: L * 1.60934e3,
                            lambda L: L / 1.60934e3]
conversion_func['nmi'] = [lambda L: L * 1852.0,
                          lambda L: L / 1852.0]
conversion_func['lat-long'] = [lambda L: L * LAT2METER,
                               lambda L: L / LAT2METER]
units['length'] = collections.OrderedDict({'m': 'meters', 'cm': 'centimeters', 
                                           'km': 'kilometers', 'miles': 'miles',
                                           'nmi': 'nautical miles', 
                                           'lat-long': 'longitude-latitude'})

# Pressure - Rigidity
conversion_func['Pa'] = [lambda P: P, lambda P: P]
conversion_func['hPa'] = [lambda P: P * 1e2, 
                          lambda P: P / 1e2]
conversion_func['KPa'] = [lambda P: P * 1e3,
                          lambda P: P / 1e3]
conversion_func['MPa'] = [lambda P: P * 1e6,
                          lambda P: P / 1e6]
conversion_func['GPa'] = [lambda P: P * 1.e9,
                          lambda P: P / 1.e9]
conversion_func['mbar'] = [lambda P: P * 1e2,
                           lambda P: P / 1e2]
conversion_func['dyne/cm^2'] = [lambda P: P * 0.1,
                                lambda P: P / 0.1]
conversion_func['dyne/m^2'] = [lambda P: P * 1.e-5,
                               lambda P: P / 1.e-5]
units['pressure'] = {'Pa': 'pascals', 'hPa': 'hectopascals', 
                     'KPa': 'kilopascals', 'MPa': 'megapascals', 
                     'GPa': 'gigapascals', 'mbar': 'millibar', 
                     'dyne/cm^2': 'Dynes/cm^2', 'dyne/m^2': 'Dynes/m^2'}

# Speeds
conversion_func['km/h'] = [lambda v: v * 1e3 / (60**2), 
                           lambda v: v * 60**2 / 1e3]
conversion_func['mph'] = [lambda v: v * 0.44704, lambda v: v / 0.44704]
conversion_func['m/s'] = [lambda v: v, lambda v: v]
conversion_func['knots'] = [lambda v: v * 0.51444444, 
                            lambda v: v / 0.51444444]
units['speed'] = {'m/s': 'meters/second', 'knots': 'knots (nm / hour)', 
                  'km/h': 'kilometers/hour', 'mph': "miles / hour"}

# Moments
conversion_func['N-m'] = [lambda M: M, lambda M: M]
conversion_func['dyne-cm'] = [lambda M: M * 1.e-7,
                              lambda M: M / 1.e-7]
units['moment'] = {'N-m': "Newton-Meters", 'dyne-cm': "Dynes-Centimeter"}

# Temperature
conversion_func['C'] = [lambda temp: temp, lambda temp: temp]
conversion_func['F'] = [lambda temp: (temp - 32.0) * 5.0 / 9.0, 
                        lambda temp: temp * 9.0 / 5.0 + 32.0]
conversion_func['K'] = [lambda temp: temp + 273.15,
                        lambda temp: temp - 273.15]
units['temperature'] = {'C': "Celsius", 'F': "Fahrenheit", 'K': "Kelvin"}

# Time
conversion_func['s'] = [lambda t: t, lambda t: t]
conversion_func['min'] = [lambda t: t * 60.0, lambda t: t / 60.0]
conversion_func['minutes'] = conversion_func['min']
conversion_func['h'] = [lambda t: t * 3600.0, lambda t: t / 3600.0]
conversion_func['hours'] = conversion_func['h']
conversion_func['days'] = [lambda t: t * 86400.0, lambda t: t / 86400.0]
units['time'] = collections.OrderedDict({'s': 'seconds', 'min': 'minutes',
                                         'h': 'hours', 'days': 'days'})

# Unit contract for GeoClaw NetCDF input.
# These are the units that Fortran assumes when reading NetCDF data.
# Python interrogators must verify and, if necessary, record a conversion
# factor before writing the descriptor so that Fortran receives data in
# these units.
GEOCLAW_NETCDF_UNITS = {
    'topo':     'm',      # bathymetry / topography elevation
    'wind_u':   'm/s',    # eastward wind component
    'wind_v':   'm/s',    # northward wind component
    'pressure': 'Pa',     # surface pressure (absolute)
    'time':     's',      # time coordinate (seconds from descriptor offset)
    # TODO: 'precipitation': ???  (unit TBD)
    # TODO: 'friction':      ???  (unit TBD)
}


# ---------------------------------------------------------------------------
# Units policy registry
# ---------------------------------------------------------------------------
# GeoClaw's units policy is stated in dev/design/units_policy.md.  This registry
# is the machine-readable form of the table in that document, and it is the
# single source of truth for both:
#
#   * tests/test_units_policy.py, which drives each real reader with a fixture
#     and asserts the behavior declared here -- so a row cannot claim
#     something the code does not do; and
#   * the generated table in dev/design/units_policy.md, rendered from these
#     rows -- so the document cannot drift from the registry.
#
# A row whose `gap` is not None does *not* yet meet the policy.  Its
# conformance test is marked xfail(strict=True), so the row is visible as a
# known hole and the marker cannot be left behind once the hole is closed.


@dataclasses.dataclass(frozen=True)
class UnitsPolicyRow:
    """One input path's declared units behavior.

    The field values are the vocabulary the conformance test understands; see
    ON_MISSING_VALUES etc. below.
    """

    key: str                   # stable identifier used by the tests
    reader: str                # human-readable entry point
    contract: str              # GeoClaw's internal unit for this quantity
    declared_in_file: bool     # can the format carry a units declaration?
    override: str              # the argument that states units explicitly
    on_missing: str            # no declaration present
    on_convertible: str        # recognized unit that is not the contract unit
    on_unrecognized: str       # unit string we cannot interpret
    magnitude_check: bool      # is the post-conversion sanity check applied?
    gap: str = ''              # non-empty => does not yet conform, and why


# Vocabulary for the behavior fields, so a typo in a row is caught rather
# than silently producing an untested case.
ON_MISSING_VALUES = frozenset({
    'raise',            # refuse to guess (the policy default)
    'warn+assume',      # assume the contract unit, but say so
    'silent-assume',    # assume the contract unit with no message (a hole)
    'format-default',   # the file format documents the unit; use it
    'n/a',              # the format has no notion of declared units
})
ON_CONVERTIBLE_VALUES = frozenset({'convert+warn', 'convert', 'n/a'})
ON_UNRECOGNIZED_VALUES = frozenset({'raise', 'n/a'})


UNITS_POLICY: tuple[UnitsPolicyRow, ...] = (
    UnitsPolicyRow(
        key='topo_netcdf',
        reader='Topography.read (topo_type=4)',
        contract='m',
        declared_in_file=True,
        override="nc_params={'assume_units': str}",
        on_missing='raise',
        on_convertible='convert+warn',
        on_unrecognized='raise',
        magnitude_check=True,
    ),
    UnitsPolicyRow(
        key='topo_ascii',
        reader='Topography.read (topo_type=1,2,3)',
        contract='m',
        declared_in_file=False,
        override='none',
        on_missing='silent-assume',
        on_convertible='n/a',
        on_unrecognized='n/a',
        magnitude_check=False,
        gap='ASCII carries no units and has no override; elevation in cm or '
            'feet is read as meters with no message and no sanity check.',
    ),
    UnitsPolicyRow(
        key='dtopo_netcdf_dz',
        reader='DTopoInspector (deformation)',
        contract='m',
        declared_in_file=True,
        override='assume_units (str)',
        on_missing='raise',
        on_convertible='convert+warn',
        on_unrecognized='raise',
        magnitude_check=False,
    ),
    UnitsPolicyRow(
        key='dtopo_netcdf_time',
        reader='DTopoInspector (time axis)',
        contract='s',
        declared_in_file=True,
        override='none',
        on_missing='warn+assume',
        on_convertible='convert',
        on_unrecognized='raise',
        magnitude_check=False,
    ),
    UnitsPolicyRow(
        key='dtopo_ascii',
        reader='DTopography.read (dtopo_type=1,2,3)',
        contract='m',
        declared_in_file=False,
        override='none',
        on_missing='silent-assume',
        on_convertible='n/a',
        on_unrecognized='n/a',
        magnitude_check=False,
        gap='ASCII dtopo carries no units and has no override.',
    ),
    UnitsPolicyRow(
        key='met_netcdf',
        reader='MetInspector (wind, pressure)',
        contract='m/s, Pa',
        declared_in_file=True,
        override='assume_units (bool), format_units (dict)',
        on_missing='raise',
        on_convertible='convert+warn',
        on_unrecognized='raise',
        magnitude_check=True,
    ),
    UnitsPolicyRow(
        key='subfault_generic',
        reader='Fault.read (columnar subfault files)',
        contract='m, Pa',
        declared_in_file=False,
        override='input_units (dict)',
        on_missing='warn+assume',
        on_convertible='convert',
        on_unrecognized='raise',
        magnitude_check=False,
    ),
    UnitsPolicyRow(
        key='subfault_csv',
        reader='CSVFault.read (units in column headings)',
        contract='m, Pa',
        declared_in_file=True,
        override='input_units (dict), overrides the heading',
        on_missing='warn+assume',
        on_convertible='convert',
        on_unrecognized='raise',
        magnitude_check=False,
    ),
)


def render_units_policy_table() -> str:
    """Render :data:`UNITS_POLICY` as the Markdown table used in the design doc.

    ``dev/design/units_policy.md`` holds this between generated-block markers
    and a test asserts the two agree, so the prose cannot drift from the code.
    """
    header = ('| Path | Contract | Declared in file | Override | Missing | '
              'Non-contract | Unrecognized | Magnitude | Conforms |')
    sep = '|' + '---|' * 9
    lines = [header, sep]
    for row in UNITS_POLICY:
        conforms = 'yes' if not row.gap else '**no** -- ' + row.gap
        lines.append(
            f"| `{row.reader}` | {row.contract} | "
            f"{'yes' if row.declared_in_file else 'no'} | {row.override} | "
            f"{row.on_missing} | {row.on_convertible} | {row.on_unrecognized} "
            f"| {'yes' if row.magnitude_check else 'no'} | {conforms} |")
    return "\n".join(lines)


def units_available():
    r"""
    Constructs a string suitable for reading detailing the units available.
    """
    output = ""
    for (measurement_type, measurement_units) in units.items():
        output = "\n".join((output, measurement_type.capitalize()))
        for (abbrv, full_name) in measurement_units.items():
            output = "\n".join((output, "  %s (%s)" % (abbrv, full_name)))

    return output


# Unit conversion function
def convert(value, old_units, new_units, verbose=False):
    r"""Convert *value* from *old_units* to *new_units*

    :Note:
    Currently this function only handles multiplicative conversions.  The
    reasoning behind not just returning this conversion factor as this function
    in the future will also support more complex unit conversions, e.g.
    converting between temperature scales.

    :Input:
     - *value* (ndarray or float) The value(s) to be converted.
     - *old_units* (string) Type of units that value comes in as.
     - *new_units* (string) Type of units that value should be converted to.
     - *verbose* (bool) Verbose output (default is False)

    :Output:
     - (ndarray or float) The converted value(s)
    """

    found_type = None
    for (measurement_type, measurement_units) in units.items():
        if old_units in measurement_units:
            found_type = measurement_type
            break

    if found_type is None:
        raise ValueError("Units %s not found in list of " % str(old_units),
                         "supported conversions." )

    if new_units not in units[found_type].keys():
        raise ValueError("Units %s not found in list of " % str(new_units),
                         "supported conversions of %s type." % found_type)

    if verbose:
        print("Convert %s %s to %s." % (value, old_units, new_units))

    return conversion_func[new_units][1](conversion_func[old_units][0](value))


if __name__ == '__main__':
    # Add commandline unit conversion capability
    if len(sys.argv) == 4:
        convert(float(sys.argv[1]), sys.argv[2], sys.argv[3])
    else:
        # Usage and available units
        print("Usage:  Convert value in units to new units.")
        print("  units <value> <old units> <new units>")
        print("where <old units> and <new units> are one of the available")
        print("units listed below.")
        print("")
        print("Available Units:")
        print("  First value is the abbreviation that should be used as an")
        print("  input unit while the second is the full name of the unit.")
        print(units_available())
