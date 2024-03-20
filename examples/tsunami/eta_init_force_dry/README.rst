
.. _geoclaw_examples_eta_init_force_dry:

Tsunami modeling example with subsidence and dry regions below sea level
========================================================================

This illustrates two features of GeoClaw available starting in Version 5.7.0:

- The ability to use a `set_eta_init` subroutine to set a spacially varying
  initial surface elevation (rather than constant `sea_level` everywhere).
  For documentation, see `Set eta init 
  <https://www.clawpack.org/set_eta_init.html>`__

- The ability to identify some regions as initially dry even though the
  topography value is below the initial eta elevation specified.
  This can be used to initialize dry land behind dikes, for example.
  For documentation, see `Force dry
  <https://www.clawpack.org/force_dry.html>`__

See the Jupyter notebooks (or rendered html versions linked below)
for more discussion of
these features and how they are used in this example.

To run this example
--------------------

Input data can be created via::

    make input

This executes the Python script `make_input_files.py`, which was created from
the Jupyter notebook `make_input_files.ipynb`.  

The input files should appear in a new subdirectory `input_files`.

After creating the input files, GeoClaw can be run and plots created as usual
via::

    make .plots

using the settings in `setrun.py`.  

In addition, the notebook `run_geoclaw.ipynb` can be used to run a sequence of
experiments and display the results, with some comments on what is
illustrated in this example.

**Rendered versions of the Jupyter notebooks** for this example can 
be viewed from the `Clawpack gallery version of this file.
<http://www.clawpack.org/gallery/_static/geoclaw/examples/tsunami/eta_init_force_dry/README.html>`__


Version
-------

- These capabilities and this example were first introduced in Clawpack v5.7.0

- The results changed slightly in v5.8.0 due to a previous bug being
  fixed (the `maxlevel` parameter in flagregions specified as ruled
  rectangles was not properly handled and so the refinement to level
  2 did not occur as expected along the coast at early times).
