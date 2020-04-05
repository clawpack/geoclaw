
.. _geoclaw_examples_eta_init_force_dry:

Tsunami modeling example with subsidence and dry regions below sea level
========================================================================

This illustrates two features of GeoClaw available starting in Version 5.7.0:

- The ability to use a `set_eta_init` subroutine to set a spacially varying
  initial surface elevation (rather than constant `sea_level` everywhere).
  For documentation, see 
  http://depts.washington.edu/clawpack/sampledocs/dev_v5.7.0/set_eta_init.html

- The ability to identify some regions as initially dry even though the
  topography value is below the initial eta elevation specified.
  This can be used to initialize dry land behind dikes, for example.
  For documentation, see 
  http://depts.washington.edu/clawpack/sampledocs/dev_v5.7.0/force_dry.html

See the Jupyter notebooks (or rendered html versions) for more discussion of
these features and how they are used in this example.

To run this example
--------------------

Input data can be created via::

    make input

This executes the Python script `make_input_files.py`, which was created from
the Jupyter notebook `make_input_files.ipynb`.  A rendered version of the
notebook is also stored as make_input_files.html, which explains the process
used to generate the input data.

The input files should appear in a new subdirectory `input_files`.

After creating the input files, GeoClaw can be run and plots created as usual
via::

    make .plots

using the settings in `setrun.py`.  

In addition, the notebook `run_geoclaw.ipynb` can be used to run a sequence of
experiments and display the results, with some comments on what is
illustrated in this example.   A rendered version of this notebook is
available as `run_geoclaw.html`.


Version
-------

- New in Version 5.7.0.
