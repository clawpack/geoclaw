
.. _geoclaw_examples_tsunami_radial-ocean-island-fgmax:

Radial ocean with an island and fgmax grids
============================================

The topography is a synthetic radially symmetric ocean with a continental
shelf.  An island is added at one location on the shelf.  

The topography was adapted from a test problem in this paper:

- The GeoClaw software for depth-averaged flows with adaptive refinement
  by M. J. Berger, D. L. George, R. J. LeVeque, and K. T. Mandli,
  Advances in Water Resources 34 (2011), pp. 1195-1206.
  `<http://faculty.washington.edu/rjl/pubs/awr10/index.html>`_

To create the input data::

    make data
    make input

Running `make input` executes `make_input_files.py`, or you can instead run
the notebook `make_input_files.ipynb`.
For sample output see `make_input_files.html <make_input_files.html>`_

After running the code, ::

    make fgmax_plots

creates the fgmax plots by running `process_fgmax.py`.
For sample output see `process_fgmax.html <process_fgmax.html>`_

Finally ::

    make plots

creates the usual plots and also adds links to the plot index for the fgmax
plots.

You can do everything with::

    make all

