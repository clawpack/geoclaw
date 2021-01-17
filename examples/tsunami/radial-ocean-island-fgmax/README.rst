
.. _geoclaw_examples_tsunami_radial-ocean-island-fgmax:

Radial ocean with an island and fgmax grids
============================================

The topography is a synthetic radially symmetric ocean with a continental
shelf.  An island is added at one location on the shelf.  

The topography was adapted from a test problem in this paper:

   M. J. Berger, D. L. George, R. J. LeVeque, and K. T. Mandli,
   The GeoClaw software for depth-averaged flows with adaptive refinement
   Advances in Water Resources 34 (2011), pp. 1195-1206.
   `[link] <http://faculty.washington.edu/rjl/pubs/awr10/index.html>`__

To create the input data::

    make data
    make input

Running `make input` executes `make_input_files.py`, or you can instead run
the notebook `make_input_files.ipynb`,
from which the `make_input_files.py` script was generated using `nbconvert`.

See sample output in the gallery version linked below.

After running the code, e.g. by::

    make .output

you can make the additional fgmax plots via::

    make fgmax_plots

creates the fgmax plots by running `process_fgmax.py`.
Instead, you can run the notebook `process_fgmax.ipynb`, 
from which the `process_fgmax.py` script was generated using `nbconvert`.

See sample output in the gallery version linked below.

Finally ::

    make plots

creates the usual plots and also adds links to the plot index for the fgmax
plots. (This can be done before or after creating the fgmax plots.)

You can do everything with::

    make all

**Rendered versions of the Jupyter notebooks** for this example can 
be viewed from the `Clawpack gallery version of this file.
<http://www.clawpack.org/gallery/_static/geoclaw/examples/tsunami/radial-ocean-island-fgmax/README.html>`__


Version
-------

- New in Version 5.7.0.
- Some Python issues fixed in Version 5.7.1.
- The results changed very slightly in v5.8.0 due to changes in the 
  transverse Riemann solver.
