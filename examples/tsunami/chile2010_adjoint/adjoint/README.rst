
.. _geoclaw_examples_tsunami_chile2010_adjoint_adjoint:

Adjoint solution for Chile 2010 with adjoint flagging
=====================================================

This directory has the adjoint solver that must be run first before
the forward solver in the enclosing directory.  

The adjoint solution is initialized with a hump of water at the location of
interest, in this case the location of
`DART buoy 32412 <http://www.ndbc.noaa.gov/station_page.php?station=32412>`__.

This is specified in the function `maketopo.qinit`.

