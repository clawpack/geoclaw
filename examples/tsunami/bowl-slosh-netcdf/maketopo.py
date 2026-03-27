#!/usr/bin/env python
"""
Module to create topo and qinit data files for this example.
"""

from pathlib import Path
import numpy as np
import clawpack.geoclaw.topotools as topotools

def maketopo(output_dir: Path) -> None:
    topography = topotools.Topography(topo_func=topo)
    topography.x = np.linspace(-2.0, 2.0, 200)
    topography.y = np.linspace(-2.0, 2.0, 200)
    topography.write(output_dir / "bowl.nc", Z_format="%22.15e")


def topo(x,y):
    """Parabolic bowl topography."""
    a = 1.0
    h0 = 0.1
    return h0*(x**2 + y**2)/a**2 - h0


if __name__=='__main__':
    maketopo(Path())
