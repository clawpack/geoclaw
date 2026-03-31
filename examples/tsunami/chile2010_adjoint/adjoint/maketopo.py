"""Download topo and create qinit data for the Chile 2010 adjoint example."""

from pathlib import Path
import os

import numpy as np

import clawpack.clawutil.data
from clawpack.geoclaw import topotools
from clawpack.geoclaw.util import haversine


try:
    CLAW = os.environ["CLAW"]
except KeyError as exc:
    raise RuntimeError("*** Must first set CLAW environment variable") from exc

# Scratch directory for storing topo and dtopo files.
scratch_dir = Path(CLAW) / "geoclaw" / "scratch"


def get_topo(output_dir: Path | None = None, makeplots: bool = False) -> None:
    """Retrieve the topo file from the GeoClaw repository."""
    topo_name = "etopo10min120W60W60S0S.asc"
    url = f"http://depts.washington.edu/clawpack/geoclaw/topo/etopo/{topo_name}"

    if output_dir is None:
        output_dir = scratch_dir
    else:
        output_dir = Path(output_dir)

    clawpack.clawutil.data.get_remote_file(
        url,
        output_dir=output_dir,
        file_name=topo_name,
        verbose=True,
    )

    if makeplots:
        import matplotlib.pyplot as plt

        topo = topotools.Topography(output_dir / topo_name, topo_type=2)
        topo.plot()
        figure_path = Path(topo_name).with_suffix(".png")
        plt.savefig(figure_path)
        print(f"Created {figure_path}")


def makeqinit(output_dir: Path | None = None, center: tuple[float, float] = (-86.392, -17.975)) -> None:
    """Create the qinit data file."""
    
    # Default value Initial data for adjoint is Gaussian hump around this
    # location DART 32412 location
    

    if output_dir is None:
        output_dir = Path()
    else:
        output_dir = Path(output_dir)

    output_file = output_dir / "hump.xyz"
    topotools.topo1writer(
        str(output_file),
        lambda x, y: qinit(x, y, center),
        center[0] - 1.5,
        center[0] + 1.5,
        center[1] - 1.5,
        center[1] + 1.5,
        201,
        201,
    )


def qinit(x, y, center):
    """Return a Gaussian hump centered at the target DART buoy location."""
    r = haversine(x, y, center[0], center[1])
    z = np.exp(-((r / 20e3) ** 2))
    return z


if __name__ == "__main__":
    get_topo()
    makeqinit()
