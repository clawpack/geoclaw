#!/usr/bin/env python

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import clawpack.pyclaw.gauges as gauges

old_gauge = gauges.GaugeSolution(0, Path() / "regression_data")
new_gauge = gauges.GaugeSolution(0, Path() / "PlaneWaveMultilayerTest_output")

print(old_gauge.q.shape)
print(new_gauge.q.shape)

for i in range(8):
    fig, ax = plt.subplots()
    ax.plot(old_gauge.t, old_gauge.q[i, :], "x", label="old")
    ax.plot(new_gauge.t, new_gauge.q[i, :], "+", label="new")
    ax.set_title(f"Field {i}")
    ax.legend()

plt.show()