#!/usr/bin/env python

import os
import numpy
import datetime
import subprocess

import batch.batch

import clawpack.geoclaw.topotools as topotools

days2seconds = lambda days: days * 60.0**2 * 24.0

# def eta(x, y, A=1.0, sigma=50e3, x0=0.0):
#     return A * numpy.exp(-(x - x0)**2 / sigma**2)


class VortexJob(batch.batch.Job):

    def __init__(self, num_cells, rp_type):

        super(VortexJob, self).__init__()

        self.rp_type = rp_type

        self.type = "vortex_example"
        self.name = rp_type
        self.prefix = f"n{str(num_cells).zfill(4)}"
        self.executable = "xgeoclaw"

        # Create base data object
        import setrun
        self.rundata = setrun.setrun()

        self.rundata.clawdata.num_cells = [num_cells, num_cells]
        self.rundata.clawdata.output_format = 'binary'
        self.rundata.amrdata.max1d = 3000


    def __str__(self):
        output = super(VortexJob, self).__str__()
        output += f"  N: {self.rundata.clawdata.num_cells[0]}\n"
        output += f"  RP: {self.rp_type}\n"
        return output


    def write_data_objects(self):
        r""""""

        # Write out all data files
        super(VortexJob, self).write_data_objects()


if __name__ == '__main__':

    jobs = []
    num_cells = [2**n for n in range(6, 12)]
    # num_cells = [50, 100, 200, 400, 800, 1600]
    for rp_type in ['simple', 'geoclaw']:
        subprocess.run(['make', 'new', f'RP={rp_type}'])
        jobs = []
        for N in num_cells:
            jobs.append(VortexJob(N, rp_type))

        controller = batch.batch.BatchController(jobs)
        controller.wait = True
        controller.plot = False
        print(controller)
        controller.run()

        print(f"Done with {rp_type}!")

    # Run convergence script
    subprocess.run(['./plot_comparison.py', num_cells])
