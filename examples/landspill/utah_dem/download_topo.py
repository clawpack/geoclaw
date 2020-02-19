#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Download topography file
"""
import os
import sys
import time

if sys.version_info.major == 2:
    from urllib import urlretrieve
elif sys.version_info.major == 3:
    from urllib.request import urlretrieve
else:
    raise ImportError("Unknown Python version.")

def reporthook(count, block_size, total_size):
    """Progress bar

    The code snippet is obtained from:
    https://blog.shichao.io/2012/10/04/progress_speed_indicator_for_urlretrieve_in_python.html
    """
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

if __name__ == "__main__":

    if os.path.isfile("./utah_dem_topo_3.txt"):
        print("utah_dem_topo_3.txt already exists. Skip downloading.")
    else:
        print("Downloading topography file: utah_dem_topo_3.txt ...")

        urlretrieve(
            "https://dl.dropboxusercontent.com/s/hhpebow2s81yzgo/utah_dem_topo_3.txt?dl=0",
            "./utah_dem_topo_3.txt", reporthook)

        print("\nFinish downloading.")
