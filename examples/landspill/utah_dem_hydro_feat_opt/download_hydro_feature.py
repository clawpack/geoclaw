#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Download hydrologic features
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

    if os.path.isfile("./hydro_feature1.asc"):
        print("hydro_feature1.asc already exists. Skip downloading.")
    else:
        print("Downloading hydrologic feature: hydro_feature1.asc ...")

        urlretrieve(
            "https://dl.dropboxusercontent.com/s/02lbf2kg84b5sij/hydro_feature1.asc?dl=0",
            "./hydro_feature1.asc", reporthook)

        print("\nFinish downloading.")

    if os.path.isfile("./hydro_feature2.asc"):
        print("hydro_feature2.asc already exists. Skip downloading.")
    else:
        print("Downloading hydrologic feature: hydro_feature2.asc ...")

        urlretrieve(
            "https://dl.dropboxusercontent.com/s/mh8eh7wx6jy1z2t/hydro_feature2.asc?dl=0",
            "./hydro_feature2.asc", reporthook)

        print("\nFinish downloading.")

    if os.path.isfile("./hydro_feature3.asc"):
        print("hydro_feature3.asc already exists. Skip downloading.")
    else:
        print("Downloading hydrologic feature: hydro_feature3.asc ...")

        urlretrieve(
            "https://dl.dropboxusercontent.com/s/gs3g7amhatn9pxf/hydro_feature3.asc?dl=0",
            "./hydro_feature3.asc", reporthook)

        print("\nFinish downloading.")
