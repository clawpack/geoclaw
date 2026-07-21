#!/usr/bin/env python

r"""Workflow tools for meteorological-forcing / storm data.

This module hosts higher-level workflow helpers that operate on storm files but
are not part of the core object model.  Currently it provides
``make_multi_structure`` for splitting a multi-track ATCF file into individual
``Storm`` objects.  ``storm.py`` keeps an import shim so the historical import
path continues to work.
"""


def make_multi_structure(path):
    r"""Create a dictionary of Storm objects for ATCF files with multiple storm tracks in them
    """
    from clawpack.geoclaw.surge.storm import Storm

    with open(path, 'r') as f:
        lines = f.readlines()
        curTime = "test"
        curTrack = "test"
        os.mkdir("Clipped_ATCFs")
        stormDict = {}
        for line in lines:
            lineArr = line.split(", ")
            if curTime in lineArr[2]:
                if curTrack in lineArr[4]:
                    fileWrite.writelines(line)
                else:
                    fileWrite.close()
                    stormDict[curTime].update({curTrack: Storm(path=os.path.join(os.path.expandvars(
                        os.getcwd()), "Clipped_ATCFs", curTime, curTrack), file_format="ATCF")})
                    curTrack = lineArr[4]
                    fileWrite = open("Clipped_ATCFs/" +
                                     curTime + "/" + curTrack, 'w')
                    fileWrite.writelines(line)
            else:
                if curTime != "test":
                    fileWrite.close()
                    stormDict[curTime].update({curTrack: Storm(path=os.path.join(os.path.expandvars(
                        os.getcwd()), "Clipped_ATCFs", curTime, curTrack), file_format="ATCF")})
                curTime = lineArr[2]
                curTrack = lineArr[4]
                stormDict[curTime] = {}
                os.mkdir("Clipped_ATCFs/" + curTime)
                fileWrite = open("Clipped_ATCFs/" +
                                 curTime + "/" + curTrack, 'w')
                fileWrite.writelines(line)
    return stormDict
