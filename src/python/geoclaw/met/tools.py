#!/usr/bin/env python

r"""Workflow tools for meteorological-forcing / storm data.

This module hosts higher-level workflow helpers that operate on storm files but
are not part of the core object model.  It provides ``iter_atcf`` for walking a
multi-track ATCF file in memory and ``make_multi_structure`` for the older
split-to-disk workflow (now a thin wrapper over ``iter_atcf``).  ``storm.py``
keeps an import shim so the historical import path continues to work.
"""

import os
import tempfile
from collections import OrderedDict


def _group_atcf_records(path):
    r"""Group an ATCF file's records by storm identity, in first-seen order.

    Identity is basin + cyclone number (the first two comma-separated fields of
    a b-deck record, e.g. ``AL`` and ``09``).  ``read_atcf`` splits on commas
    and strips, so this matches that.
    """
    storm_records = OrderedDict()
    with open(path, 'r') as data_file:
        for line in data_file:
            fields = line.split(",")
            if len(fields) < 2:
                continue
            storm_id = fields[0].strip() + fields[1].strip()
            storm_records.setdefault(storm_id, []).append(line)
    return storm_records


def iter_atcf(path):
    r"""Iterate the storms in a multi-storm ATCF file, without writing to disk.

    The in-memory counterpart to :func:`make_multi_structure`, and the companion
    of ``met.track.iter_hurdat`` / ``iter_ibtracs`` for the other two archive
    formats.

    :Input:
     - *path* (path-like) - the (possibly multi-storm) ATCF file to read.

    :Output:
     - Generator of ``(storm_id, Storm)`` pairs, where *storm_id* is basin plus
       cyclone number (e.g. ``"AL09"``), in the order the storms first appear.

    ``read_atcf`` takes a path rather than an open handle, so each storm's
    records are staged through a temporary file that is removed immediately
    afterwards; nothing is left in the working directory.
    """
    from clawpack.geoclaw.met.storm import Storm

    for storm_id, records in _group_atcf_records(path).items():
        with tempfile.TemporaryDirectory() as scratch:
            storm_path = os.path.join(scratch, "%s.storm" % storm_id)
            with open(storm_path, 'w') as out_file:
                out_file.writelines(records)
            yield storm_id, Storm(path=storm_path, file_format="ATCF")


def make_multi_structure(path, output_dir="Clipped_ATCFs"):
    r"""Split a multi-storm ATCF file into individual :class:`Storm` objects.

    Some ATCF files bundle several storms in one file, distinguished by their
    basin and cyclone number (the first two comma-separated fields of each
    b-deck record, e.g. ``AL`` and ``09``).  This groups the records by that
    identity, writes each storm's records to its own file under *output_dir*,
    and reads them back as :class:`Storm` objects.

    Use :func:`iter_atcf` instead when the per-storm files are not wanted; this
    function exists for the workflow where they are.

    :Input:
     - *path* (path-like) - the multi-storm ATCF file to split.
     - *output_dir* (path-like) - directory the per-storm ATCF files are
       written into (created if it does not exist).

    :Output:
     - (OrderedDict) mapping each storm id (basin + cyclone number, e.g.
       ``"AL09"``) to a :class:`Storm` read from that storm's records, in the
       order the storms first appear in the file.
    """
    from clawpack.geoclaw.met.storm import Storm

    os.makedirs(output_dir, exist_ok=True)

    # Write each storm's records to its own file and read it back as a Storm.
    storms = OrderedDict()
    for storm_id, records in _group_atcf_records(path).items():
        storm_path = os.path.join(output_dir, "%s.storm" % storm_id)
        with open(storm_path, 'w') as out_file:
            out_file.writelines(records)
        storms[storm_id] = Storm(path=storm_path, file_format="ATCF")

    return storms
