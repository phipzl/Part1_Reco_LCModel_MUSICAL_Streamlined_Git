#!/usr/bin/env python3
"""What has already been applied to a reconstructed CombinedCSI.mat.

Steps after the reconstruction read the same file, and one of them correcting
the field shift is invisible to the next. walinet_clean_csi.py writes this
record and run_deepmrsi.py reads it, so the correction is applied once however
many stages the configuration chains together.

The record carries the size and modification time of the file it describes.
A reconstruction rewrites CombinedCSI.mat, and a record left by an earlier run
in the same output directory is then about data that no longer exists; the
mismatch makes it ignored rather than believed.
"""

import json
import os

# Filesystems disagree on modification time resolution, and a bind mount can
# round. One second is far below the gap between two runs.
_MTIME_TOLERANCE_S = 1.0


def record_path(combined):
    return os.path.splitext(combined)[0] + "_processing.json"


def write_record(combined, **applied):
    """Record what was applied to combined, keyed as the readers ask for it."""
    stat = os.stat(combined)
    record = dict(applied)
    record["describes"] = {"size": stat.st_size, "mtime": stat.st_mtime}
    with open(record_path(combined), "w") as f:
        json.dump(record, f, indent=2)
    return record


def read_record(combined):
    """The record for this exact file, or an empty one when it does not apply."""
    path = record_path(combined)
    if not os.path.isfile(path) or not os.path.isfile(combined):
        return {}
    try:
        with open(path) as f:
            record = json.load(f)
    except (OSError, ValueError):
        return {}
    stat = os.stat(combined)
    describes = record.get("describes") or {}
    if describes.get("size") != stat.st_size:
        return {}
    if abs(float(describes.get("mtime", 0)) - stat.st_mtime) > _MTIME_TOLERANCE_S:
        return {}
    return record


def _main(argv):
    """Record from a shell step, which knows what ran but not how to write it.

        python processing_record.py <tmp_dir> b0_corrected=true

    The MATLAB reconstruction aligns the frequency inside itself, so nothing
    python-side sees it happen; run_matlab.sh calls this once it has, and the
    deep fitting then reads the same record the Julia route leaves.
    """
    import combined_csi_io

    if len(argv) < 2:
        raise SystemExit(_main.__doc__)
    combined = combined_csi_io.combined_path(argv[0])
    applied = {}
    for pair in argv[1:]:
        key, _, value = pair.partition("=")
        applied[key] = value.lower() in ("true", "1", "yes")
    write_record(combined, **applied)
    print(f"processing_record: {record_path(combined)} says {applied}")


if __name__ == "__main__":
    import sys

    _main(sys.argv[1:])
