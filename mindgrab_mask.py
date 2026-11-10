#!/usr/bin/env python3
"""Brain mask of the reconstruction's own reference scan, with MindGrab.

    python mindgrab_mask.py <tmp_dir>

This is `-m mindgrab`, the masking of the online FIRE route: MindGrab on the
magnitude of the reference scan, with that route's partial-volume threshold.

Reads  <out_path>/CombinedCSI.mat    image_FullFID.Data, csi.RecoPar.PatientAffineLPS
Writes <tmp_dir>/mask_mindgrab.raw   on the CSI grid; Part1 reads it like a mask given with -m
"""

import os
import sys

import combined_csi_io as csi_io
import h5py
import numpy as np

AFFINE = "csi/RecoPar/PatientAffineLPS"


def load_patient_affine(combined):
    """The voxel-to-patient affine of the CSI grid (LPS), or None when the file has none."""
    with h5py.File(combined, "r") as f:
        if AFFINE not in f:
            return None
        # h5py reads MATLAB's column-major matrix transposed.
        affine = np.asarray(f[AFFINE][()], dtype=np.float64).T
    return affine if affine.shape == (4, 4) else None


def mindgrab_brain_mask(reference, affine_lps):
    """Boolean mask on the grid of reference (x, y, z, t, cha), by the online route's masking."""
    import deep_crt_mrsi.package_config as fire_conf
    from deep_crt_mrsi.bet_mask import bet_mask
    from deep_crt_mrsi.masking.affine import lps_to_ras

    info = {
        "masking_method": "mindgrab",
        "mindgrab_pv": fire_conf.PACKAGE_CONFIG.mindgrab_pv,
        "mindgrab_model_dir": fire_conf.PACKAGE_CONFIG.mindgrab_model_dir,
        "makehomogeneous_sigma": fire_conf.PACKAGE_CONFIG.makehomogeneous_sigma,
    }
    mask, _, _ = bet_mask(reference, info, affine=lps_to_ras(affine_lps))
    return np.asarray(mask, dtype=bool)


def write_mask_raw(path, mask):
    """The mask as float32, column-major, the layout of Part1's own mask_brain.raw."""
    np.asarray(mask, dtype=np.float32).flatten(order="F").tofile(path)


def main(tmp_dir):
    combined = csi_io.combined_path(tmp_dir)

    reference, _ = csi_io.load_reference_and_dwelltime(combined)
    if reference is None:
        sys.exit("ERROR: -m mindgrab needs the reference scan, but CombinedCSI.mat carries none.")
    affine_lps = load_patient_affine(combined)
    if affine_lps is None:
        sys.exit(
            "ERROR: -m mindgrab needs the patient affine of the CSI grid, and CombinedCSI.mat "
            "carries none: MindGrab leaves the brain with a wrong one. Reconstruct with -S "
            "and a current MRSI.jl."
        )

    mask = mindgrab_brain_mask(reference, affine_lps)
    print(f"mindgrab_mask: reference {reference.shape}, {int(mask.sum())} voxels in the mask")
    if not mask.any():
        sys.exit("ERROR: the MindGrab mask is empty.")

    path = os.path.join(tmp_dir, "mask_mindgrab.raw")
    write_mask_raw(path, mask)
    print(f"mindgrab_mask: wrote {path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
