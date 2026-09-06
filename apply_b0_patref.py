#!/usr/bin/env python3
"""Correct the field shift on a reconstruction, from the reference scan.

    python apply_b0_patref.py <tmp_dir>

This is `-A Patref`, the estimator the online FIRE route uses: the map comes
from the phase evolution across the uncombined PATREFSCAN, so it needs no second
acquisition. MRSI_Reconstruction.m runs the same estimator inside the MATLAB
reconstruction; here it runs on the reconstruction the Julia route produced,
which is the equivalent point in the chain and, like MATLAB's, is before
anything downstream reads the file. A WALINET removal in particular expects
corrected data, because its model was trained on it.

Reads  <out_path>/CombinedCSI.mat   csi.Data and image_FullFID.Data
       <tmp_dir>/mask_brain.raw     the same mask MRSI_Reconstruction uses
Writes <out_path>/CombinedCSI.mat   csi.Data corrected in place
       <out_path>/AlignFreq/B0map_Hz.raw   the map, for comparison with MATLAB's
       <out_path>/CombinedCSI_processing.json  so nothing corrects it again
"""

import sys

import combined_csi_io as csi_io
import numpy as np
import processing_record


def main(tmp_dir):
    out_path = csi_io.read_out_path(tmp_dir)
    combined = csi_io.combined_path(tmp_dir)

    if processing_record.read_record(combined).get("b0_corrected"):
        print("apply_b0_patref: the reconstruction is already corrected, nothing to do")
        return

    reference, dwelltime = csi_io.load_reference_and_dwelltime(combined)
    if reference is None:
        sys.exit("ERROR: -A Patref needs the reference scan, but CombinedCSI.mat carries "
                 "no image_FullFID or RecoPar. Reconstruct with a current Part1.")

    csi = csi_io.load_csi(combined)
    mask = csi_io.load_mask(tmp_dir, csi.shape[:3], "apply_b0_patref")
    print(f"apply_b0_patref: csi {csi.shape}, reference {reference.shape}, "
          f"{int(mask.sum())} voxels in the mask")

    from deep_crt_mrsi.b0_correction import B0_correct_fids, calculate_B0

    b0 = calculate_B0(reference, dwelltime)
    in_mask = b0[mask]
    print(f"apply_b0_patref: B0 map median {np.median(in_mask):.2f} Hz, "
          f"5-95% {np.percentile(in_mask, 5):.1f}..{np.percentile(in_mask, 95):.1f} Hz")

    corrected = B0_correct_fids(csi, b0, dwelltime, brainmask=mask)
    csi_io.save_csi(combined, corrected)
    print(f"apply_b0_patref: wrote the corrected FIDs to {combined}")
    print(f"apply_b0_patref: wrote {csi_io.write_b0_map_raw(out_path, b0)}")

    # The deep quantification reads this same file and would otherwise correct it
    # a second time, from a reference scan this step does not modify.
    processing_record.write_record(combined, b0_corrected=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
