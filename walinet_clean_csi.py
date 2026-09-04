#!/usr/bin/env python3
"""Remove water and lipids from a reconstructed CSI, for LCModel to fit.

    python walinet_clean_csi.py <tmp_dir> <model>

WALINET normally sits inside the deepmrsi quantification, which fits with dlfit or
gpufit and replaces LCModel. This runs only the removal and writes the result
back beside the reconstruction, so the ordinary LCModel path fits cleaned spectra
instead. The reconstruction stores FIDs and the removal takes spectra, so this
transforms in both directions with deepmrsi.py's convention.

The model expects field corrected data. Correcting is `-A Patref`, which runs
before this step and leaves its own record, so nothing is done about it here.

Reads  <out_path>/CombinedCSI.mat   csi.Data, as run_julia_reco.jl writes it
       <tmp_dir>/mask_brain.raw     the same mask MRSI_Reconstruction uses
Writes <out_path>/CombinedCSI.mat   csi.Data replaced by the cleaned FIDs
       <out_path>/CombinedCSI_beforeWalinet.mat  the original, kept once
"""

import os
import shutil
import sys

import combined_csi_io as csi_io
import numpy as np
import processing_record


def main(tmp_dir, model):
    out_path = csi_io.read_out_path(tmp_dir)
    combined = csi_io.combined_path(tmp_dir)

    csi = csi_io.load_csi(combined)
    mask = csi_io.load_mask(tmp_dir, csi.shape[:3], "walinet_clean_csi")
    print(f"walinet_clean_csi: csi {csi.shape}, {int(mask.sum())} voxels in the mask")

    was_corrected = bool(processing_record.read_record(combined).get("b0_corrected"))
    print(f"walinet_clean_csi: field corrected before the removal: {was_corrected}")

    clean = run_walinet(
        csi, mask, model, larmor_hz=csi_io.read_larmor_hz(combined, tmp_dir)
    )

    backup = os.path.join(out_path, "CombinedCSI_beforeWalinet.mat")
    if not os.path.isfile(backup):
        shutil.copy2(combined, backup)
        print(f"walinet_clean_csi: kept the original as {backup}")
    csi_io.save_csi(combined, clean)
    print(f"walinet_clean_csi: wrote the cleaned FIDs to {combined}")
    # Rewriting the file invalidates the record, which names the size and time of
    # the file it describes, so it is restated for whatever reads the file next.
    processing_record.write_record(combined, b0_corrected=was_corrected)


def run_walinet(csi, mask, model, larmor_hz=None):
    import walinet.package_config as conf
    import walinet.remove_water_and_lipids as rw

    if model not in (None, "", "off"):
        conf.PACKAGE_CONFIG.model_relative_path = conf.resolve_model_relative_path(model)
    print(f"walinet_clean_csi: model {conf.PACKAGE_CONFIG.model_relative_path}")

    # Cropping and the field check live in walinet so both routes into it refuse
    # the same things. "-L WALINET,3T" on a 7T scan is the case they exist for:
    # the 3T model merely makes the acquisition look long, and without the check
    # it is cut from 840 points to 288 and processed anyway.
    csi = rw.crop_to_supported_length(csi, larmor_hz=larmor_hz)

    # remove_water_and_lipids takes spectra: it inverts this transform to recover
    # the FID it feeds the network, so handing it a FID would run the network on
    # time-domain data. Same convention as deepmrsi.py.
    spectra = np.fft.fftshift(np.fft.fft(csi, axis=-1), axes=-1)
    clean_spectra = rw.remove_water_and_lipids(spectra, mask)
    clean = np.fft.ifft(
        np.fft.ifftshift(clean_spectra, axes=-1), axis=-1
    ).astype(np.complex64)

    # isfinite on the complex array directly: .view(np.float32) needs a
    # contiguous last axis and raises on anything WALINET returns as a view.
    if not np.all(np.isfinite(clean)):
        sys.exit("ERROR: WALINET returned non-finite values.")
    if np.allclose(clean[mask], csi[mask], rtol=1e-4, atol=1e-4):
        sys.exit("ERROR: WALINET left the masked voxels unchanged, inference did not run.")
    return clean


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
