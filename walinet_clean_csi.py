#!/usr/bin/env python3
"""Remove water and lipids from a reconstructed CSI, for LCModel to fit.

    python walinet_clean_csi.py <tmp_dir> <model>

WALINET normally sits inside the deepmrsi quantification, which fits with dlfit or
gpufit and replaces LCModel. This runs only the removal, in the time domain, and
writes the result back beside the reconstruction so the ordinary LCModel path
fits cleaned spectra instead. Nothing is transformed to the frequency domain and
back, so no FFT convention has to agree with anyone.

Reads  <out_path>/CombinedCSI.mat   csi.Data, as run_julia_reco.jl writes it
       <tmp_dir>/mask_brain.raw     the same mask MRSI_Reconstruction uses
Writes <out_path>/CombinedCSI.mat   csi.Data replaced by the cleaned FIDs
       <out_path>/CombinedCSI_beforeWalinet.mat  the original, kept once
"""

import json
import os
import shutil
import sys

import h5py
import numpy as np


def main(tmp_dir, model):
    out_path = read_out_path(tmp_dir)
    combined = os.path.join(out_path, "CombinedCSI.mat")
    if not os.path.isfile(combined):
        sys.exit(f"ERROR: {combined} not found, run the reconstruction first.")

    csi = load_csi(combined)
    mask = load_mask(tmp_dir, csi.shape[:3])
    print(f"walinet_clean_csi: csi {csi.shape}, {int(mask.sum())} voxels in the mask")

    clean = run_walinet(csi, mask, model)

    backup = os.path.join(out_path, "CombinedCSI_beforeWalinet.mat")
    if not os.path.isfile(backup):
        shutil.copy2(combined, backup)
        print(f"walinet_clean_csi: kept the original as {backup}")
    save_csi(combined, clean)
    print(f"walinet_clean_csi: wrote the cleaned FIDs to {combined}")


def read_out_path(tmp_dir):
    par_file = os.path.join(tmp_dir, "InitialParameters.json")
    if not os.path.isfile(par_file):
        sys.exit(f"ERROR: {par_file} not found.")
    with open(par_file) as f:
        out_path = json.load(f).get("out_path", "")
    if not out_path:
        sys.exit(f"ERROR: out_path missing in {par_file}.")
    return out_path


def load_csi(path):
    """csi.Data as (x, y, z, t).

    Both pipelines write MAT v7.3, which is HDF5 and which scipy cannot read at
    all. h5py can, but stores the dimensions in the opposite order and complex
    numbers as a real/imag compound, so both are undone here.
    """
    with h5py.File(path, "r") as f:
        if "csi/Data" not in f:
            sys.exit(f"ERROR: {path} has no csi/Data.")
        raw = f["csi/Data"][()]
    if raw.dtype.names != ("real", "imag"):
        sys.exit(f"ERROR: csi/Data in {path} is {raw.dtype}, expected complex.")
    data = (raw["real"] + 1j * raw["imag"]).astype(np.complex64)
    return data.transpose(range(data.ndim)[::-1])


def save_csi(path, data):
    """Write the cleaned FIDs back into the same dataset.

    In place, so every MATLAB attribute and everything else in the file survives
    untouched. Only the values change, and the shape has to match what was read.
    """
    packed = np.empty(data.transpose(range(data.ndim)[::-1]).shape,
                      dtype=[("real", "<f4"), ("imag", "<f4")])
    flipped = data.transpose(range(data.ndim)[::-1])
    packed["real"] = flipped.real
    packed["imag"] = flipped.imag
    with h5py.File(path, "r+") as f:
        dset = f["csi/Data"]
        if dset.shape != packed.shape:
            sys.exit(f"ERROR: cleaned data is {packed.shape}, the file holds {dset.shape}.")
        dset[...] = packed


def load_mask(tmp_dir, shape):
    """Same mask and same fallback as MRSI_Reconstruction and julia_write_lcm_files."""
    mask_path = os.path.join(tmp_dir, "mask_brain.raw")
    if not os.path.isfile(mask_path):
        print("walinet_clean_csi: no mask_brain.raw, cleaning every voxel")
        return np.ones(shape, dtype=bool)
    mask = np.fromfile(mask_path, dtype=np.float32)
    if mask.size != int(np.prod(shape)):
        print(f"walinet_clean_csi: mask has {mask.size} values, expected "
              f"{int(np.prod(shape))}; cleaning every voxel instead")
        return np.ones(shape, dtype=bool)
    # MATLAB wrote it column-major.
    mask = mask.reshape(shape, order="F") > 0
    return mask if mask.any() else np.ones(shape, dtype=bool)


def run_walinet(csi, mask, model):
    import walinet.package_config as conf
    import walinet.remove_water_and_lipids as rw

    if model not in (None, "", "off"):
        conf.PACKAGE_CONFIG.model_relative_path = conf.resolve_model_relative_path(model)
    print(f"walinet_clean_csi: model {conf.PACKAGE_CONFIG.model_relative_path}")

    clean = rw.remove_water_and_lipids(csi, mask)
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
