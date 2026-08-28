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

import numpy as np
import scipy.io as sio


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
    """csi.Data as (x, y, z, t), the layout run_julia_reco.jl writes."""
    mat = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
    if "csi" not in mat:
        sys.exit(f"ERROR: {path} has no csi variable.")
    data = np.asarray(mat["csi"].Data)
    if not np.iscomplexobj(data):
        sys.exit(f"ERROR: csi.Data in {path} is {data.dtype}, expected complex.")
    return data


def save_csi(path, data):
    sio.savemat(path, {"csi": {"Data": data}}, do_compression=True)


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
    if not np.all(np.isfinite(clean.view(np.float32))):
        sys.exit("ERROR: WALINET returned non-finite values.")
    if np.allclose(clean[mask], csi[mask], rtol=1e-4, atol=1e-4):
        sys.exit("ERROR: WALINET left the masked voxels unchanged, inference did not run.")
    return clean


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
