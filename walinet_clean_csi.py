#!/usr/bin/env python3
"""Remove water and lipids from a reconstructed CSI, for LCModel to fit.

    python walinet_clean_csi.py <tmp_dir> <model> [--b0]

WALINET normally sits inside the deepmrsi quantification, which fits with dlfit or
gpufit and replaces LCModel. This runs only the removal and writes the result
back beside the reconstruction, so the ordinary LCModel path fits cleaned spectra
instead. The reconstruction stores FIDs and the removal takes spectra, so this
transforms in both directions with deepmrsi.py's convention.

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


def main(tmp_dir, model, b0=False):
    out_path = read_out_path(tmp_dir)
    combined = os.path.join(out_path, "CombinedCSI.mat")
    if not os.path.isfile(combined):
        sys.exit(f"ERROR: {combined} not found, run the reconstruction first.")

    csi = load_csi(combined)
    mask = load_mask(tmp_dir, csi.shape[:3])
    print(f"walinet_clean_csi: csi {csi.shape}, {int(mask.sum())} voxels in the mask")

    if b0:
        reference, dwelltime = load_reference_and_dwelltime(combined)
        if reference is None:
            sys.exit("ERROR: -0 asked for the B0 correction but CombinedCSI.mat carries no "
                     "image_FullFID or RecoPar. Reconstruct with a current Part1.")
        csi = apply_b0(csi, reference, dwelltime, mask)

    clean = run_walinet(csi, mask, model, larmor_hz=read_larmor_hz(tmp_dir))

    backup = os.path.join(out_path, "CombinedCSI_beforeWalinet.mat")
    if not os.path.isfile(backup):
        shutil.copy2(combined, backup)
        print(f"walinet_clean_csi: kept the original as {backup}")
    save_csi(combined, clean)
    print(f"walinet_clean_csi: wrote the cleaned FIDs to {combined}")


def read_larmor_hz(tmp_dir):
    """The acquisition frequency, so the model's field can be checked against it.

    Best effort: without it the length guard still catches a grossly wrong model,
    so a missing value degrades the check rather than blocking the run.
    """
    par_file = os.path.join(tmp_dir, "InitialParameters.json")
    if not os.path.isfile(par_file):
        return None
    try:
        with open(par_file) as f:
            par = json.load(f)
    except (OSError, ValueError):
        return None
    for key in ("LarmorFreq", "larmor_frequency", "larmor_freq"):
        value = par.get(key)
        if value:
            return float(value)
    return None


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

    Everything else in the file survives untouched. The FID axis may be shorter
    than what was read, because the model bounds the length it accepts; a fixed
    size dataset cannot shrink, so that case is replaced rather than assigned
    into, carrying the MATLAB attributes over so the file stays readable there.
    """
    flipped = data.transpose(range(data.ndim)[::-1])
    packed = np.empty(flipped.shape, dtype=[("real", "<f4"), ("imag", "<f4")])
    packed["real"] = flipped.real
    packed["imag"] = flipped.imag
    with h5py.File(path, "r+") as f:
        dset = f["csi/Data"]
        if dset.shape == packed.shape:
            dset[...] = packed
            return
        if dset.shape[1:] != packed.shape[1:] or dset.shape[0] < packed.shape[0]:
            sys.exit(f"ERROR: cleaned data is {packed.shape}, the file holds {dset.shape}.")
        attrs = dict(dset.attrs)
        del f["csi/Data"]
        new = f["csi"].create_dataset("Data", data=packed)
        for name, value in attrs.items():
            new.attrs[name] = value


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


def load_reference_and_dwelltime(path):
    """The uncombined reference and the dwelltime, for the B0 map.

    run_julia_reco.jl stores both beside the CSI: image_FullFID.Data is the
    uncombined PATREFSCAN, which is what the phase difference is measured on,
    and RecoPar.Dwelltimes is in nanoseconds.
    """
    with h5py.File(path, "r") as f:
        if "image_FullFID/Data" not in f or "csi/RecoPar/Dwelltimes" not in f:
            return None, None
        raw = f["image_FullFID/Data"][()]
        dwelltime = float(np.ravel(f["csi/RecoPar/Dwelltimes"][()])[0])
    if raw.dtype.names == ("real", "imag"):
        ref = (raw["real"] + 1j * raw["imag"]).astype(np.complex64)
    else:
        ref = np.asarray(raw, dtype=np.complex64)
    return ref.transpose(range(ref.ndim)[::-1]), dwelltime


def apply_b0(csi, reference, dwelltime, mask):
    """Correct the field shift, in the order the online route uses.

    deepmrsi corrects the combined FIDs immediately before handing them to the
    removal, so the same is done here rather than anywhere more convenient.
    """
    from deep_crt_mrsi.b0_correction import B0_correct_fids, calculate_B0

    b0 = calculate_B0(reference, dwelltime)
    print(f"walinet_clean_csi: B0 map median {np.median(b0):.2f} Hz, "
          f"5-95% {np.percentile(b0, 5):.1f}..{np.percentile(b0, 95):.1f} Hz")
    return B0_correct_fids(csi, b0, dwelltime, brainmask=mask)


def run_walinet(csi, mask, model, larmor_hz=None):
    import walinet.package_config as conf
    import walinet.remove_water_and_lipids as rw

    if model not in (None, "", "off"):
        conf.PACKAGE_CONFIG.model_relative_path = conf.resolve_model_relative_path(model)
    print(f"walinet_clean_csi: model {conf.PACKAGE_CONFIG.model_relative_path}")

    # Cropping and the field check live in walinet so both routes into it refuse
    # the same things. Doing the crop here was how "-L WALINET,3T" on a 7T scan
    # stopped being an error: the 3T model merely made the acquisition look long,
    # and it was quietly cut from 840 points to 288 and processed anyway.
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
    Args = [a for a in sys.argv[1:] if a != "--b0"]
    main(Args[0], Args[1] if len(Args) > 1 else None, b0="--b0" in sys.argv)
