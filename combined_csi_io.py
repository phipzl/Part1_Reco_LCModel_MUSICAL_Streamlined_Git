#!/usr/bin/env python3
"""Reading and writing the reconstruction's CombinedCSI.mat.

Both post-reconstruction steps on the Julia route work on the same file, so the
h5py handling lives here rather than in each of them. Both pipelines write MAT
v7.3, which is HDF5 and which scipy cannot read at all.
"""

import json
import os
import sys

import h5py
import numpy as np


def read_out_path(tmp_dir):
    par_file = os.path.join(tmp_dir, "InitialParameters.json")
    if not os.path.isfile(par_file):
        sys.exit(f"ERROR: {par_file} not found.")
    with open(par_file) as f:
        out_path = json.load(f).get("out_path", "")
    if not out_path:
        sys.exit(f"ERROR: out_path missing in {par_file}.")
    return out_path


def combined_path(tmp_dir):
    """The reconstruction this run produced, or an exit if it is not there."""
    combined = os.path.join(read_out_path(tmp_dir), "CombinedCSI.mat")
    if not os.path.isfile(combined):
        sys.exit(f"ERROR: {combined} not found, run the reconstruction first.")
    return combined


def load_csi(path):
    """csi.Data as (x, y, z, t).

    h5py stores the dimensions in the opposite order and complex numbers as a
    real/imag compound, so both are undone here.
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
    """Write FIDs back into the same dataset.

    Everything else in the file survives untouched. The FID axis may be shorter
    than what was read, because a model bounds the length it accepts; a fixed
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
            sys.exit(f"ERROR: data is {packed.shape}, the file holds {dset.shape}.")
        attrs = dict(dset.attrs)
        del f["csi/Data"]
        new = f["csi"].create_dataset("Data", data=packed)
        for name, value in attrs.items():
            new.attrs[name] = value


def load_mask(tmp_dir, shape, label):
    """Same mask and same fallback as MRSI_Reconstruction and julia_write_lcm_files."""
    mask_path = os.path.join(tmp_dir, "mask_brain.raw")
    if not os.path.isfile(mask_path):
        print(f"{label}: no mask_brain.raw, using every voxel")
        return np.ones(shape, dtype=bool)
    mask = np.fromfile(mask_path, dtype=np.float32)
    if mask.size != int(np.prod(shape)):
        print(f"{label}: mask has {mask.size} values, expected "
              f"{int(np.prod(shape))}; using every voxel instead")
        return np.ones(shape, dtype=bool)
    # MATLAB wrote it column-major.
    mask = mask.reshape(shape, order="F") > 0
    return mask if mask.any() else np.ones(shape, dtype=bool)


def load_reference_and_dwelltime(path):
    """The uncombined reference and the dwelltime, for the field map.

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


def read_larmor_hz(combined, tmp_dir):
    """The acquisition frequency in Hz, so a model's field can be checked.

    The reconstruction stores it beside the data it wrote, as csi.RecoPar
    LarmorFreq, and that copy describes these FIDs. InitialParameters.json is
    only a fallback: write_InitialParameters.sh does not put the frequency there.
    """
    try:
        with h5py.File(combined, "r") as f:
            if "csi/RecoPar/LarmorFreq" in f:
                value = float(np.ravel(f["csi/RecoPar/LarmorFreq"][()])[0])
                if value:
                    return value
    except OSError:
        pass
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


def write_b0_map_raw(out_path, b0_hz):
    """The field map as float32, column-major, beside the MATLAB route's own.

    Both reconstructions write this file from the same estimator, so the two are
    directly comparable; validation/compare_b0_maps.py does that comparison.
    """
    align_dir = os.path.join(out_path, "AlignFreq")
    os.makedirs(align_dir, exist_ok=True)
    path = os.path.join(align_dir, "B0map_Hz.raw")
    np.asarray(b0_hz, dtype=np.float32).flatten(order="F").tofile(path)
    return path
