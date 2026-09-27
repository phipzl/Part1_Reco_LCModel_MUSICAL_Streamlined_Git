#!/usr/bin/env python3
"""Scale the fitters' results to the water reference, outside the fitters (-w W3).

    python water_scaling.py <tmp_dir>

Per voxel c = C_fit * F / A_w: F the fitter's result, A_w the water amplitude at the
water FID's first sample and C_fit the fitter's constant (fitter_calibration.py). The
unit is the institutional unit of pure water, 55556 mM with two protons; relaxation,
saturation and the tissue water content are not corrected.

Reads  <out_path>/WaterReference.mat   the W1 water pass, on the metabolite grid
       <out_path>/spectra/*.table      LCModel, fitted without water scaling
       <out_path>/maps/Orig_<Fitter>/  the deep fitters' <met>_amp_map.nii.gz
Writes <out_path>/maps/Orig_<Fitter>/<met>_iu_map.nii.gz
       <out_path>/maps/WaterReference/Water_{amp,fwhm_Hz,shift_Hz,residual}_map.nii.gz
       <out_path>/maps/WaterReference/water_scaling.json   C_fit with its inputs
"""

import glob
import json
import os
import re
import sys

import combined_csi_io as csi_io
import fitter_calibration as calibration
import h5py
import nibabel as nib
import numpy as np
import scipy.io
from scipy.optimize import least_squares

# The first points of the reconstructed water are distorted on both routes.
FIT_FIRST_POINT = 3
FIT_SECONDS = 0.16
DEEP_FITTERS = {"SpatialRegu": "gpufit", "PHIVE": "dlfit"}


def fit_water(fid, dt, t0=0.0):
    """Voigt fit of one water FID: amplitude at t0, FWHM and frequency in Hz, relative residual.

    Model a exp(2 pi i f t - pi L t - (g t)^2) on the points from FIT_FIRST_POINT past t0 on.
    """
    first = FIT_FIRST_POINT + int(np.ceil(t0 / dt - 1e-6))
    last = min(fid.size, max(first + 16, int(round((t0 + FIT_SECONDS) / dt))))
    if last - first < 16:
        raise ValueError(f"a water FID of {fid.size} points leaves too few to fit past t0 = {t0} s")
    n = np.arange(first, last)
    t, y = n * dt, np.asarray(fid, complex)[n]
    f0 = np.angle(np.sum(y[1:10] * np.conj(y[:9]))) / (2 * np.pi * dt)
    lw0 = max(1.0, -np.polyfit(t[:30], np.log(np.abs(y[:30]) + 1e-30), 1)[0] / np.pi)

    def model(p):
        return np.exp(-np.pi * p[0] * t - (p[2] * t) ** 2 + 2j * np.pi * p[1] * t)

    def residual(p):
        m = model(p)
        r = (y - np.vdot(m, y) / np.vdot(m, m) * m) / np.abs(y).max()
        return np.concatenate([r.real, r.imag])

    lw, f, g = least_squares(residual, [lw0, f0, 10.0], x_scale=[5.0, 5.0, 10.0]).x
    m = model((lw, f, g))
    a = np.vdot(m, y) / np.vdot(m, m)
    gauss_fwhm = 2 * abs(g) * np.sqrt(np.log(2)) / np.pi
    return {
        "amplitude": abs(a) * np.exp(-np.pi * lw * t0 - (g * t0) ** 2),
        "fwhm_hz": 0.5346 * lw + np.sqrt(0.2166 * lw**2 + gauss_fwhm**2),
        "frequency_hz": f,
        "residual": np.linalg.norm(y - a * m) / np.linalg.norm(y),
    }


def water_amplitude(fid, dt):
    return fit_water(fid, dt)["amplitude"]


def load_water(out_path):
    """The water FIDs on the metabolite grid as (x, y, z, t), their source and first sample time.

    The MATLAB water pass leaves WaterReferenceOnCSIGrid, delayed by TE through its
    first-order phase correction; the Julia pass leaves WaterReferenceCSI, undelayed.
    """
    path = os.path.join(out_path, "WaterReference.mat")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: {path} not found; -w W3 needs the water pass.")
    if h5py.is_hdf5(path):
        with h5py.File(path, "r") as f:
            source = (
                "WaterReferenceOnCSIGrid" if "WaterReferenceOnCSIGrid" in f else "WaterReferenceCSI"
            )
            raw = f[source if source == "WaterReferenceOnCSIGrid" else "WaterReferenceCSI/Data"][()]
        water = (raw["real"] + 1j * raw["imag"]).transpose(range(raw.ndim)[::-1])
    else:
        names = [name for name, _, _ in scipy.io.whosmat(path)]
        source = (
            "WaterReferenceOnCSIGrid" if "WaterReferenceOnCSIGrid" in names else "WaterReferenceCSI"
        )
        loaded = scipy.io.loadmat(path, variable_names=[source])[source]
        water = loaded if source == "WaterReferenceOnCSIGrid" else loaded["Data"][0, 0]
    first_sample = 0.0
    if source == "WaterReferenceOnCSIGrid":
        par = scipy.io.loadmat(os.path.join(out_path, "Parameters.mat"), struct_as_record=False)
        te_us = float(np.ravel(par["Par"][0, 0].CSI[0, 0].TEs)[0])
        first_sample = te_us * 1e-6 if te_us < 10000 else 0.0
    return np.asarray(water, np.complex64), source, first_sample


def recopar(combined):
    """Dwell time (s), Larmor frequency (Hz) and point count of the reconstruction."""
    with h5py.File(combined, "r") as f:
        dwell = float(np.ravel(f["csi/RecoPar/Dwelltimes"][()])[0]) * 1e-9
        larmor = float(np.ravel(f["csi/RecoPar/LarmorFreq"][()])[0])
        points = int(f["csi/Data"].shape[0])
    return dwell, larmor, points


def concentration_correction(shape):
    """Per-voxel factor from institutional units to mM (saturation, tissue water, CSF): none yet."""
    return np.ones(shape, np.float32), None


def write_map(folder, name, volume, template):
    os.makedirs(folder, exist_ok=True)
    image = nib.Nifti1Image(np.asarray(volume, np.float32), template.affine, template.header)
    image.set_data_dtype(np.float32)
    nib.save(image, os.path.join(folder, name + ".nii.gz"))


def lcmodel_amplitudes(spectra_dir, shape):
    """{name: map} of the LCModel concentrations, from the tables named <...>_x<i>_y<j>_z<k>;
    correction terms such as -CrCH2 left out."""
    maps = {}
    for path in glob.glob(os.path.join(spectra_dir, "*.table")):
        index = re.search(r"_x(\d+)_y(\d+)_z(\d+)\.table$", path)
        if index is None:
            continue
        voxel = tuple(int(i) - 1 for i in index.groups())
        for name, (conc, _) in calibration.read_lcmodel_table(path).items():
            if not name.startswith("-"):
                maps.setdefault(name, np.zeros(shape, np.float32))[voxel] = conc
    return maps


def deep_amplitudes(folder):
    """{name: map} of a deep fitter's <name>_amp_map.nii.gz, ratios to tCr left out."""
    maps = {}
    for path in sorted(glob.glob(os.path.join(folder, "*_amp_map.nii.gz"))):
        name = os.path.basename(path)[: -len("_amp_map.nii.gz")]
        if not name.endswith("_tCr"):
            image = nib.load(path)
            maps[name] = np.asarray(image.dataobj, np.float32).reshape(image.shape[:3])
    return maps


def main(tmp_dir):
    out_path = csi_io.read_out_path(tmp_dir)
    combined = csi_io.combined_path(tmp_dir)
    maps_dir = os.path.join(out_path, "maps")
    template_path = os.path.join(maps_dir, "csi_template.nii.gz")
    template = nib.load(template_path if os.path.isfile(template_path) else template_path[:-3])
    with open(os.path.join(tmp_dir, "InitialParameters.json")) as f:
        par = json.load(f)
    dwell, larmor, points = recopar(combined)

    water, source, first_sample = load_water(out_path)
    shape = water.shape[:3]
    mask = csi_io.load_mask(tmp_dir, shape, "water_scaling")
    print(
        f"water_scaling: {source} {water.shape}, first sample at {first_sample * 1e3:.3f} ms, "
        f"{int(mask.sum())} voxels in the mask"
    )
    qc = {
        key: np.zeros(shape, np.float32)
        for key in ("amplitude", "fwhm_hz", "frequency_hz", "residual")
    }
    for voxel in zip(*np.nonzero(mask & np.any(water != 0, axis=-1))):
        for key, value in fit_water(water[voxel], dwell, first_sample).items():
            qc[key][voxel] = value
    a_w = qc["amplitude"]
    water_dir = os.path.join(maps_dir, "WaterReference")
    for key, name in (
        ("amplitude", "Water_amp_map"),
        ("fwhm_hz", "Water_fwhm_Hz_map"),
        ("frequency_hz", "Water_shift_Hz_map"),
        ("residual", "Water_residual_map"),
    ):
        write_map(water_dir, name, qc[key], template)
    median_water = float(np.median(a_w[mask]))
    print(
        f"water_scaling: water amplitude median {median_water:.4g}, FWHM median "
        f"{np.median(qc['fwhm_hz'][mask]):.2f} Hz"
    )

    correction, correction_source = concentration_correction(shape)
    scale = np.where(mask & (a_w > 0), correction / np.where(a_w > 0, a_w, 1), 0)
    record = {
        "method": "W3",
        "formula": "c = C_fit * F / A_w",
        "unit": "institutional units: pure water 55556 mM with 2 protons, no relaxation, "
        "saturation or water-content correction",
        "correction": correction_source,
        "water": {
            "file": os.path.join(out_path, "WaterReference.mat"),
            "source": source,
            "water_scan": par.get("WaterReference_MethodAndFile", "").split(",")[-1].strip(),
            "estimator": "Voigt a exp(2 pi i f t - pi L t - (g t)^2), magnitude at first sample",
            "fit_first_point": FIT_FIRST_POINT,
            "fit_seconds": FIT_SECONDS,
            "first_sample_s": first_sample,
            "dwell_s": dwell,
            "larmor_hz": larmor,
            "voxels": int(mask.sum()),
            "median_amplitude": median_water,
            "median_fwhm_hz": float(np.median(qc["fwhm_hz"][mask])),
            "shift_convention": "frequency of exp(2 pi i f t) in the stored FIDs",
        },
        "fitters": {},
    }

    work = os.path.join(tmp_dir, "water_scaling")
    tables = sorted(glob.glob(os.path.join(out_path, "spectra", "*.table")))
    if tables:
        controls = sorted(glob.glob(os.path.join(out_path, "spectra", "*.control")))
        if not controls:
            sys.exit("ERROR: LCModel tables without their controls, which the calibration needs.")
        lcmodel = par.get("LCM_Path") or "lcmodel"
        c_fit, inputs = calibration.calibrate_lcmodel(
            controls[0], lcmodel, os.path.join(work, "LCModel"), water_amplitude, median_water
        )
        record["fitters"]["LCModel"] = dict(C_fit=c_fit, **inputs)
        fitted = lcmodel_amplitudes(os.path.join(out_path, "spectra"), shape)
        write_scaled("LCModel", fitted, c_fit, scale, maps_dir, template)
    bases = par.get("basis_path") or []
    for label, fitter in DEEP_FITTERS.items():
        folder = os.path.join(maps_dir, f"Orig_{label}")
        fitted = deep_amplitudes(folder)
        if not fitted:
            continue
        info = {"larmor_frequency": larmor, "dwelltime": dwell * 1e9}
        try:
            c_fit, inputs = calibration.calibrate_deep(
                fitter,
                info,
                water_amplitude,
                median_water,
                points,
                sources=[bases] if isinstance(bases, str) else bases,
            )
        except (ValueError, RuntimeError) as reason:
            print(f"water_scaling: no calibration for {label}, its maps stay unscaled: {reason}")
            continue
        record["fitters"][label] = dict(C_fit=c_fit, **inputs)
        print(f"water_scaling: {label} a_p {inputs['a_p']:.6g} ({inputs['a_p_method']})")
        write_scaled(label, fitted, c_fit, scale, maps_dir, template)

    with open(os.path.join(water_dir, "water_scaling.json"), "w") as f:
        json.dump(record, f, indent=1)
    print(f"water_scaling: done, {os.path.join(water_dir, 'water_scaling.json')}")


def write_scaled(label, fitted, c_fit, scale, maps_dir, template):
    folder = os.path.join(maps_dir, f"Orig_{label}")
    for name, amplitude in fitted.items():
        write_map(folder, f"{name}_iu_map", c_fit * amplitude * scale, template)
    print(f"water_scaling: {label} C_fit {c_fit:.5g}, {len(fitted)} maps in {folder}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
