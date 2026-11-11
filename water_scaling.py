#!/usr/bin/env python3
"""Scale the fitters' results to the water reference, outside the fitters (-w W3).

    python water_scaling.py <tmp_dir> [--tissue GM,WM,CSF --b1 <B1 map>]

Per voxel c = C_fit * F / A_w: F the fitter's result, A_w the water amplitude at the
water FID's first sample and C_fit the fitter's constant (fitter_calibration.py; PHIVE
has one per metabolite group). The unit is the institutional unit of pure water,
55556 mM with two protons; relaxation, saturation and the tissue water content are not
corrected. With the tissue partial volumes and a B1 map, the groups tNAA, tCr, tCho, Glu
and mIns also get mmol/kg tissue water (concentration_model.py).

Reads  <out_path>/WaterReference.mat   the W1 water pass, on the metabolite grid
       <out_path>/spectra/*.table      LCModel, fitted without water scaling
       <out_path>/maps/Orig_<Fitter>/  the deep fitters' <met>_amp_map.nii.gz
       <out_path>/maps/SpecMap_PHIVEInput_{real,imag}.nii.gz   the spectra PHIVE fitted
Writes <out_path>/maps/Orig_<Fitter>/<met>_iu_map.nii.gz
       <out_path>/maps/WaterReference/Water_{amp,fwhm_Hz,shift_Hz,residual}_map.nii.gz
       <out_path>/maps/Orig_<Fitter>/<met>_mM_map.nii.gz   with --tissue and --b1
       <out_path>/maps/WaterReference/{GM,WM,CSF}_fraction_map, B1_relative_map, ...
       <out_path>/maps/WaterReference/water_scaling.json   C_fit with its inputs
"""

import argparse
import glob
import json
import os
import re
import sys

import combined_csi_io as csi_io
import concentration_model as model
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
# The slab profile is not modelled: partitions centred closer than this to the VOI edge get
# no mmol/kg.
EDGE_MARGIN_MM = 15.0
SA2RAGE_PREFIX = "sa2rage:"


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


def water_route(out_path):
    """How the Julia water pass brought the water onto the metabolite grid
    (WaterReferenceCSI Par.Route, "metabolite_grid" for another slab), else None."""
    path = os.path.join(out_path, "WaterReference.mat")
    if h5py.is_hdf5(path):
        with h5py.File(path, "r") as f:
            route = f.get("WaterReferenceCSI/Par/Route")
            return None if route is None else "".join(map(chr, np.ravel(route[()])))
    names = [name for name, _, _ in scipy.io.whosmat(path)]
    if "WaterReferenceCSI" not in names:
        return None
    par = scipy.io.loadmat(path, variable_names=["WaterReferenceCSI"])["WaterReferenceCSI"]
    if "Par" not in par.dtype.names or "Route" not in par["Par"][0, 0].dtype.names:
        return None
    return str(np.ravel(par["Par"][0, 0]["Route"][0, 0])[0])


def recopar(combined):
    """Dwell time (s), Larmor frequency (Hz) and point count of the reconstruction."""
    with h5py.File(combined, "r") as f:
        dwell = float(np.ravel(f["csi/RecoPar/Dwelltimes"][()])[0]) * 1e-9
        larmor = float(np.ravel(f["csi/RecoPar/LarmorFreq"][()])[0])
        points = int(f["csi/Data"].shape[0])
    return dwell, larmor, points


def kspace_window(combined):
    """The reconstruction's k-space window: csi.RecoPar.KSpaceWindow of the Julia route, else
    Hamming, MATLAB's."""
    with h5py.File(combined, "r") as f:
        value = f.get("csi/RecoPar/KSpaceWindow")
        return "hamming" if value is None else "".join(map(chr, np.ravel(value[()])))


def load_b1(path):
    """Relative B1, its affine and record: a Siemens flip angle map (DICOM), a NIfTI of
    relative B1 with a JSON beside it giving the reference amplitude it was measured with,
    or "sa2rage:<NIfTI>", the B1map series of an SA2RAGE (model.sa2rage_b1_map)."""
    if path.startswith(SA2RAGE_PREFIX):
        return model.sa2rage_b1_map(path[len(SA2RAGE_PREFIX) :])
    if not path.endswith((".nii", ".nii.gz")):
        return model.siemens_b1_map(path)
    sidecar = re.sub(r"\.nii(\.gz)?$", ".json", path)
    if not os.path.isfile(sidecar):
        sys.exit(f"ERROR: {sidecar} with the B1 map's reference_amplitude_v is missing.")
    with open(sidecar) as f:
        record = json.load(f)
    if "reference_amplitude_v" not in record:
        sys.exit(f"ERROR: {sidecar} has no reference_amplitude_v.")
    image = nib.load(path)
    return (
        np.asarray(image.dataobj, float).reshape(image.shape[:3]),
        image.affine,
        dict(record, file=path),
    )


def load_tissue(paths):
    """GM, WM, CSF partial volumes (3, ...) and their affine; all three on one grid."""
    images = [nib.load(path) for path in paths]
    if len(images) != 3:
        sys.exit(f"ERROR: --tissue needs GM,WM,CSF, got {len(images)} files.")
    for path, image in zip(paths[1:], images[1:]):
        if image.shape[:3] != images[0].shape[:3] or not np.allclose(
            image.affine, images[0].affine, atol=1e-3
        ):
            sys.exit(f"ERROR: {path} is not on the grid of {paths[0]}.")
    pve = np.stack([np.asarray(image.dataobj, float).reshape(image.shape[:3]) for image in images])
    return pve, images[0].affine


def water_psf(met, water, shape, hamming, route=None, window="hamming"):
    """The water scan's k-space weight on the metabolite grid when it was zero-filled onto it
    (same FOV, matrix not larger), else None: the water was matched to the metabolite PSF.
    route "metabolite_grid" (water_route): reconstructed with the metabolite PSF, so None."""
    if route == "metabolite_grid":
        return None
    same_fov = (
        met["matrix"] is not None
        and water["matrix"] is not None
        and None not in met["fov_mm"] + water["fov_mm"]
        and np.allclose(met["fov_mm"], water["fov_mm"], atol=0.1)
        and abs(met["slab_thickness_mm"] - water["slab_thickness_mm"]) < 0.1
        and np.allclose(met["slab_position"], water["slab_position"], atol=0.1)
    )
    if not same_fov or any(a > n for a, n in zip(water["matrix"], shape)):
        return None
    return model.mrsi_kspace_weight(
        shape, hamming=hamming, hamming_z=True, acquired=water["matrix"], window=window
    )


def concentration_correction(
    template,
    mask,
    tissue,
    b1,
    met_file,
    water_file,
    edge_margin_mm,
    min_brain=0.9,
    hamming=True,
    route=None,
    window="hamming",
):
    """{group: factor} from institutional units to mmol/kg tissue water (concentration_model),
    the maps it rests on and its record. Voxels outside the inner partitions of the VOI, with
    less than min_brain brain tissue in the excited volume of either signal, or no B1, get 0.
    hamming: the reconstruction's in-plane Hamming (-h); MRSI.jl applies its kz window always.
    route: the water pass's (water_route). window: the k-space window of both reconstructions."""
    shape, csi_affine = template.shape[:3], template.affine
    met, water = model.read_protocol(met_file), model.read_protocol(water_file)
    pve, pve_affine = load_tissue(tissue)
    weight = model.mrsi_kspace_weight(
        shape, hamming=hamming, hamming_z=True, acquired=met["matrix"], window=window
    )
    amounts, brain, excited = model.tissue_on_grid(
        pve, pve_affine, csi_affine, shape, weight, model.voi_profile(met, shape[2])
    )
    weight_water = water_psf(met, water, shape, hamming, route, window)
    if weight_water is None:
        amounts_water, brain_water = amounts, brain
    else:
        amounts_water, brain_water, _ = model.tissue_on_grid(
            pve, pve_affine, csi_affine, shape, weight_water, model.voi_profile(water, shape[2])
        )
    total = amounts.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        volume = np.where(total > 0, amounts / np.where(total > 0, total, 1), 0)
    b1_map, b1_affine, b1_record = load_b1(b1)
    b1_grid, b1_support = model.b1_on_grid(b1_map, b1_affine, csi_affine, shape)
    inner = model.inner_partitions(met, shape[2], edge_margin_mm)
    valid = (
        mask
        & inner[None, None, :]
        & (brain >= min_brain)
        & (brain_water >= min_brain)
        & (amounts[0] + amounts[1] > 0)
        & (b1_support > 0.5)
    )
    factors = model.correction_factors(
        amounts,
        b1_grid,
        met,
        water,
        b1_record["reference_amplitude_v"],
        volume_water=None if weight_water is None else amounts_water,
    )
    factors = {
        g: np.where(valid & np.isfinite(f), f, 0).astype(np.float32) for g, f in factors.items()
    }
    maps = {
        "GM_fraction": volume[0],
        "WM_fraction": volume[1],
        "CSF_fraction": volume[2],
        "brain_fraction": brain,
        "brain_fraction_water": brain_water,
        "excited_fraction": excited,
        "B1_relative": b1_grid,
        "flip_deg_metabolite": np.degrees(
            model.actual_flip(met, b1_grid, b1_record["reference_amplitude_v"])
        ),
        "flip_deg_water": np.degrees(
            model.actual_flip(water, b1_grid, b1_record["reference_amplitude_v"])
        ),
        "molal_valid": valid.astype(np.float32),
        **{f"molal_factor_{g}": f for g, f in factors.items()},
    }
    field = model.field_of(met["b0_t"])
    record = {
        "unit": "mmol/kg tissue water (molal)",
        "formula": "[M] = c_iu * 55510/55556 * sum_t V^w_t d_t R_W,t / sum_GM,WM V_t d_t R_M,t, "
        "V and V^w the tissue amounts the metabolite and the water signal see (with one PSF "
        "this is sum_t f_t R_W,t / (f_GM R_M,GM + f_WM R_M,WM), f_t = v_t d_t / sum v d), "
        "R = sin a (1 - E1) / (1 - cos a E1) * exp(-TE / T2), "
        "a = nominal flip * B1 * reference amplitude / B1 map's",
        "method": "Gasparovic 2006 with CSF correction, as Near 2021 eq. 4 and 5",
        "protons": "in C_fit: each fitter's basis carries its metabolites' protons, the water 2",
        "field": field,
        "metabolite_scan": met,
        "water_scan": water,
        "b1": b1_record,
        "tissue": {
            "files": list(tissue),
            "order": "GM, WM, CSF",
            "grid": "partial volumes x the VOI box along the partitions, filtered with the "
            "reconstruction's k-space weight (disk per partition, its window in-plane as -h, "
            "Hamming in kz unless none); fractions of the brain for the maps, amounts for the "
            "factors",
            "hamming_in_plane": bool(hamming),
            "kspace_window": window,
            "metabolite_acquired_matrix": met["matrix"],
            "water_psf": (
                f"its own, acquired {water['matrix']} zero-filled onto {list(shape)}"
                if weight_water is not None
                else (
                    "the metabolite's: reconstructed on the metabolite grid (route "
                    "metabolite_grid)"
                    if route == "metabolite_grid"
                    else "the metabolite's (matched)"
                )
            ),
        },
        "valid": {
            "edge_margin_mm": edge_margin_mm,
            "inner_partitions": np.flatnonzero(inner).tolist(),
            "min_brain_fraction": min_brain,
            "voxels": int(valid.sum()),
        },
        "constants": {
            "relaxation": model.CONSTANTS[field],
            "water_content": model.WATER_CONTENT,
            "water_molal_mmol_per_kg": model.H2O_MOLAL,
            "iu_water_mM": model.IU_WATER_MM,
        },
        "sources": model.SOURCES,
        "median_factor": {
            g: float(np.median(f[valid])) if valid.any() else None for g, f in factors.items()
        },
    }
    return factors, maps, record


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


def fitter_input(maps_dir, label, mask):
    """The spectra a deep fitter fitted, in the mask, and their ppm axis (run_deepmrsi.py)."""
    paths = [
        os.path.join(maps_dir, f"SpecMap_{label}Input_{part}.nii.gz") for part in ("real", "imag")
    ]
    if not all(os.path.isfile(path) for path in paths):
        raise RuntimeError(f"{paths[0]} or its imaginary part is missing")
    real, imag = (nib.load(path) for path in paths)
    spectra = np.asarray(real.dataobj, np.float32) + 1j * np.asarray(imag.dataobj, np.float32)
    header = real.header
    ppm = float(header["toffset"]) + np.arange(spectra.shape[-1]) * float(header["pixdim"][4])
    return spectra[mask], ppm


def main(tmp_dir, tissue=None, b1=None, edge_margin_mm=EDGE_MARGIN_MM):
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

    scale = np.where(mask & (a_w > 0), 1 / np.where(a_w > 0, a_w, 1), 0)
    factors, correction = {}, None
    if tissue and b1:
        met_file = par.get("csi_path") or [""]
        met_file = met_file[0] if isinstance(met_file, list) else met_file
        water_file = par.get("WaterReference_MethodAndFile", "").split(",")[-1].strip()
        factors, correction_maps, correction = concentration_correction(
            template,
            mask,
            tissue,
            b1,
            met_file,
            water_file,
            edge_margin_mm,
            hamming=par.get("hamming_flag", 0) == 1,
            route=water_route(out_path),
            window=kspace_window(combined),
        )
        for name, volume in correction_maps.items():
            write_map(water_dir, f"{name}_map", volume, template)
        inner = correction["valid"]["inner_partitions"]
        print(
            "water_scaling: mmol/kg tissue water on "
            f"{correction['valid']['voxels']} voxels, partitions "
            f"{f'{inner[0]}-{inner[-1]}' if inner else 'none'}, median factors "
            + ", ".join(f"{g} {v:.3f}" for g, v in correction["median_factor"].items() if v)
        )
    record = {
        "method": "W3",
        "formula": "c = C_fit * F / A_w",
        "unit": "institutional units: pure water 55556 mM with 2 protons, no relaxation, "
        "saturation or water-content correction; <met>_mM_map: see correction",
        "correction": correction,
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
        write_scaled("LCModel", fitted, c_fit, {}, scale, maps_dir, template, factors)
    bases = par.get("basis_path") or []
    for label, fitter in DEEP_FITTERS.items():
        folder = os.path.join(maps_dir, f"Orig_{label}")
        fitted = deep_amplitudes(folder)
        if not fitted:
            continue
        info = {"larmor_frequency": larmor, "dwelltime": dwell * 1e9}
        try:
            invivo = fitter_input(maps_dir, label, mask) if fitter == "dlfit" else None
            c_fit, inputs = calibration.calibrate_deep(
                fitter,
                info,
                water_amplitude,
                median_water,
                points,
                invivo=invivo,
                sources=[bases] if isinstance(bases, str) else bases,
            )
        except (ValueError, RuntimeError) as reason:
            print(f"water_scaling: no calibration for {label}, its maps stay unscaled: {reason}")
            continue
        record["fitters"][label] = dict(C_fit=c_fit, **inputs)
        per_group = inputs["C_per_metabolite"] if inputs["per_group"] else {}
        print(f"water_scaling: {label} a_p {inputs['a_p']:.6g} ({inputs['a_p_method']})")
        if "frequency_offset" in inputs:
            offset = inputs["frequency_offset"]
            print(
                f"water_scaling: {label} synthetic at {offset['hz']:+.2f} Hz, the offset of its "
                f"input (NAA {offset['NAA_hz']:+.2f}, Cr {offset['Cr_hz']:+.2f} Hz)"
            )
        write_scaled(label, fitted, c_fit, per_group, scale, maps_dir, template, factors)

    with open(os.path.join(water_dir, "water_scaling.json"), "w") as f:
        json.dump(record, f, indent=1)
    print(f"water_scaling: done, {os.path.join(water_dir, 'water_scaling.json')}")


def write_scaled(label, fitted, c_fit, per_group, scale, maps_dir, template, factors=None):
    """The <met>_iu_maps, with a group's own constant where per_group has one, and the
    <met>_mM_maps of the groups with a molal factor."""
    folder = os.path.join(maps_dir, f"Orig_{label}")
    molal = []
    for name, amplitude in fitted.items():
        c = per_group.get(calibration.group_of(name), c_fit)
        iu = c * amplitude * scale
        write_map(folder, f"{name}_iu_map", iu, template)
        group = model.group_of(name)
        if factors and group in factors:
            write_map(folder, f"{name}_mM_map", iu * factors[group], template)
            molal.append(name)
    if molal:
        print(f"water_scaling: {label} mmol/kg maps of {', '.join(molal)}")
    groups = "".join(f", {group} {c:.5g}" for group, c in per_group.items())
    print(f"water_scaling: {label} C_fit {c_fit:.5g}{groups}, {len(fitted)} maps in {folder}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("tmp_dir")
    parser.add_argument(
        "--tissue", help="GM,WM,CSF partial volume NIfTIs (e.g. FAST pve_1,pve_2,pve_0)"
    )
    parser.add_argument(
        "--b1",
        help="Siemens flip angle map (DICOM), relative B1 NIfTI with a JSON, or "
        "sa2rage:<B1map series NIfTI> (value 900 = nominal flip at 270 V)",
    )
    parser.add_argument(
        "--edge-margin-mm",
        type=float,
        default=EDGE_MARGIN_MM,
        help="partitions closer than this to the VOI edge get no mmol/kg",
    )
    args = parser.parse_args()
    if bool(args.tissue) != bool(args.b1):
        sys.exit("ERROR: mmol/kg needs both --tissue and --b1.")
    main(
        args.tmp_dir, args.tissue.split(",") if args.tissue else None, args.b1, args.edge_margin_mm
    )
