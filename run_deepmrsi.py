#!/usr/bin/env python3
"""
run_deepmrsi.py: the deep fitters of the online FIRE route as one step of Part1.

Called from step 7 of Part1_ProcessMRSI.sh for -l PHIVE, SpatialRegu or Both:
    python run_deepmrsi.py <tmp_dir> <out_path> [--fitting X] [--mask part1|deepmrsi]

Takes the spectra Part1's earlier steps left in CombinedCSI.mat (or a
deepmrsi_inputs/ directory), fits them with deep_crt_mrsi's fitter and writes, in
the layout of Part1's DeepLearning route:

  maps/Orig_<Fitter>/<met>_amp_map.nii        amplitude, as fitted
  maps/Orig_<Fitter>/<met>_tCr_amp_map.nii    ratio to tCr
  maps/Orig_<Fitter>/<met>_sd_map.nii         SpatialRegu posterior SD, percent
  maps/Orig_<Fitter>/<met>_crlb_map.nii       PHIVE CRLB, percent
  maps/SpecMap_<Fitter>{Input,Fit,Baseline}_{real,imag}.nii.gz
  maps/SpecMap_{Raw,Processed}_{real,imag}.nii.gz   before and after Part1's steps
  AlignFreq/B0map_Hz.nii                      the measured field map, with -A Patref

<Fitter> is SpatialRegu or PHIVE. Every map the fitter makes is written; the online
route's removal of individual maps does not apply here. Geometry from
maps/csi_template.nii(.gz); spectral axes as ppm = toffset + index * pixdim[4].
"""

import argparse
import json
import os
import sys

import numpy as np
import processing_record

parser = argparse.ArgumentParser(description="Quantify the reconstructed MRSI data with deepmrsi.")
parser.add_argument("tmp_dir", help="temporary directory of the current run")
parser.add_argument("output_dir", help="Part1's output directory; the maps go to <output_dir>/maps")
parser.add_argument("--fitting", choices=("dlfit", "gpufit", "both", "off"), default=None,
                    help="fitting backend, deepmrsi decides if it is not given")
parser.add_argument("--mask", choices=("part1", "deepmrsi"), default="part1",
                    help="which brain mask to fit: Part1's, cut to the excited volume "
                         "(default), or the one deepmrsi derives from the reference scan")
args = parser.parse_args()

tmp_dir = args.tmp_dir
output_dir = args.output_dir


def find_inputs_dir(tmp_dir):
    """Find <out_path>/deepmrsi_inputs/, starting at tmp_dir."""
    candidate = os.path.join(tmp_dir, "deepmrsi_inputs")
    if os.path.isdir(candidate):
        return candidate
    # Common layout: out_path/tmp_*/
    candidate = os.path.join(os.path.dirname(tmp_dir), "deepmrsi_inputs")
    if os.path.isdir(candidate):
        return candidate
    # Otherwise take out_path out of the parameter file. write_InitialParameters.sh
    # writes it as JSON next to the MATLAB one, so nothing has to be scanned here.
    par_file = os.path.join(tmp_dir, "InitialParameters.json")
    if os.path.isfile(par_file):
        with open(par_file) as f:
            out_path = json.load(f).get("out_path", "")
        if out_path:
            candidate = os.path.join(out_path, "deepmrsi_inputs")
            if os.path.isdir(candidate):
                return candidate
    raise FileNotFoundError(
        f"No deepmrsi_inputs directory found near tmp_dir={tmp_dir!r}. "
        "The reconstruction has to write the deep learning inputs first."
    )


try:
    inputs_dir = find_inputs_dir(tmp_dir)
except FileNotFoundError as missing_inputs:
    inputs_dir = None
    print(f"run_deepmrsi: no deepmrsi_inputs directory ({missing_inputs})")
else:
    print(f"run_deepmrsi: reading the inputs from {inputs_dir}")


def find_combined_csi(tmp_dir):
    """The reconstruction's CombinedCSI.mat, which carries the same data."""
    par_file = os.path.join(tmp_dir, "InitialParameters.json")
    if os.path.isfile(par_file):
        with open(par_file) as f:
            out_path = json.load(f).get("out_path", "")
        if out_path:
            candidate = os.path.join(out_path, "CombinedCSI.mat")
            if os.path.isfile(candidate):
                return candidate
    for candidate in (os.path.join(tmp_dir, "CombinedCSI.mat"),
                      os.path.join(os.path.dirname(tmp_dir), "CombinedCSI.mat")):
        if os.path.isfile(candidate):
            return candidate
    return None


meta = {}
if inputs_dir is not None:
    meta_path = os.path.join(inputs_dir, "deepmrsi_metadata.json")
    if not os.path.isfile(meta_path):
        print(f"ERROR: deepmrsi_metadata.json not found in {inputs_dir}", file=sys.stderr)
        sys.exit(1)
    with open(meta_path) as f:
        meta = json.load(f)
    info = {
        "dwelltime": meta["dwelltime"],                # ms
        "larmor_frequency": meta["larmor_frequency"],  # Hz
        "inplane_res": meta["inplane_res"],            # mm
        "fov_slice": meta["fov_slice"],                # mm
    }
else:
    info = None  # filled from CombinedCSI.mat below, together with the arrays
try:
    import nibabel as nib
except ImportError:
    print("ERROR: nibabel is not installed. Run: pip install nibabel", file=sys.stderr)
    sys.exit(1)


def load_nifti_complex(path):
    """Load a NIfTI file as a complex array.

    Three storage conventions are supported, the same ones as in deepmrsi.py:
    a complex data type, real and imaginary part as last dimension of size two,
    or real valued data (the imaginary part is then zero).
    """
    img = nib.load(path)
    data = np.asarray(img.dataobj)
    if np.iscomplexobj(data):
        return data.astype(np.complex64)
    if data.shape[-1] == 2:
        return (data[..., 0] + 1j * data[..., 1]).astype(np.complex64)
    return data.astype(np.complex64)


def require(path, what):
    if not os.path.isfile(path):
        print(f"ERROR: {what} not found at {path}", file=sys.stderr)
        sys.exit(1)
    return path


prescan = None
combined_read_from = None
if inputs_dir is not None:
    fid = load_nifti_complex(require(os.path.join(inputs_dir, "csi.nii.gz"), "csi.nii.gz"))
    patref = load_nifti_complex(require(os.path.join(inputs_dir, "musical.nii.gz"), "musical.nii.gz"))
    prescan_path = os.path.join(inputs_dir, "prescan.nii.gz")
    if os.path.isfile(prescan_path):
        prescan = load_nifti_complex(prescan_path)
        print(f"run_deepmrsi: prescan shape = {prescan.shape}")
    else:
        print("run_deepmrsi: no prescan.nii.gz, deepmrsi masks on the patref instead")
else:
    # deep_crt_mrsi already knows how to read a CombinedCSI.mat: it takes the
    # combined spectra from csi.Data, the uncombined water reference from
    # image_FullFID.Data and the acquisition parameters from csi.RecoPar. There
    # is no prescan in that file, so masking falls back to the patref.
    from deep_crt_mrsi.combined_csi import load_combined_csi

    combined = find_combined_csi(tmp_dir)
    if combined is None:
        print(
            "ERROR: neither a deepmrsi_inputs directory nor a CombinedCSI.mat "
            "was found. The reconstruction has to leave one of them behind.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"run_deepmrsi: reading the inputs from {combined}")
    fid, patref, info = load_combined_csi(combined)
    combined_read_from = combined

print(f"run_deepmrsi: fid shape = {fid.shape}")
print(f"run_deepmrsi: patref shape = {patref.shape}")

# Settings the fitters read, from the metadata when a reconstruction wrote one.
for key in ("bet_f", "bet_g", "fitting", "makehomogeneous_sigma"):
    if key in meta:
        info[key] = meta[key]
if args.fitting is not None:
    info["fitting"] = args.fitting

# Part1's own steps (-A for the field, -L for water and lipids) ran on this data
# before. The fitter only needs to know whether the field is corrected, since that
# decides its landmark search. The record describes CombinedCSI.mat, not a
# deepmrsi_inputs directory written separately.
if combined_read_from is not None and processing_record.read_record(
    combined_read_from
).get("b0_corrected"):
    print("run_deepmrsi: the reconstruction is B0 corrected")
    info["b0_corrected"] = True

try:
    import deep_crt_mrsi.deepmrsi as fire
    import deep_crt_mrsi.package_config as fire_conf
except ImportError as e:
    print(f"ERROR: could not import deep_crt_mrsi: {e}", file=sys.stderr)
    sys.exit(1)

info["online"] = False
for key, default in (("bet_f", fire_conf.PACKAGE_CONFIG.bet_f),
                     ("bet_g", fire_conf.PACKAGE_CONFIG.bet_g),
                     ("makehomogeneous_sigma", fire_conf.PACKAGE_CONFIG.makehomogeneous_sigma),
                     ("fitting", fire_conf.PACKAGE_CONFIG.fitting)):
    info.setdefault(key, default)

out_path = os.path.abspath(output_dir)
maps_dir = os.path.join(out_path, "maps")
template_path = next((p for p in (os.path.join(maps_dir, "csi_template.nii.gz"),
                                  os.path.join(maps_dir, "csi_template.nii"))
                      if os.path.isfile(p)), None)
if template_path is None:
    print(f"ERROR: no csi_template.nii(.gz) in {maps_dir}; the maps take their geometry from it",
          file=sys.stderr)
    sys.exit(1)
template = nib.load(template_path)
grid = tuple(int(n) for n in fid.shape[:3])
if tuple(template.shape[:3]) != grid:
    print(f"ERROR: csi_template is {template.shape[:3]}, the data {grid}", file=sys.stderr)
    sys.exit(1)


def read_part1_raw(path):
    """A Part1 .raw volume on the spectroscopy grid, x fastest, as (x, y, z)."""
    values = np.fromfile(path, dtype="<f4")
    if values.size != int(np.prod(grid)):
        return None
    return values.reshape(grid, order="F")


# Part1 masked the brain on the anatomical and cut the mask to the excited volume,
# so the deep fitters fit the voxels LCModel fits. Without one, the reference scan
# is masked the way the online route does it.
mask = None
if args.mask == "part1":
    for candidate in (os.path.join(maps_dir, "mask.raw"), os.path.join(tmp_dir, "mask_brain.raw")):
        if os.path.isfile(candidate):
            values = read_part1_raw(candidate)
            if values is not None:
                mask = values > 0.5
                print(f"run_deepmrsi: fitting the {int(mask.sum())} voxels of Part1's mask, {candidate}")
                break
if mask is None:
    print("run_deepmrsi: no Part1 mask for this grid, masking the reference scan with BET")
    mask, _, _ = fire.bet_mask(patref, info)
if not mask.any():
    print("run_deepmrsi: the mask is empty, fitting every voxel")
    mask[:] = True

fid4 = fire.ensure_4d(fid) if fire.is_already_combined(fid) else fire.coil_combine(fid, patref)


def write_map(folder, name, volume):
    os.makedirs(folder, exist_ok=True)
    image = nib.Nifti1Image(np.asarray(volume, dtype=np.float32).reshape(grid),
                            template.affine, template.header)
    image.set_data_dtype(np.float32)
    nib.save(image, os.path.join(folder, name + ".nii"))


def ppm_axis(model_grid, n):
    """ppm of each sample: the parent grid's centre sample at 4.7 ppm, rising with
    the index, which is how the viewer places the online route's spectra."""
    parent = int(model_grid["signal_length"])
    first = int(model_grid["interval_bounds"][0])
    step = 1.0 / (model_grid["dwelltime_s"] * parent * model_grid["reference_frequency_mhz"])
    return 4.7 + (first + np.arange(n) - (parent >> 1)) * step, step


def write_spectral_map(name, spectra, model_grid):
    """Real and imaginary part as two 4D files, the ppm axis in the header
    (toffset + index * pixdim[4]), as Part1's DeepLearning route writes them."""
    spectra = np.asarray(spectra)
    ppm, step = ppm_axis(model_grid, spectra.shape[3])
    for part, values in (("real", spectra.real), ("imag", spectra.imag)):
        header = template.header.copy()
        header.set_data_dtype(np.float32)
        image = nib.Nifti1Image(values.astype(np.float32), template.affine, header)
        image.header["pixdim"][4] = step
        image.header["toffset"] = float(ppm[0])
        image.header["descrip"] = b"ppm = toffset + index * pixdim[4]"
        nib.save(image, os.path.join(maps_dir, f"SpecMap_{name}_{part}.nii.gz"))


def full_band(fid_array):
    """A FID cube as spectra on its own full grid, in the online route's convention."""
    n = int(fid_array.shape[3])
    model_grid = {"signal_length": n, "interval_bounds": [0, n],
                  "dwelltime_s": float(info["dwelltime"]) * 1e-9,
                  "reference_frequency_mhz": float(info["larmor_frequency"]) * 1e-6}
    return fire._fid_to_spec(np.conj(fid_array)), model_grid


def map_file_name(fire_name):
    """The fitter's map name, after adjust_naming, in Part1's scheme."""
    if fire_name.endswith("_Hz"):
        return fire_name
    if fire_name.endswith("/tCr"):
        return fire_name[: -len("/tCr")] + "_tCr_amp_map"
    if fire_name.endswith(("_sd", "_crlb")):
        return fire_name + "_map"
    return fire_name + "_amp_map"


def map_values(fire_name, volume):
    """Amplitudes and ratios as fitted; uncertainties, a fraction from the fitter,
    in percent like LCModel's %SD; field maps in Hz."""
    if fire_name.endswith(("_sd", "_crlb")):
        return np.asarray(volume) * 100.0
    return volume


# The spectra before and after Part1's own steps: the reconstruction as it came out,
# when the Julia route left it behind, and what the fitter fits.
spectra, model_grid = full_band(fid4)
write_spectral_map("Processed", spectra, model_grid)
raw_path = os.path.join(out_path, "julia_csi.raw")
if os.path.isfile(raw_path):
    raw = np.fromfile(raw_path, dtype=np.complex64)
    n_raw = raw.size // int(np.prod(grid))
    if n_raw * int(np.prod(grid)) == raw.size:
        spectra, model_grid = full_band(raw.reshape(grid + (n_raw,), order="F"))
        write_spectral_map("Raw", spectra, model_grid)

b0_raw = os.path.join(out_path, "AlignFreq", "B0map_Hz.raw")
if os.path.isfile(b0_raw):
    b0 = read_part1_raw(b0_raw)
    if b0 is not None:
        write_map(os.path.join(out_path, "AlignFreq"), "B0map_Hz", b0)

resolved = fire._resolve_fitting(info)
fitter_names = {"dlfit": ["dlfit"], "gpufit": ["gpufit"], "both": ["dlfit", "gpufit"]}.get(resolved, [])
if not fitter_names:
    print(f"run_deepmrsi: fitting resolved to {resolved!r}, no maps written")
for fitter in fitter_names:
    make = fire._dlfit_fitter if fitter == "dlfit" else fire._gpufit_fitter
    label = fire.OUTWARD_NAMES[fitter]
    print(f"run_deepmrsi: fitting with {label}")
    images, names, fits, baselines, windowed, fit_grid = make(info, False)(fid4, mask)
    images, names = fire._drop_empty_landmark_map(images, names)
    names = list(names)
    fire.adjust_naming(names)
    folder = os.path.join(maps_dir, f"Orig_{label}")
    for volume, name in zip(images, names):
        write_map(folder, map_file_name(name), map_values(name, volume))
    for data, kind in ((windowed, "Input"), (fits, "Fit"), (baselines, "Baseline")):
        write_spectral_map(f"{label}{kind}", data, fit_grid)
    print(f"run_deepmrsi: {len(names)} maps in {folder}")

print(f"run_deepmrsi: done, maps in {maps_dir}")
