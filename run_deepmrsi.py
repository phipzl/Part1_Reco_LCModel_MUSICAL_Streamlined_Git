#!/usr/bin/env python3
"""
run_deepmrsi.py: bridge between the MRSI pipeline and the deepmrsi quantification.

Called from step 7 of Part1_ProcessMRSI.sh when the -Q flag is set:
    python run_deepmrsi.py <tmp_dir> <output_dir> [--fitting X] [--walinet_model Y]

Takes its input from whichever the reconstruction left behind:

  <out_path>/deepmrsi_inputs/   NIfTI files, when a reconstruction writes them
  <out_path>/CombinedCSI.mat    otherwise, read through deep_crt_mrsi.combined_csi

then calls process_deep_mrsi_offline() and writes the metabolite maps as NIfTI to
<output_dir>. Nothing in Part1 currently writes deepmrsi_inputs/, so in practice
the second is the path taken; the first is kept because it is the richer input
(it can carry a separate prescan) and costs nothing to keep.
"""

import argparse
import json
import os
import sys

import numpy as np
import processing_record

parser = argparse.ArgumentParser(description="Quantify the reconstructed MRSI data with deepmrsi.")
parser.add_argument("tmp_dir", help="temporary directory of the current run")
parser.add_argument("output_dir", help="directory the metabolic maps are written to")
parser.add_argument("--fitting", choices=("dlfit", "gpufit", "off"), default=None,
                    help="fitting backend, deepmrsi decides if it is not given")
# The names come from walinet itself rather than a copy here, which is how this
# list fell behind: it still offered legacy_7T/final_7T/final_3T after they were
# renamed, so the 7T and 3T that Part1's own usage text documents were rejected.
try:
    from walinet.package_config import MODEL_CHOICES as WALINET_MODEL_CHOICES
except ImportError:
    WALINET_MODEL_CHOICES = None

parser.add_argument("--walinet_model", choices=WALINET_MODEL_CHOICES, default=None,
                    help="WALINET model for the lipid suppression, deepmrsi decides if it is not given")
parser.add_argument("--b0_correction", choices=("true", "false"), default=None,
                    help="correct the field shift before fitting; Part1 passes -0 through here. "
                         "Ignored when the reconstruction already carries the correction")
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

# Optional settings that deepmrsi understands, only passed on if they are there
for key in ("bet_f", "bet_g", "walinet", "walinet_model", "fitting",
            "lipidSuppression_beta", "use_prescan_for_masking",
            "writeWithoutSuppression", "makehomogeneous_sigma", "b0_correction"):
    if key in meta:
        info[key] = meta[key]

# The command line wins over the metadata. If neither sets them, deepmrsi keeps
# its own defaults.
if args.fitting is not None:
    info["fitting"] = args.fitting
if args.walinet_model is not None:
    info["walinet_model"] = args.walinet_model
if args.b0_correction is not None:
    info["b0_correction"] = args.b0_correction == "true"

# A step between the reconstruction and here may already have corrected the
# field: the WALINET removal does, because its model is trained on corrected
# data. deepmrsi is told so and skips its own correction rather than shifting
# the FIDs a second time, and the fitter still learns the input is corrected.
# Only for the file that was actually read: the record describes CombinedCSI.mat
# and says nothing about a deepmrsi_inputs directory written separately.
if combined_read_from is not None and processing_record.read_record(
    combined_read_from
).get("b0_corrected"):
    print("run_deepmrsi: the reconstruction is already B0 corrected, not correcting again")
    info["b0_corrected"] = True

os.makedirs(output_dir, exist_ok=True)

try:
    from deep_crt_mrsi.deepmrsi import process_deep_mrsi_offline
except ImportError as e:
    print(f"ERROR: could not import deep_crt_mrsi: {e}", file=sys.stderr)
    print("Install it first: pip install -e [path to deep_crt_mrsi]", file=sys.stderr)
    sys.exit(1)

print(f"run_deepmrsi: calling process_deep_mrsi_offline, output to {output_dir}")
print(f"  dwelltime={info['dwelltime']} ms, larmor_frequency={info['larmor_frequency']} Hz")
print(f"  inplane_res={info['inplane_res']} mm, fov_slice={info['fov_slice']} mm")
print(f"  fitting={info.get('fitting', 'deepmrsi default')}, "
      f"walinet_model={info.get('walinet_model', 'deepmrsi default')}")
print(f"  b0_correction={info.get('b0_correction', 'deepmrsi default')}, "
      f"input already corrected={info.get('b0_corrected', False)}")

process_deep_mrsi_offline(fid, patref, info, output_dir, uncomb_prescan=prescan)

print(f"run_deepmrsi: done, the metabolite maps are in {output_dir}")
