#!/usr/bin/env python3
"""
run_deepmrsi.py: bridge between the MRSI pipeline and the deepmrsi quantification.

Called from step 7 of Part1_ProcessMRSI.sh when the -Q flag is set:
    python run_deepmrsi.py <tmp_dir> <output_dir>

Reads the NIfTI files that the reconstruction wrote to <out_path>/deepmrsi_inputs/,
calls process_deep_mrsi_offline() of the deep_crt_mrsi package and writes the
metabolite maps as NIfTI to <output_dir>.
"""

import json
import os
import sys

import numpy as np

if len(sys.argv) < 3:
    print(f"Usage: python {sys.argv[0]} <tmp_dir> <output_dir>", file=sys.stderr)
    sys.exit(1)

tmp_dir = sys.argv[1]
output_dir = sys.argv[2]


def find_inputs_dir(tmp_dir):
    """Find <out_path>/deepmrsi_inputs/, starting at tmp_dir."""
    candidate = os.path.join(tmp_dir, "deepmrsi_inputs")
    if os.path.isdir(candidate):
        return candidate
    # Common layout: out_path/tmp_*/
    candidate = os.path.join(os.path.dirname(tmp_dir), "deepmrsi_inputs")
    if os.path.isdir(candidate):
        return candidate
    # Otherwise take out_path out of the parameter file
    par_file = os.path.join(tmp_dir, "InitialParameters.m")
    if os.path.isfile(par_file):
        with open(par_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("out_path"):
                    # e.g.  out_path = '/some/path';
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        val = parts[1].strip().rstrip(";").strip().strip("'").strip('"')
                        candidate = os.path.join(val, "deepmrsi_inputs")
                        if os.path.isdir(candidate):
                            return candidate
    raise FileNotFoundError(
        f"No deepmrsi_inputs directory found near tmp_dir={tmp_dir!r}. "
        "The reconstruction has to write the deep learning inputs first."
    )


inputs_dir = find_inputs_dir(tmp_dir)
print(f"run_deepmrsi: reading the inputs from {inputs_dir}")

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
# Optional settings that deepmrsi understands, only passed on if they are there
for key in ("bet_f", "bet_g", "walinet", "lipidSuppression_beta",
            "use_prescan_for_masking", "writeWithoutSuppression",
            "makehomogeneous_sigma"):
    if key in meta:
        info[key] = meta[key]

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


csi_path = os.path.join(inputs_dir, "csi.nii.gz")
if not os.path.isfile(csi_path):
    print(f"ERROR: csi.nii.gz not found in {inputs_dir}", file=sys.stderr)
    sys.exit(1)
fid = load_nifti_complex(csi_path)
print(f"run_deepmrsi: fid shape = {fid.shape}")

musical_path = os.path.join(inputs_dir, "musical.nii.gz")
if not os.path.isfile(musical_path):
    print(f"ERROR: musical.nii.gz not found in {inputs_dir}", file=sys.stderr)
    sys.exit(1)
patref = load_nifti_complex(musical_path)
print(f"run_deepmrsi: patref shape = {patref.shape}")

prescan = None
prescan_path = os.path.join(inputs_dir, "prescan.nii.gz")
if os.path.isfile(prescan_path):
    prescan = load_nifti_complex(prescan_path)
    print(f"run_deepmrsi: prescan shape = {prescan.shape}")
else:
    print("run_deepmrsi: no prescan.nii.gz, deepmrsi masks on the patref instead")

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

process_deep_mrsi_offline(fid, patref, info, output_dir, uncomb_prescan=prescan)

print(f"run_deepmrsi: done, the metabolite maps are in {output_dir}")
