#!/usr/bin/env python3
"""Tests of mindgrab_mask.py and of the calls Part1 makes into deep_crt_mrsi.

    python -m unittest discover -s tests

The tests that call deep_crt_mrsi are skipped where it is not installed.
"""

import importlib.util
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import combined_csi_io as csi_io  # noqa: E402
import mindgrab_mask  # noqa: E402

HAS_DEEP = importlib.util.find_spec("deep_crt_mrsi") is not None
GRID = (6, 5, 3)
AFFINE = np.array([[0, 2, 0, -5], [3, 0, 0, 7], [0, 0, 4, 1], [0, 0, 0, 1]], dtype=np.float64)


def as_matlab(data):
    """An array the way MATLAB v7.3 stores it: axes reversed, complex as a compound."""
    flipped = np.asarray(data).transpose(range(np.ndim(data))[::-1])
    if not np.iscomplexobj(flipped):
        return flipped
    packed = np.empty(flipped.shape, dtype=[("real", "<f4"), ("imag", "<f4")])
    packed["real"], packed["imag"] = flipped.real, flipped.imag
    return packed


def write_run(folder, affine=AFFINE):
    """A run directory with a CombinedCSI.mat; returns its temporary directory."""
    tmp_dir = os.path.join(folder, "tmp1")
    os.makedirs(tmp_dir)
    with open(os.path.join(tmp_dir, "InitialParameters.json"), "w") as f:
        json.dump({"out_path": folder}, f)
    rng = np.random.default_rng(0)
    reference = rng.normal(size=GRID + (4, 2)) + 1j * rng.normal(size=GRID + (4, 2))
    with h5py.File(os.path.join(folder, "CombinedCSI.mat"), "w") as f:
        f["csi/Data"] = as_matlab(reference[..., 0])
        f["csi/RecoPar/Dwelltimes"] = np.array([[360000.0]])
        f["image_FullFID/Data"] = as_matlab(reference)
        if affine is not None:
            f["csi/RecoPar/PatientAffineLPS"] = as_matlab(affine)
    return tmp_dir


class MindgrabMaskFile(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)

    def test_affine_is_read_as_written(self):
        tmp_dir = write_run(self.folder.name)
        affine = mindgrab_mask.load_patient_affine(csi_io.combined_path(tmp_dir))
        np.testing.assert_array_equal(affine, AFFINE)

    def test_mask_file_has_the_layout_of_mask_brain_raw(self):
        tmp_dir = write_run(self.folder.name)
        mask = np.zeros(GRID, dtype=bool)
        mask[1:4, 2, 0:2] = True
        with mock.patch.object(mindgrab_mask, "mindgrab_brain_mask", return_value=mask) as masker:
            mindgrab_mask.main(tmp_dir)
        self.assertEqual(masker.call_args[0][0].shape, GRID + (4, 2))
        np.testing.assert_array_equal(masker.call_args[0][1], AFFINE)
        os.replace(
            os.path.join(tmp_dir, "mask_mindgrab.raw"), os.path.join(tmp_dir, "mask_brain.raw")
        )
        np.testing.assert_array_equal(csi_io.load_mask(tmp_dir, GRID, "test"), mask)

    def test_no_affine_stops_the_run(self):
        tmp_dir = write_run(self.folder.name, affine=None)
        with mock.patch.object(mindgrab_mask, "mindgrab_brain_mask") as masker:
            with self.assertRaises(SystemExit) as stop:
                mindgrab_mask.main(tmp_dir)
        self.assertIn("patient affine", str(stop.exception))
        masker.assert_not_called()

    def test_empty_mask_stops_the_run(self):
        tmp_dir = write_run(self.folder.name)
        empty = np.zeros(GRID, dtype=bool)
        with mock.patch.object(mindgrab_mask, "mindgrab_brain_mask", return_value=empty):
            with self.assertRaises(SystemExit):
                mindgrab_mask.main(tmp_dir)
        self.assertFalse(os.path.isfile(os.path.join(tmp_dir, "mask_mindgrab.raw")))


@unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi is not installed")
class DeepCrtMrsiCalls(unittest.TestCase):
    def test_mindgrab_gets_the_affine_in_ras(self):
        reference = np.ones(GRID + (4, 1), dtype=np.complex64)
        mask = np.ones(GRID, dtype=bool)
        with mock.patch("deep_crt_mrsi.bet_mask.bet_mask", return_value=(mask, None, None)) as bet:
            result = mindgrab_mask.mindgrab_brain_mask(reference, AFFINE)
        np.testing.assert_array_equal(result, mask)
        self.assertEqual(bet.call_args[0][1]["masking_method"], "mindgrab")
        ras = np.diag([-1.0, -1.0, 1.0, 1.0]) @ AFFINE
        np.testing.assert_array_equal(bet.call_args[1]["affine"], ras)

    def test_patref_field_map_is_the_mean_phase_step(self):
        import apply_b0_patref

        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        tmp_dir = write_run(folder.name)
        combined = csi_io.combined_path(tmp_dir)
        offset_hz, dwell_s = 12.5, 360e-6
        t = np.arange(4) * dwell_s
        reference = np.ones(GRID + (1, 2)) * np.exp(2j * np.pi * offset_hz * t)[:, None]
        with h5py.File(combined, "r+") as f:
            del f["image_FullFID/Data"]
            f["image_FullFID/Data"] = as_matlab(reference.astype(np.complex64))
        apply_b0_patref.main(tmp_dir)
        b0 = np.fromfile(os.path.join(folder.name, "AlignFreq", "B0map_Hz.raw"), dtype="<f4")
        np.testing.assert_allclose(b0, offset_hz, rtol=1e-3)


if __name__ == "__main__":
    unittest.main()
