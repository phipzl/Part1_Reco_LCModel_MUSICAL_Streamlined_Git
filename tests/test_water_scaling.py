#!/usr/bin/env python3
"""Tests of water_scaling.py and fitter_calibration.py.

    python -m unittest discover -s tests

The round trips need LCModel with a basis (LCMODEL_BASIS, default the 3T one of the
container) and the deep fitters; they are skipped where those are missing.
"""

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import scipy.io

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import fitter_calibration as calibration  # noqa: E402
import water_scaling  # noqa: E402

DT = 1.06e-3
LCMODEL_BASIS = os.environ.get("LCMODEL_BASIS", "/opt/lcmodel_basis/3T/fid_1.000000ms.basis")
LCMODEL = shutil.which(os.environ.get("LCM_Path", "lcmodel")) or shutil.which("lcmodel")
HAS_DEEP = importlib.util.find_spec("deep_crt_mrsi") is not None
CONTROL = """ $LCMODL
 OWNER='test'
 key = 210387309
 Title='NoName_x1_y1_z1'
 HZPPPM=1.232469e+02, DELTAT=1.060000e-03, NUNFIL=288
 FILBAS='{basis}'
 DOREFS(1) = T
 DOREFS(2) = T
 NREFPK(2) = 1
 PPMREF(2,1) = 2.01
 DKNTMN = 0.3
 DOWS = T
 FILH2O='{dir}/NoName_x1_y1_z1_watref.RAW'
 PPMST = 4.2
 PPMEND = 1.8
 NOMIT = 15
 CHOMIT(1) = 'Cho'
 CHOMIT(2) = 'Act'
 CHOMIT(3) = 'mm3'
 CHOMIT(4) = 'mm4'
 CHOMIT(5) = 'Glc_B'
 CHOMIT(6) = 'Lip_c'
 CHOMIT(7) = 'Ala'
 CHOMIT(8) = 'Asp'
 CHOMIT(9) = 'ETA'
 CHOMIT(10) = 'Gly'
 CHOMIT(11) = 'Scyllo'
 CHOMIT(12) = 'Tau'
 CHOMIT(13) = 'TwoHG'
 CHOMIT(14) = 'Lac'
 CHOMIT(15) = 'Glc'
 LTABLE = 7
 FILTAB='{dir}/NoName_x1_y1_z1.table'
 FILRAW='{dir}/NoName_x1_y1_z1.RAW'
 FILPS='{dir}/NoName_x1_y1_z1.ps'
 $END
"""
RAW_HEADER = " $SEQPAR\n hzpppm=1.232469e+02\n $END\n $NMID\n fmtdat='(2e14.5)'\n $END\n"


def voigt(amplitude, lorentz_hz, gauss_fwhm_hz, freq_hz, n=288, dt=DT, delay=0.0):
    t = np.arange(n) * dt - delay
    g = np.pi * gauss_fwhm_hz / (2 * np.sqrt(np.log(2)))
    return amplitude * np.exp(
        -np.pi * lorentz_hz * t - (g * t) ** 2 + 2j * np.pi * freq_hz * t + 0.7j
    )


def matlab_first_order_delay(fid, te, dt):
    """The first-order phase correction MRSI_Reconstruction.m applies to the W1 water pass."""
    n = fid.size
    omega_step = 2 * np.pi / ((n - 1) * dt)
    correction = np.exp(-1j * te * omega_step * np.arange(-(n // 2), n - n // 2))
    return np.fft.ifft(np.fft.ifftshift(np.fft.fftshift(np.fft.fft(fid)) * correction))


class VoigtEstimator(unittest.TestCase):
    def test_lorentzian(self):
        q = water_scaling.fit_water(voigt(5.0, 8.0, 0.0, -12.0), DT)
        self.assertAlmostEqual(q["amplitude"], 5.0, delta=5e-5)
        self.assertAlmostEqual(q["fwhm_hz"], 8.0, delta=0.01)
        self.assertAlmostEqual(q["frequency_hz"], -12.0, delta=0.01)
        self.assertLess(q["residual"], 1e-6)

    def test_voigt_amplitude_and_width(self):
        fid = voigt(5.0, 6.0, 7.0, 3.0)
        q = water_scaling.fit_water(fid, DT)
        self.assertAlmostEqual(q["amplitude"], 5.0, delta=5e-5)
        # FWHM against the absorption spectrum of a long, finely sampled FID
        long_fid = voigt(1.0, 6.0, 7.0, 0.0, n=2**17, dt=1e-4) * np.exp(-0.7j)
        long_fid[0] /= 2
        spectrum = np.fft.fft(long_fid).real
        half = np.nonzero(spectrum >= spectrum.max() / 2)[0]
        half = np.where(half > 2**16, half - 2**17, half)
        half.sort()
        self.assertAlmostEqual(q["fwhm_hz"], (half[-1] - half[0] + 1) / (2**17 * 1e-4), delta=0.15)

    def test_distorted_first_points(self):
        fid = voigt(5.0, 10.0, 6.0, 0.0)
        fid[:3] = [0.1, 2.0, -1.0]
        self.assertAlmostEqual(water_scaling.fit_water(fid, DT)["amplitude"], 5.0, delta=5e-4)

    def test_too_short(self):
        with self.assertRaises(ValueError):
            water_scaling.fit_water(voigt(5.0, 7.0, 8.0, 10.0, n=32, dt=3.6e-4), 3.6e-4, t0=9e-3)

    def test_water_delayed_by_the_matlab_route(self):
        for dt, te, n in ((1.06e-3, 1.0e-3, 288), (3.6e-4, 1.3e-3, 840)):
            fid = matlab_first_order_delay(voigt(5.0, 7.0, 8.0, 10.0, n=n, dt=dt), te, dt)
            with self.subTest(dt=dt, te=te):
                amplitude = water_scaling.fit_water(fid, dt, t0=te)["amplitude"]
                self.assertAlmostEqual(amplitude, 5.0, delta=0.05)

    def test_noise(self):
        rng = np.random.default_rng(3)
        errors = []
        for _ in range(20):
            fid = voigt(1.0, 7.0, 6.0, 5.0)
            fid = fid + 0.005 * (rng.standard_normal(fid.size) + 1j * rng.standard_normal(fid.size))
            errors.append(water_scaling.fit_water(fid, DT)["amplitude"] - 1.0)
        self.assertLess(abs(np.mean(errors)), 0.005)
        self.assertLess(np.std(errors), 0.01)


class Inputs(unittest.TestCase):
    def test_tissue_on_one_grid(self):
        import nibabel as nib

        with tempfile.TemporaryDirectory() as out:
            paths = [os.path.join(out, f"{t}.nii.gz") for t in ("gm", "wm", "csf", "flipped")]
            for path in paths[:3]:
                nib.save(nib.Nifti1Image(np.zeros((4, 5, 6), np.float32), np.eye(4)), path)
            nib.save(
                nib.Nifti1Image(np.zeros((4, 5, 6), np.float32), np.diag([-1, 1, 1, 1])), paths[3]
            )
            pve, affine = water_scaling.load_tissue(paths[:3])
            self.assertEqual(pve.shape, (3, 4, 5, 6))
            with self.assertRaises(SystemExit):
                water_scaling.load_tissue([paths[0], paths[1], paths[3]])

    def test_sa2rage_b1(self):
        import json

        import nibabel as nib

        with tempfile.TemporaryDirectory() as out:
            path = os.path.join(out, "B1map.nii.gz")
            nib.save(nib.Nifti1Image(np.full((4, 5, 6), 720, np.int16), np.eye(4)), path)
            b1, affine, record = water_scaling.load_b1("sa2rage:" + path)
            np.testing.assert_allclose(b1, 0.8)
            self.assertEqual((record["unit"], record["reference_amplitude_v"]), (900.0, 270.0))
            with self.assertRaises(ValueError):
                water_scaling.load_b1("sa2rage:" + os.path.join(out, "B1map.json"))
            # Without the prefix a NIfTI is relative B1 and needs its own JSON
            with self.assertRaises(SystemExit):
                water_scaling.load_b1(path)
            with open(os.path.join(out, "B1map.json"), "w") as f:
                json.dump({"TxRefAmp": 250}, f)
            self.assertEqual(
                water_scaling.load_b1("sa2rage:" + path)[2]["reference_amplitude_v"], 250.0
            )
            # A map of relative B1 is not in the SA2RAGE unit
            relative = os.path.join(out, "relative.nii.gz")
            nib.save(nib.Nifti1Image(np.full((4, 5, 6), 0.8, np.float32), np.eye(4)), relative)
            with self.assertRaises(ValueError):
                water_scaling.load_b1("sa2rage:" + relative)

    def test_water_psf(self):
        met = dict(
            matrix=[50, 50, 31],
            fov_mm=[220.0, 220.0],
            slab_thickness_mm=140.0,
            slab_position=[0.0, 0.0, 0.0],
        )
        water = dict(met, matrix=[36, 36, 23])
        weight = water_scaling.water_psf(met, water, (50, 50, 31), True)
        self.assertEqual(int(weight.any(axis=(0, 1)).sum()), 23)
        # Another slab (7T) or a larger matrix: the water was matched to the metabolite PSF
        self.assertIsNone(
            water_scaling.water_psf(met, dict(water, slab_thickness_mm=110.0), (50, 50, 31), True)
        )
        self.assertIsNone(
            water_scaling.water_psf(met, dict(water, matrix=[64, 64, 23]), (50, 50, 31), True)
        )
        # A header without the phase FOV
        self.assertIsNone(
            water_scaling.water_psf(met, dict(water, fov_mm=[220.0, None]), (50, 50, 31), True)
        )
        # Reconstructed on the metabolite grid: the metabolite PSF, whatever its header says
        self.assertIsNone(
            water_scaling.water_psf(met, water, (50, 50, 31), True, route="metabolite_grid")
        )

    def test_water_route(self):
        import h5py

        water = np.zeros((4, 4, 1, 32), np.complex64)
        with tempfile.TemporaryDirectory() as out:
            path = os.path.join(out, "WaterReference.mat")
            # MAT.jl's v7.3 layout: char arrays as uint16 with their MATLAB class
            with h5py.File(path, "w") as f:
                route = f.create_dataset(
                    "WaterReferenceCSI/Par/Route",
                    data=np.array([[ord(c)] for c in "metabolite_grid"], np.uint16).T[..., None],
                )
                route.attrs["MATLAB_class"] = np.bytes_("char")
            self.assertEqual(water_scaling.water_route(out), "metabolite_grid")
            with h5py.File(path, "w") as f:
                f.create_dataset("WaterReferenceCSI/Par/AcquiredMatrix", data=[36.0, 36.0, 23.0])
            self.assertIsNone(water_scaling.water_route(out))
            for par, expected in (
                ({"Route": "metabolite_grid"}, "metabolite_grid"),
                ({"AcquiredMatrix": [36, 36, 23]}, None),
                (None, None),
            ):
                with self.subTest(par=par):
                    content = (
                        {"WaterReferenceOnCSIGrid": water}
                        if par is None
                        else {"WaterReferenceCSI": {"Data": water, "Par": par}}
                    )
                    scipy.io.savemat(path, content)
                    self.assertEqual(water_scaling.water_route(out), expected)

    def test_kspace_window(self):
        import h5py

        with tempfile.TemporaryDirectory() as out:
            path = os.path.join(out, "CombinedCSI.mat")
            with h5py.File(path, "w") as f:
                f.create_dataset("csi/RecoPar/Dwelltimes", data=[[360000.0]])
            # MATLAB's reconstruction writes none: Hamming
            self.assertEqual(water_scaling.kspace_window(path), "hamming")
            with h5py.File(path, "a") as f:
                f.create_dataset(
                    "csi/RecoPar/KSpaceWindow",
                    data=np.array([[ord(c)] for c in "hann"], np.uint16).T[..., None],
                )
            self.assertEqual(water_scaling.kspace_window(path), "hann")

    def test_single_slice_water(self):
        water = np.full((4, 4, 1, 32), 1 + 1j, np.complex64)
        for name, content in (
            ("WaterReferenceOnCSIGrid", water),
            ("WaterReferenceCSI", {"Data": water}),
        ):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as out:
                scipy.io.savemat(os.path.join(out, "WaterReference.mat"), {name: content})
                scipy.io.savemat(
                    os.path.join(out, "Parameters.mat"), {"Par": {"CSI": {"TEs": 1300.0}}}
                )
                loaded, source, first_sample = water_scaling.load_water(out)
                self.assertEqual((loaded.shape, source), ((4, 4, 1, 32), name))
                self.assertAlmostEqual(
                    first_sample, 1.3e-3 if name == "WaterReferenceOnCSIGrid" else 0
                )

    def test_lcmodel_groups_sum_their_members(self):
        table = {
            "NAA": (8.0, 3),
            "NAAG": (1.0, 9),
            "Cr": (3.0, 4),
            "PCr": (4.5, 5),
            "GPC": (1.0, 6),
            "PCh": (0.5, 9),
            "Cho": (0.3, 20),
            "Glu": (7.0, 5),
            "Ins": (5.0, 5),
            "Tau": (1.0, 9),
        }
        present = {name: 1.0 for name in ("NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", "Glu", "Ins")}
        self.assertEqual(
            calibration.lcmodel_groups(table, present),
            {"tNAA": 9.0, "tCr": 7.5, "tCho": 1.8, "Glu": 7.0, "mIns": 5.0},
        )
        del table["PCr"]
        self.assertNotIn("tCr", calibration.lcmodel_groups(table, present))


A_P = 0.004
LINES = {
    "Cr": [(3.027, 3), (3.913, 2)],
    "PCr": [(3.029, 3), (3.930, 2)],
    "NAA": [(2.008, 3), (2.49, 1), (2.67, 1)],
}


def write_basis(path, lines, hz=297.22, dt=1e-4, n=16384):
    """An LCModel .BASIS of lines {name: [(ppm, protons)]} with 8/s decay, 1 mM, a_p A_P."""
    t = np.arange(n) * dt
    with open(path, "w") as f:
        f.write(
            f" $SEQPAR\n HZPPPM = {hz}\n $END\n $BASIS1\n BADELT = {dt}\n NDATAB = {n}\n $END\n"
        )
        for name, peaks in lines.items():
            fid = sum(h * np.exp(2j * np.pi * (4.65 - ppm) * hz * t) for ppm, h in peaks)
            spectrum = np.fft.fft(A_P * fid * np.exp(-8.0 * t))
            f.write(
                f" $NMUSED\n $END\n $BASIS\n METABO = '{name}'\n CONC = 1.0\n ISHIFT = 0\n $END\n"
            )
            f.write("".join(f" {v.real:.12e} {v.imag:.12e}\n" for v in spectrum))


class SourceBasis(unittest.TestCase):
    """a_p of a deep basis made from an LCModel basis: the conversion's scale x the source's."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.source = os.path.join(self.dir, "source.basis")
        write_basis(self.source, LINES)
        header, entries = calibration.read_lcmodel_basis(self.source)
        t_source = np.arange(int(header[2]) // 2) * header[1]
        t = np.arange(840) * 3.6e-4
        self.names = list(LINES)
        columns = []
        for name in self.names:
            fid = np.fft.ifft(entries[name][2])[: t_source.size]
            lin = np.interp(t, t_source, fid.real) + 1j * np.interp(t, t_source, fid.imag)
            columns.append(25.0 * np.conj(lin) * np.exp(-22.0 * t))
        self.fids = np.stack(columns, axis=1)
        self.source_a_p = calibration.lcmodel_a_p(header, entries)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def deep_a_p(self, sources):
        return calibration.deep_a_p(self.fids, self.names, 3.6e-4, 297.22, sources)

    def test_source(self):
        self.assertAlmostEqual(self.source_a_p / A_P, 1.0, delta=1e-3)
        a_p, record = self.deep_a_p([self.source])
        self.assertEqual(record["a_p_method"], "source basis")
        self.assertAlmostEqual(a_p / (25.0 * self.source_a_p), 1.0, delta=1e-6)

    def test_no_or_other_source(self):
        other = os.path.join(self.dir, "other.basis")
        write_basis(other, dict(LINES, NAA=[(2.008, 3), (2.52, 1), (2.67, 1)]))
        for sources in ([], [other]):
            with self.subTest(sources=sources):
                a_p, record = self.deep_a_p(sources)
                self.assertEqual(record["a_p_method"], "singlet fit")
                self.assertAlmostEqual(a_p / (25.0 * A_P), 1.0, delta=0.05)


class FrequencyOffset(unittest.TestCase):
    def test_shifted_singlets(self):
        mhz, dt, n = 297.22, 3.6e-4, 840
        t = np.arange(n) * dt
        fids = np.stack([np.exp(2j * np.pi * (ppm - 4.7) * mhz * t) for ppm in (2.008, 3.027)], 1)
        # PHIVE's window as run_deepmrsi.py writes it: highest ppm first
        ppm = (4.7 + (196 + np.arange(252) - 512) / (dt * 1024 * mhz))[::-1]
        rng = np.random.default_rng(0)
        shape = np.exp(-np.pi * 8.0 * t - (np.pi * 2.0 * t) ** 2 / (4 * np.log(2)))
        voxels = [a * (fids @ [1.5, 1.0]) * shape * np.exp(2j * np.pi * 9.5 * t) for a in range(60)]
        spectra = np.exp(-2j * np.pi * np.outer((ppm - 4.7) * mhz, t)) @ np.stack(voxels, 1)
        spectra = spectra.T + 2.0 * (
            rng.standard_normal((60, 252)) + 1j * rng.standard_normal((60, 252))
        )
        spectra[:5] = 0
        record = calibration.frequency_offset(spectra, ppm, fids, ["NAA", "Cr"], dt, mhz * 1e6, n)
        self.assertAlmostEqual(record["hz"], 9.5, delta=0.2)
        self.assertEqual(record["voxels"], 27)


class RoundTrip(unittest.TestCase):
    """C_fit from one synthetic spectrum recovers the concentrations of another, with another
    mix, amplitude and noise."""

    FACTORS = dict(
        zip(calibration.TRUTH, (1.2, 0.8, 1.1, 0.9, 1.3, 0.7, 1.0, 1.25, 0.75, 1.1, 0.9, 1.0, 1.0))
    )
    MIX = {name: c * f for (name, c), f in zip(calibration.TRUTH.items(), FACTORS.values())}

    def check(self, fitted, estimated, present, c_fit, tolerance):
        """c_fit one constant, or one per group."""
        known = calibration.group_truth(present)
        for group, value in fitted.items():
            recovered = (c_fit[group] if isinstance(c_fit, dict) else c_fit) * value / estimated
            with self.subTest(group=group):
                self.assertAlmostEqual(recovered / known[group], 1.0, delta=tolerance)

    @unittest.skipUnless(
        LCMODEL and os.path.isfile(LCMODEL_BASIS), "LCModel or its basis not available"
    )
    def test_lcmodel(self):
        with tempfile.TemporaryDirectory() as work:
            control = os.path.join(work, "NoName_x1_y1_z1.control")
            with open(control, "w") as f:
                f.write(CONTROL.format(basis=LCMODEL_BASIS, dir=work))
            with open(os.path.join(work, "NoName_x1_y1_z1.RAW"), "w") as f:
                f.write(RAW_HEADER)
            c_fit, _ = calibration.calibrate_lcmodel(
                control, LCMODEL, os.path.join(work, "a"), water_scaling.water_amplitude, 30.0
            )
            with mock.patch.dict(calibration.TRUTH, self.MIX):
                result = calibration.fit_lcmodel_synthetic(
                    control,
                    LCMODEL,
                    os.path.join(work, "b"),
                    water_scaling.water_amplitude,
                    300.0,
                    seed=1,
                )
            self.check(*result[:3], c_fit, 0.02)

    def deep(self, fitter, info, points, tolerance, offset_hz=0.0, sources=()):
        """The test spectrum at offset_hz; the calibration measures the offset in the spectra
        the fitter fitted."""
        point = calibration.operating_point(fitter, offset_hz)
        point["mix"] = {name: c * self.FACTORS.get(name, 1.0) for name, c in point["mix"].items()}
        fitted, estimated, present, _, invivo = calibration.fit_deep_synthetic(
            fitter,
            info,
            water_scaling.water_amplitude,
            300.0,
            points,
            point,
            seed=1,
            sources=sources,
        )
        c_fit, inputs = calibration.calibrate_deep(
            fitter, info, water_scaling.water_amplitude, 30.0, points, invivo, sources
        )
        constants = inputs["C_per_metabolite"] if inputs["per_group"] else c_fit
        self.check(fitted, estimated, present, constants, tolerance)
        return inputs

    @unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi not installed")
    def test_spatialregu_3t(self):
        info = {"larmor_frequency": 123.246948e6, "dwelltime": 1.06e6}
        inputs = self.deep("gpufit", info, 288, 0.02, sources=[LCMODEL_BASIS])
        if os.path.isfile(LCMODEL_BASIS):
            self.assertEqual(inputs["a_p_method"], "source basis")

    @unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi not installed")
    def test_phive_7t(self):
        # In vivo spectra sit about 10 Hz above PHIVE's basis; constants built at 0 Hz read
        # these groups 1-7% high. The encoder's mix dependence leaves up to 3%.
        info = {"larmor_frequency": 297.22e6, "dwelltime": 3.6e5}
        inputs = self.deep("dlfit", info, 840, 0.03, offset_hz=10.0)
        self.assertAlmostEqual(inputs["frequency_offset"]["hz"], 10.0, delta=1.0)
        self.assertEqual(inputs["a_p_method"], "source basis")

    @unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi not installed")
    def test_spatialregu_basis_off_its_dwell(self):
        # 5% off the 3T basis dwell: gpu_fit fits a resampled basis, the synthetic is built on it
        from forD_gpu_fit import GPUFIT

        if len(GPUFIT.BASIS_SETS["3T"]) < 3:
            self.skipTest("this gpufit does not resample its basis")
        info = {"larmor_frequency": 123.246948e6, "dwelltime": 1.113e6, "gpufit_basis": "3T"}
        fids, _, dwell, _, path = calibration.deep_basis("gpufit", info, 288)
        self.assertAlmostEqual(dwell, 1.113e-3)
        self.assertEqual(fids.shape[0], 288)
        self.assertEqual(path, GPUFIT.BASIS_SETS["3T"][2])
        info["dwelltime"] = 1.06e6
        self.assertEqual(calibration.deep_basis("gpufit", info, 288)[2], 1.06e-3)


if __name__ == "__main__":
    unittest.main()
