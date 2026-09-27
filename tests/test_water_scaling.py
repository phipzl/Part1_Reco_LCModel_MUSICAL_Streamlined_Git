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


class RoundTrip(unittest.TestCase):
    """C_fit from one synthetic spectrum recovers the concentrations of another, with another
    mix, amplitude and noise."""

    MIX = {
        name: c * f
        for (name, c), f in zip(
            calibration.TRUTH.items(),
            (1.2, 0.8, 1.1, 0.9, 1.3, 0.7, 1.0, 1.25, 0.75, 1.1, 0.9, 1.0, 1.0),
        )
    }

    def check(self, fitted, estimated, present, c_fit, tolerance):
        known = calibration.group_truth(present)
        for group, value in fitted.items():
            recovered = c_fit * value / estimated
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

    def deep(self, fitter, info, points, tolerance, sources=()):
        c_fit, inputs = calibration.calibrate_deep(
            fitter, info, water_scaling.water_amplitude, 30.0, points, sources
        )
        with mock.patch.dict(calibration.TRUTH, self.MIX):
            result = calibration.fit_deep_synthetic(
                fitter,
                info,
                water_scaling.water_amplitude,
                300.0,
                points,
                seed=1,
                sources=sources,
            )
        self.check(*result[:3], c_fit, tolerance)
        return inputs

    @unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi not installed")
    def test_spatialregu_3t(self):
        info = {"larmor_frequency": 123.246948e6, "dwelltime": 1.06e6}
        inputs = self.deep("gpufit", info, 288, 0.02, sources=[LCMODEL_BASIS])
        if os.path.isfile(LCMODEL_BASIS):
            self.assertEqual(inputs["a_p_method"], "source basis")

    @unittest.skipUnless(HAS_DEEP, "deep_crt_mrsi not installed")
    def test_phive_7t(self):
        # The network's bias depends on the mix by a few percent per metabolite.
        self.deep("dlfit", {"larmor_frequency": 297.22e6, "dwelltime": 3.6e5}, 840, 0.05)


if __name__ == "__main__":
    unittest.main()
