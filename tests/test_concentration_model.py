#!/usr/bin/env python3
"""Tests of concentration_model.py (mmol/kg tissue water for -w W3).

    python -m unittest discover -s tests

The header tests on real data need MRSI_3T_HEADER_DIR, a folder with a 3T CRT metabolite and
water dat pair and its Siemens B1 map (B1-Map-Siemens_18_MR/1.dcm); they are skipped without it.
"""

import os
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import concentration_model as model  # noqa: E402

MRSI_3T_HEADER_DIR = os.environ.get("MRSI_3T_HEADER_DIR", "")
MET_DAT = "meas_MID00049_FID24649_csi_ViennaDeepCRT_MRSI_50.dat"
WATER_DAT = "meas_MID00050_FID24650_water_reference_36.dat"


def protocol(tr, flip, te=1e-3, reference=263.0, b0=2.89362):
    return dict(tr_s=tr, flip_deg=flip, te_s=te, reference_amplitude_v=reference, b0_t=b0)


def steady_state_by_pulses(flip, tr, t1, n=3000):
    """Signal per M0 after n pulses, by stepping Mz through pulse and recovery."""
    mz = 1.0
    for _ in range(n):
        signal = mz * np.sin(flip)
        mz = 1 - (1 - mz * np.cos(flip)) * np.exp(-tr / t1)
    return signal


class Saturation(unittest.TestCase):
    def test_against_a_pulse_train(self):
        for flip, tr, t1 in (
            (np.deg2rad(50), 0.45, 1.35),
            (np.deg2rad(27), 0.26, 4.3),
            (np.deg2rad(5), 0.01, 0.8),
        ):
            self.assertAlmostEqual(
                model.saturation(flip, tr, t1), steady_state_by_pulses(flip, tr, t1), places=6
            )

    def test_limits(self):
        flip = np.deg2rad(50)
        self.assertAlmostEqual(model.saturation(flip, 1e3, 1.3), np.sin(flip), places=9)
        flips = np.deg2rad(np.linspace(1, 89, 8801))
        ernst = flips[np.argmax(model.saturation(flips, 0.45, 1.3))]
        self.assertAlmostEqual(ernst, np.arccos(np.exp(-0.45 / 1.3)), places=3)


class ActualFlip(unittest.TestCase):
    def test_each_scan_scales_by_its_own_reference_amplitude(self):
        b1 = np.array([0.8, 1.0, 1.2])
        met, water = protocol(0.46, 39, reference=246.7), protocol(0.26, 27, reference=295.0)
        np.testing.assert_allclose(
            model.actual_flip(met, b1, 270.0), np.deg2rad(39) * b1 * 246.7 / 270.0, rtol=1e-12
        )
        np.testing.assert_allclose(
            model.actual_flip(water, b1, 270.0)
            / model.actual_flip(dict(water, reference_amplitude_v=270.0), b1, 270.0),
            295.0 / 270.0,
            rtol=1e-12,
        )
        # The factors take the water scan's amplitude from its header, not the metabolite's
        volume = np.array([[0.5], [0.4], [0.1]])
        factors = [
            model.correction_factors(volume, np.ones(1), met, w, 270.0, field="7T")["tNAA"][0]
            for w in (water, dict(water, reference_amplitude_v=246.7))
        ]
        self.assertNotAlmostEqual(factors[0] / factors[1], 1.0, places=3)


class MolalFactor(unittest.TestCase):
    content = {"GM": 0.78, "WM": 0.65, "CSF": 0.97}

    def test_pure_tissue_ignores_water_content(self):
        for volume in ([1, 0, 0], [0, 1, 0]):
            f = model.molal_factor(
                np.array(volume, float), [0.4, 0.3, 0.2], [0.5, 0.25], self.content
            )
            r_water, r_met = (0.4, 0.5) if volume[0] else (0.3, 0.25)
            self.assertAlmostEqual(
                f, model.H2O_MOLAL / model.IU_WATER_MM * r_water / r_met, places=9
            )

    def test_csf(self):
        # Equal relaxation everywhere: the factor is 1 / (1 - f_CSF) of the water mole fraction
        for v_csf in (0.0, 0.3, 0.9):
            volume = np.array([1 - v_csf, 0, v_csf])
            f = model.molal_factor(volume, [1, 1, 1], [1, 1], self.content)
            f_csf = v_csf * 0.97 / ((1 - v_csf) * 0.78 + v_csf * 0.97)
            self.assertAlmostEqual(f, model.H2O_MOLAL / model.IU_WATER_MM / (1 - f_csf), places=9)
        self.assertFalse(
            np.isfinite(model.molal_factor(np.array([0, 0, 1.0]), [1, 1, 1], [1, 1], self.content))
        )

    def test_forward_model(self):
        """Known molal concentrations through signal formation and back, per voxel: tissue
        mixtures with CSF, the water signal seeing another mixture than the metabolite signal
        (its own PSF), B1 0.7-1.3, metabolite and water scans with other TR, flip and
        transmitter reference, as at 7T."""
        rng = np.random.default_rng(0)
        n = 500
        volume = 0.98 * rng.dirichlet([3, 3, 1], n).T
        volume_water = 0.95 * rng.dirichlet([3, 3, 1], n).T
        b1 = rng.uniform(0.7, 1.3, n)
        reference_b1 = 300.0
        for field, met, water in (
            ("3T", protocol(0.45, 50, reference=263.0), protocol(0.45, 50, reference=263.0)),
            (
                "7T",
                protocol(0.46, 39, 1.3e-3, 303.6, 6.98),
                protocol(0.26, 27, 1.3e-3, 295.0, 6.98),
            ),
        ):
            table = model.values(model.CONSTANTS[field])
            content = model.values(model.WATER_CONTENT)
            a_met = np.deg2rad(met["flip_deg"]) * b1 * met["reference_amplitude_v"] / reference_b1
            a_water = (
                np.deg2rad(water["flip_deg"]) * b1 * water["reference_amplitude_v"] / reference_b1
            )
            protons, water_molal = 3, model.H2O_MOLAL
            s_water = sum(
                2
                * water_molal
                * volume_water[k]
                * content[t]
                * steady_state_by_pulses(a_water, water["tr_s"], table["water_T1"][t])
                * np.exp(-water["te_s"] / table["water_T2"][t])
                for k, t in enumerate(model.TISSUES)
            )
            factors = model.correction_factors(
                volume, b1, met, water, reference_b1, field=field, volume_water=volume_water
            )
            for group in model.GROUPS:
                truth = rng.uniform(1, 12)
                s_met = sum(
                    protons
                    * truth
                    * volume[k]
                    * content[t]
                    * steady_state_by_pulses(a_met, met["tr_s"], table["T1"][group][t])
                    * np.exp(-met["te_s"] / table["T2"][group][t])
                    for k, t in enumerate(("GM", "WM"))
                )
                c_iu = model.IU_WATER_MM * 2 / protons * s_met / s_water
                np.testing.assert_allclose(
                    c_iu * factors[group], truth, rtol=1e-6, err_msg=f"{field} {group}"
                )

    def test_one_psf_is_the_mole_fraction_form(self):
        volume = np.random.default_rng(1).dirichlet([3, 3, 1], 50).T
        r_water, r_met = np.reshape([0.4, 0.3, 0.2], (3, 1)), [0.5, 0.25]
        f = model.water_mole_fractions(volume, self.content)
        expected = (
            model.H2O_MOLAL
            / model.IU_WATER_MM
            * (f * r_water).sum(0)
            / (f[0] * r_met[0] + f[1] * r_met[1])
        )
        for amounts in (volume, 0.9 * volume):
            for volume_water in (None, amounts):
                np.testing.assert_allclose(
                    model.molal_factor(amounts, r_water, r_met, self.content, volume_water),
                    expected,
                    rtol=1e-12,
                )

    def test_bounds_of_the_budget(self):
        volume = np.array([[0.6], [0.3], [0.1]])
        met = water = protocol(0.45, 50)
        central = model.correction_factors(volume, np.ones(1), met, water, 263.0)
        longer = model.correction_factors(
            volume,
            np.ones(1),
            met,
            water,
            263.0,
            which={("T1", "tNAA", "GM"): "high", ("T1", "tNAA", "WM"): "high"},
        )
        # A longer metabolite T1 saturates the metabolite more, so its concentration is higher
        self.assertGreater(longer["tNAA"][0], central["tNAA"][0])
        self.assertEqual(longer["tCr"][0], central["tCr"][0])


class Grid(unittest.TestCase):
    shape = (20, 20, 11)

    def setUp(self):
        self.csi = np.diag([4.4, 4.4, 140 / 11, 1.0])
        self.csi[:3, 3] = [-44, -44, -70]
        self.fine = np.diag([1.0, 1.0, 1.0, 1.0])
        self.fine[:3, 3] = [-60, -60, -90]
        z = np.arange(180) - 90
        self.x = (np.arange(120) - 60)[:, None, None] * np.ones((1, 120, 180))
        self.z = z[None, None, :] * np.ones((120, 120, 1))

    def test_kspace_weight(self):
        weight = model.mrsi_kspace_weight(self.shape)
        self.assertEqual(weight[0, 0, 0], 1)
        self.assertTrue(np.all(weight >= 0) and np.all(weight <= 1))
        # The outer partitions hold only the two innermost circles
        nz = self.shape[2]
        self.assertEqual(
            np.count_nonzero(weight[:, :, nz // 2]),
            np.count_nonzero(
                np.fft.fftfreq(20)[:, None] ** 2 + np.fft.fftfreq(20)[None, :] ** 2 < (2 / 20) ** 2
            ),
        )

    def test_kspace_weight_zero_filled(self):
        weight = model.mrsi_kspace_weight(self.shape, acquired=(12, 12, 7))
        planes = np.flatnonzero(weight.any(axis=(0, 1)))
        index = np.fft.fftfreq(self.shape[2], 1 / self.shape[2])
        self.assertEqual(sorted(index[planes]), list(range(-3, 4)))
        # The central partition's six circles reach 6 / 12 cycles per acquired voxel
        kx = np.fft.fftfreq(20, 1 / 20)
        self.assertEqual(np.abs(kx[np.flatnonzero(weight[:, 0, 0])]).max(), 5)

    def test_kspace_weight_windows(self):
        # -h off: no in-plane Hamming, the kz window stays
        weight = model.mrsi_kspace_weight(self.shape, hamming=False, hamming_z=True)
        central = weight[:, :, 0][weight[:, :, 0] > 0]
        np.testing.assert_allclose(central, 1)
        self.assertLess(weight[0, 0, 1], 1)
        # A truncated grid keeps its own partitions only, nothing wraps around
        truncated = model.mrsi_kspace_weight(self.shape, acquired=(24, 24, 15))
        self.assertEqual(int(truncated.any(axis=(0, 1)).sum()), self.shape[2])
        self.assertTrue(np.all(truncated <= 1))

    def test_kspace_weight_protocol_windows(self):
        hamming = model.mrsi_kspace_weight(self.shape)
        hann = model.mrsi_kspace_weight(self.shape, window="hann")
        none = model.mrsi_kspace_weight(self.shape, window="none")
        # Hann in-plane, Hamming along kz as MRSI.jl's kernel
        r = np.abs(np.fft.fftfreq(20)[3])
        self.assertAlmostEqual(hann[3, 0, 0], 0.5 + 0.5 * np.cos(2 * np.pi * r))
        self.assertAlmostEqual(hann[0, 0, 1], hamming[0, 0, 1])
        # none: neither window
        np.testing.assert_allclose(none[none > 0], 1)
        np.testing.assert_array_equal(none > 0, hamming > 0)
        with self.assertRaises(ValueError):
            model.mrsi_kspace_weight(self.shape, window="tukey")

    def test_tissue_on_grid_segmentation_fov(self):
        # GM everywhere, but segmented only for x >= 0: the rest counts as no brain
        t1 = np.diag([1.0, 1.0, 1.0, 1.0])
        t1[:3, 3] = [0, -60, -90]
        gm = np.ones((60, 120, 180))
        pve = np.stack([gm, 0 * gm, 0 * gm])
        weight = model.mrsi_kspace_weight(self.shape)
        amounts, brain, excited = model.tissue_on_grid(
            pve, t1, self.csi, self.shape, weight, lambda z: np.ones_like(z)
        )
        np.testing.assert_allclose(excited, 1, atol=1e-9)
        self.assertGreater(brain[15, 10, 5], 0.98)  # x = 22 mm
        self.assertLess(brain[3, 10, 5], 0.02)  # x = -31 mm
        np.testing.assert_allclose(amounts[0], brain, atol=1e-9)
        # Without any window the PSF rings at the edge; the amounts stay non-negative
        flat = model.mrsi_kspace_weight(self.shape, hamming=False, hamming_z=False)
        amounts, _, _ = model.tissue_on_grid(
            pve, t1, self.csi, self.shape, flat, lambda z: np.ones_like(z)
        )
        self.assertGreaterEqual(amounts.min(), 0)

    def test_uniform_and_edge(self):
        uniform = np.ones((120, 120, 180))
        half = (self.x >= 0).astype(float)
        out = model.psf_average(np.stack([uniform, half]), self.fine, self.csi, self.shape)
        np.testing.assert_allclose(out[0], 1, atol=1e-9)
        # Voxel 10 is centred at x = 0, the edge of the half space; the fine sample on the
        # edge counts to the inside, so a little above 0.5
        self.assertAlmostEqual(out[1, 10, 10, 5], 0.5, delta=0.1)
        self.assertLess(out[1, 6, 10, 5], 0.02)
        self.assertGreater(out[1, 14, 10, 5], 0.98)

    def test_voi_box(self):
        p = dict(
            slab_thickness_mm=140.0,
            voi_thickness_mm=90.0,
            slab_position=[0, 0, 0],
            slab_normal=[0, 0, 1],
            voi_position=[0, 0, 0],
        )
        offsets = model.partition_offsets_mm(p, 11)
        np.testing.assert_allclose(offsets, (np.arange(11) - 5) * 140 / 11)
        inner = model.inner_partitions(p, 11, 140 / 11)
        np.testing.assert_array_equal(np.flatnonzero(inner), [3, 4, 5, 6, 7])
        moved = dict(p, voi_position=[0, 0, 140 / 11])
        np.testing.assert_allclose(model.partition_offsets_mm(moved, 11), offsets - 140 / 11)


class Headers(unittest.TestCase):
    def test_twix_layout(self):
        """A file in the VD layout: raid header, two measurements, MeasYaps of the last one read."""
        text = "\n".join(
            [
                "### ASCCONV BEGIN ###",
                "alTR[0]\t = 450000",
                "adFlipAngleDegree[0]\t = 50.0",
                "alTE[0]\t = 1000",
                "sTXSPEC.asNucleusInfo[0].flReferenceAmplitude\t = 263.0",
                "sTXSPEC.asNucleusInfo[0].lFrequency\t = 123246948",
                "sProtConsistencyInfo.flNominalB0\t = 2.89362001419",
                "sSliceArray.asSlice[0].dThickness\t = 140.0",
                "sSpecPara.sVoI.dThickness\t = 90.0",
                "sKSpace.lPartitions\t = 31",
                "sSliceArray.asSlice[0].sNormal.dTra\t = 1.0",
                "### ASCCONV END ###",
                "",
            ]
        ).encode("latin1")

        def measurement(yaps):
            buffers = b"MeasYaps\x00" + struct.pack("<I", len(yaps)) + yaps
            return struct.pack("<2I", 8 + len(buffers), 1) + buffers

        first = measurement(text.replace(b"450000", b"999"))
        last = measurement(text)
        head = 8 + 2 * 152
        entries = b""
        for offset, block in ((head, first), (head + len(first), last)):
            entries += struct.pack("<2I2Q", 1, 1, offset, len(block)) + b"\x00" * 128
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "test.dat")
            with open(path, "wb") as f:
                f.write(struct.pack("<2I", 0, 2) + entries + first + last)
            p = model.read_protocol(path)
        self.assertAlmostEqual(p["tr_s"], 0.45)
        self.assertEqual((p["flip_deg"], p["te_s"], p["partitions"]), (50.0, 1e-3, 31))
        self.assertEqual(p["slab_position"], [0.0, 0.0, 0.0])
        self.assertEqual(p["slab_normal"], [0.0, 0.0, 1.0])

    @unittest.skipUnless(
        os.path.isfile(os.path.join(MRSI_3T_HEADER_DIR, MET_DAT)), "MRSI_3T_HEADER_DIR not set"
    )
    def test_real_3t(self):
        met = model.read_protocol(os.path.join(MRSI_3T_HEADER_DIR, MET_DAT))
        water = model.read_protocol(os.path.join(MRSI_3T_HEADER_DIR, WATER_DAT))
        for p, partitions, base in ((met, 31, 50), (water, 23, 36)):
            self.assertEqual(p["matrix"], [base, base, partitions])
            self.assertEqual(p["fov_mm"], [220.0, 220.0])
            self.assertAlmostEqual(p["tr_s"], 0.45)
            self.assertEqual(p["flip_deg"], 50.0)
            self.assertAlmostEqual(p["te_s"], 1e-3)
            self.assertEqual(p["reference_amplitude_v"], 263.0)
            self.assertEqual(
                (p["slab_thickness_mm"], p["voi_thickness_mm"], p["partitions"]),
                (140.0, 90.0, partitions),
            )
            self.assertEqual(model.field_of(p["b0_t"]), "3T")
        np.testing.assert_allclose(met["slab_normal"], [0.039831, 0.479084, 0.876865])
        self.assertEqual(met["voi_position"], met["slab_position"])

    @unittest.skipUnless(
        os.path.isfile(os.path.join(MRSI_3T_HEADER_DIR, "B1-Map-Siemens_18_MR", "1.dcm")),
        "MRSI_3T_HEADER_DIR not set",
    )
    def test_real_3t_b1_map(self):
        b1, affine, record = model.siemens_b1_map(
            os.path.join(MRSI_3T_HEADER_DIR, "B1-Map-Siemens_18_MR", "1.dcm")
        )
        self.assertEqual(b1.shape, (64, 64, 30))
        self.assertEqual(record["nominal_flip_deg"], 80.0)
        self.assertAlmostEqual(record["reference_amplitude_v"], 269.409820557)
        np.testing.assert_allclose(np.linalg.norm(affine[:3, :3], axis=0), [4, 4, 6], atol=0.01)
        inside = b1[(b1 > 0.3) & (b1 < 2)]
        self.assertTrue(0.85 < np.median(inside) < 1.05)


class DicomRescale(unittest.TestCase):
    def test_per_frame_shared_and_top_level(self):
        from pydicom.dataset import Dataset

        def group(slope, intercept):
            g, t = Dataset(), Dataset()
            t.RescaleSlope, t.RescaleIntercept = slope, intercept
            g.PixelValueTransformationSequence = [t]
            return g

        ds = Dataset()
        ds.PerFrameFunctionalGroupsSequence = [group(1, 0), group(2, -5)]
        self.assertEqual(model.frame_rescale(ds, 1), (2.0, -5.0))
        ds = Dataset()
        ds.PerFrameFunctionalGroupsSequence = [Dataset(), Dataset()]
        ds.SharedFunctionalGroupsSequence = [group(3, 1)]
        self.assertEqual(model.frame_rescale(ds, 1), (3.0, 1.0))
        ds = Dataset()
        ds.PerFrameFunctionalGroupsSequence = [Dataset()]
        self.assertEqual(model.frame_rescale(ds, 0), (1.0, 0.0))


class Names(unittest.TestCase):
    def test_groups(self):
        for name, group in (
            ("NAA+NAAG", "tNAA"),
            ("tNAA", "tNAA"),
            ("Cr+PCr", "tCr"),
            ("GPC+PCh", "tCho"),
            ("Ins", "mIns"),
            ("mIns", "mIns"),
            ("Glu", "Glu"),
            ("Glu+Gln", None),
            ("GSH", None),
        ):
            self.assertEqual(model.group_of(name), group, name)

    def test_every_constant_has_a_source(self):
        def walk(node):
            if "source" in node:
                yield node
            else:
                for entry in node.values():
                    yield from walk(entry)

        for c in list(walk(model.CONSTANTS)) + list(walk(model.WATER_CONTENT)):
            for key in [c["source"], *c["spread_sources"]]:
                self.assertIn(key, model.SOURCES)
            self.assertLessEqual(c["low"], c["value"])
            self.assertLessEqual(c["value"], c["high"])


if __name__ == "__main__":
    unittest.main()
