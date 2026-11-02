#!/usr/bin/env python3
"""From water-scaled institutional units to mmol/kg tissue water (molal), for -w W3.

The form is Gasparovic et al. 2006 with CSF correction, as written in the consensus of
Near et al. 2021 (eq. 4 and 5):

    [M] = c_iu * (H2O_MOLAL / IU_WATER_MM) * sum_t f_t R_W,t / (f_GM R_M,GM + f_WM R_M,WM)

c_iu the W3 result (pure water 55556 mM, 2 protons; the metabolite's protons are in each
fitter's basis and so in C_fit), f_t the water mole fractions f_t = v_t d_t / sum v d of
the volume fractions v_t and water contents d_t, R = saturation * exp(-TE / T2) per tissue,
with the saturation of a spoiled steady state at the flip angle the voxel gets:
nominal flip * B1 map * (reference amplitude of the scan / of the B1 map). No metabolite
in CSF; GM and WM share one molal concentration per voxel.
"""

import re
import struct

import numpy as np

H2O_MOLAL = 55510.0  # mmol/kg, the inverse of water's molar mass (Near 2021)
IU_WATER_MM = 55556.0  # the water concentration the institutional units are relative to
TISSUES = ("GM", "WM", "CSF")
GROUPS = ("tNAA", "tCr", "tCho", "Glu", "mIns")

SOURCES = {
    "gasparovic2006": dict(
        citation="Gasparovic C et al. Magn Reson Med 2006;55:1219-1226",
        doi="10.1002/mrm.20901",
        location="method; form checked in the ISMRM 2006 abstract 3224 and in Near 2021 eq. 4",
        verified=False,
        note="journal full text not accessible",
    ),
    "near2021": dict(
        citation="Near J et al. NMR Biomed 2021;34:e4257",
        doi="10.1002/nbm.4257",
        location="eq. 4 (molal form), eq. 5 (water mole fractions), text: d_GM 0.78, d_WM "
        "0.65, d_CSF 0.97 'often assumed', [H2O] 55.51 mol/kg",
        verified=True,
        note="glossary gives 55.49 mol/kg; d values cited there to Kreis, Ernst, Ross 1993",
    ),
    "ernst1993": dict(
        citation="Ernst T, Kreis R, Ross BD. Absolute quantitation of water and metabolites in "
        "the human brain. I. Compartments and water. J Magn Reson B 1993;102:1-8",
        doi=None,
        location="water contents 0.78 / 0.65 / 0.97 as cited by Near 2021, Gasparovic 2006 "
        "and the LCModel manual (WCONC 35880 = 55556 x 0.646)",
        verified=False,
        note="title, volume and pages seen in a reference list; DOI not found; values verified "
        "only as quoted",
    ),
    "marques2010": dict(
        citation="Marques JP et al. NeuroImage 2010;49:1271-1281",
        doi="10.1016/j.neuroimage.2009.10.002",
        location="abstract: 3T WM 0.81+-0.03, GM 1.35+-0.05 s; 7T WM 1.15+-0.06, GM 1.92+-0.16 s",
        verified=True,
    ),
    "wansapura1999": dict(
        citation="Wansapura JP, Holland SK, Dunn RS, Ball WS. NMR relaxation times in the human "
        "brain at 3.0 tesla. J Magn Reson Imaging 1999;9:531-538",
        doi="10.1002/(SICI)1522-2586(199904)9:4<531::AID-JMRI4>3.0.CO;2-L",
        location="abstract: T1 GM 1331, WM 832 ms; T2 printed 'GM and WM 80 and 110 ms'",
        verified=True,
        note="T2 tissue order in the abstract looks swapped; used as GM 110, WM 80 ms "
        "(as Osprey does); table not seen",
    ),
    "wright2008": dict(
        citation="Wright PJ et al. MAGMA 2008;21:121-130",
        doi="10.1007/s10334-008-0104-8",
        location="abstract: 3T WM 838+-50, GM 1607+-112 ms; 7T WM 1126+-97, GM 1939+-149 ms",
        verified=True,
    ),
    "stanisz2005": dict(
        citation="Stanisz GJ et al. Magn Reson Med 2005;54:507-512",
        doi="10.1002/mrm.20605",
        location="T1 GM 1820+-114, WM 1084+-45 ms; T2 GM 99+-7, WM 69+-3 ms",
        verified=False,
        note="secondary only (Osprey source comments); used for the spread only",
    ),
    "rooney2007": dict(
        citation="Rooney WD et al. Magn Reson Med 2007;57:308-318",
        doi="10.1002/mrm.21122",
        location="abstract: CSF T1 4.3 (+-0.2) s, independent of B0; 7T GM 2130, WM 1220 ms "
        "as quoted in Marjanska 2012",
        verified=True,
        note="7T tissue values secondary (Marjanska 2012 methods)",
    ),
    "lu2005": dict(
        citation="Lu H et al. J Magn Reson Imaging 2005;22:13-22",
        doi="10.1002/jmri.20356",
        location="CSF T1 3817+-424 ms",
        verified=False,
        note="secondary only (Osprey source comment); spread only",
    ),
    "spijkerman2018": dict(
        citation="Spijkerman JM et al. T2 mapping of cerebrospinal fluid: 3 T versus 7 T. "
        "MAGMA 2018;31:415-424",
        doi="10.1007/s10334-017-0659-3",
        location="results: ventricular CSF T2 3T 2.03 s (1.73-2.16), 7T 1.05 s",
        verified=True,
    ),
    "piechnik2009": dict(
        citation="Piechnik SK et al. Magn Reson Med 2009;61:579-586",
        doi="10.1002/mrm.21897",
        location="abstract: T2,csf 503.0+-64.3 ms, T2,brain 61.0+-2 ms (biexponential, 3T)",
        verified=True,
        note="partial-volume component; spread only",
    ),
    "bartha2002": dict(
        citation="Bartha R, Michaeli S, Merkle H, et al. In vivo 1H2O T2+ measurement in the human "
        "occipital lobe at 4T and 7T by Carr-Purcell MRI: detection of microscopic "
        "susceptibility contrast. Magn Reson Med 2002;47:742-750",
        doi="10.1002/mrm.10112",
        location="7T T2(GM) 50, T2(WM) 55, T2(CSF) 141 ms as quoted in Marjanska 2012 methods",
        verified=False,
        note="secondary (Marjanska 2012); own abstract confirms GM 9% below WM",
    ),
    "mlynarik2001": dict(
        citation="Mlynarik V, Gruber S, Moser E. NMR Biomed 2001;14:325-331",
        doi="10.1002/nbm.713",
        location="ISMRM 2001 abstract 1662 tables 1 (T1) and 2 (T2), 3T occipital GM / WM; "
        "journal abstract T1 ranges GM 0.97-1.47, WM 0.87-1.35 s match",
        verified=True,
    ),
    "li2008": dict(
        citation="Li Y et al. Comparison of T1 and T2 metabolite relaxation times in glioma "
        "and normal brain at 3T. J Magn Reson Imaging 2008;28:342-350",
        doi="10.1002/jmri.21453",
        location="normal volunteers (WM): T1 NAA 1.38+-0.13, Cr 1.38+-0.13, Cho 1.06+-0.11 s",
        verified=True,
        note="spread only (PMC3033227)",
    ),
    "ganji2012": dict(
        citation="Ganji SK et al. T2 measurement of J-coupled metabolites in the human brain "
        "at 3T. NMR Biomed 2012;25:523-529",
        doi="10.1002/nbm.1767",
        location="table 1 (GM-rich / WM-rich voxels): Glu 181/180, mIns 197/196 ms (T2+)",
        verified=True,
        note="PMC3852663",
    ),
    "kirov2008": dict(
        citation="Kirov II et al. The age dependence of regional proton metabolites T2 relaxation "
        "times in the human brain at 3 Tesla. Magn Reson Med 2008;60:790-795",
        doi="10.1002/mrm.21715",
        location="table 1, young adults: NAA 269/374, Cr 156/179, Cho 226/238 ms (GM/WM)",
        verified=True,
        note="spread only (PMC2631566)",
    ),
    "xin2013": dict(
        citation="Xin L, Schaller B, Mlynarik V, Lu H, Gruetter R. Magn Reson Med 2013;69:931-936",
        doi="10.1002/mrm.24352",
        location="table 2 (GM-rich / WM-rich voxels)",
        verified=True,
    ),
    "an2017": dict(
        citation="An L, Li S, Shen J. Simultaneous determination of metabolite concentrations, "
        "T1 and T2 relaxation times. Magn Reson Med 2017;78:2072-2081",
        doi="10.1002/mrm.26612",
        location="table 1, 7T frontal GM-rich (66%) / WM-rich (81%) voxels",
        verified=True,
        note="PMC5545076",
    ),
    "marjanska2012": dict(
        citation="Marjanska M et al. Localized 1H NMR spectroscopy in different regions of "
        "human brain in vivo at 7 T: T2 relaxation times and concentrations of cerebral "
        "metabolites. NMR Biomed 2012;25:332-339",
        doi="10.1002/nbm.1754",
        location="table 1, occipital / motor cortex voxels (not tissue-separated)",
        verified=True,
    ),
}


def constant(value, low, high, source, spread_sources=()):
    return dict(value=value, low=low, high=high, source=source, spread_sources=list(spread_sources))


# Times in s. low/high bound the spread the uncertainty budget varies over.
CONSTANTS = {
    "3T": {
        "water_T1": {
            "GM": constant(
                1.35, 1.331, 1.82, "marques2010", ("wansapura1999", "wright2008", "stanisz2005")
            ),
            "WM": constant(
                0.81, 0.81, 1.084, "marques2010", ("wansapura1999", "wright2008", "stanisz2005")
            ),
            "CSF": constant(4.3, 3.817, 4.5, "rooney2007", ("lu2005",)),
        },
        "water_T2": {
            "GM": constant(0.110, 0.099, 0.110, "wansapura1999", ("stanisz2005",)),
            "WM": constant(0.080, 0.069, 0.080, "wansapura1999", ("stanisz2005",)),
            "CSF": constant(2.03, 0.503, 2.16, "spijkerman2018", ("piechnik2009",)),
        },
        "T1": {
            "tNAA": {
                "GM": constant(1.47, 1.35, 1.59, "mlynarik2001", ("li2008",)),
                "WM": constant(1.35, 1.23, 1.47, "mlynarik2001", ("li2008",)),
            },
            "tCr": {
                "GM": constant(1.46, 1.39, 1.53, "mlynarik2001", ("li2008",)),
                "WM": constant(1.24, 1.14, 1.38, "mlynarik2001", ("li2008",)),
            },
            "tCho": {
                "GM": constant(1.30, 1.24, 1.36, "mlynarik2001", ("li2008",)),
                "WM": constant(1.08, 1.02, 1.14, "mlynarik2001", ("li2008",)),
            },
            "Glu": {
                "GM": constant(1.27, 1.17, 1.37, "mlynarik2001"),
                "WM": constant(1.17, 1.09, 1.25, "mlynarik2001"),
            },
            "mIns": {
                "GM": constant(1.23, 1.14, 1.32, "mlynarik2001"),
                "WM": constant(1.01, 0.92, 1.10, "mlynarik2001"),
            },
        },
        "T2": {
            "tNAA": {
                "GM": constant(0.247, 0.247, 0.269, "mlynarik2001", ("kirov2008",)),
                "WM": constant(0.295, 0.295, 0.374, "mlynarik2001", ("kirov2008",)),
            },
            "tCr": {
                "GM": constant(0.152, 0.152, 0.156, "mlynarik2001", ("kirov2008",)),
                "WM": constant(0.156, 0.156, 0.179, "mlynarik2001", ("kirov2008",)),
            },
            "tCho": {
                "GM": constant(0.207, 0.207, 0.226, "mlynarik2001", ("kirov2008",)),
                "WM": constant(0.187, 0.187, 0.238, "mlynarik2001", ("kirov2008",)),
            },
            "Glu": {
                "GM": constant(0.181, 0.165, 0.197, "ganji2012"),
                "WM": constant(0.180, 0.168, 0.192, "ganji2012"),
            },
            "mIns": {
                "GM": constant(0.197, 0.183, 0.211, "ganji2012"),
                "WM": constant(0.196, 0.179, 0.213, "ganji2012"),
            },
        },
    },
    "7T": {
        "water_T1": {
            "GM": constant(1.92, 1.92, 2.13, "marques2010", ("wright2008", "rooney2007")),
            "WM": constant(1.15, 1.126, 1.22, "marques2010", ("wright2008", "rooney2007")),
            "CSF": constant(4.3, 4.1, 4.5, "rooney2007"),
        },
        "water_T2": {
            "GM": constant(0.050, 0.036, 0.050, "bartha2002", ("marjanska2012",)),
            "WM": constant(0.055, 0.036, 0.055, "bartha2002", ("marjanska2012",)),
            "CSF": constant(1.05, 0.141, 1.05, "spijkerman2018", ("bartha2002",)),
        },
        "T1": {
            "tNAA": {
                "GM": constant(1.83, 1.43, 1.88, "xin2013", ("an2017",)),
                "WM": constant(1.90, 1.54, 1.96, "xin2013", ("an2017",)),
            },
            "tCr": {
                "GM": constant(1.74, 1.31, 1.80, "xin2013", ("an2017",)),
                "WM": constant(1.78, 1.53, 1.82, "xin2013", ("an2017",)),
            },
            "tCho": {
                "GM": constant(1.51, 1.17, 1.62, "xin2013", ("an2017",)),
                "WM": constant(1.32, 1.19, 1.38, "xin2013", ("an2017",)),
            },
            "Glu": {
                "GM": constant(1.61, 1.27, 1.70, "xin2013", ("an2017",)),
                "WM": constant(1.75, 1.55, 1.79, "xin2013", ("an2017",)),
            },
            "mIns": {
                "GM": constant(1.28, 1.20, 1.36, "xin2013"),
                "WM": constant(1.19, 1.12, 1.26, "xin2013"),
            },
        },
        "T2": {
            "tNAA": {
                "GM": constant(0.203, 0.132, 0.220, "an2017", ("marjanska2012",)),
                "WM": constant(0.169, 0.153, 0.185, "an2017", ("marjanska2012",)),
            },
            "tCr": {
                "GM": constant(0.127, 0.095, 0.133, "an2017", ("marjanska2012",)),
                "WM": constant(0.109, 0.102, 0.116, "an2017", ("marjanska2012",)),
            },
            "tCho": {
                "GM": constant(0.197, 0.152, 0.208, "an2017", ("marjanska2012",)),
                "WM": constant(0.130, 0.116, 0.144, "an2017", ("marjanska2012",)),
            },
            "Glu": {
                "GM": constant(0.184, 0.093, 0.202, "an2017", ("marjanska2012",)),
                "WM": constant(0.139, 0.098, 0.147, "an2017", ("marjanska2012",)),
            },
            "mIns": {
                "GM": constant(0.095, 0.093, 0.097, "marjanska2012"),
                "WM": constant(0.100, 0.085, 0.115, "marjanska2012"),
            },
        },
    },
}
WATER_CONTENT = {
    "GM": constant(0.78, 0.75, 0.81, "near2021", ("ernst1993",)),
    "WM": constant(0.65, 0.62, 0.68, "near2021", ("ernst1993",)),
    "CSF": constant(0.97, 0.94, 1.0, "near2021", ("ernst1993",)),
}
# No second verified source for the water contents: their spread is an assumed +-0.03.


def field_of(b0_tesla):
    """The CONSTANTS key of a nominal field, or an error for a field without a table."""
    for key, b0 in (("3T", 2.89), ("7T", 6.98)):
        if abs(b0_tesla - b0) < 0.3:
            return key
    raise ValueError(f"no relaxation constants for a {b0_tesla:.2f} T field")


def values(table, which="value"):
    """A CONSTANTS entry tree with each constant replaced by its value, low or high."""
    if "source" in table:
        return table[which]
    return {key: values(entry, which) for key, entry in table.items()}


# Twix headers


def _twix_buffers(f):
    header_length, n_buffers = struct.unpack("<2I", f.read(8))
    buffers = {}
    for _ in range(n_buffers):
        name = b""
        while (c := f.read(1)) not in (b"\x00", b""):
            name += c
        (length,) = struct.unpack("<I", f.read(4))
        buffers[name.decode("latin1")] = f.read(length)
    return buffers


def read_ascconv(path):
    """The ASCCONV protocol text of a twix file's last measurement (VD and later, or VB)."""
    with open(path, "rb") as f:
        first, n_meas = struct.unpack("<2I", f.read(8))
        if first == 0 and 0 < n_meas <= 64:
            offsets = []
            for _ in range(n_meas):
                _, _, offset, _ = struct.unpack("<2I2Q", f.read(24))
                f.read(128)  # patient and protocol names, not read
                offsets.append(offset)
            f.seek(offsets[-1])
        else:
            f.seek(0)
        text = _twix_buffers(f).get("MeasYaps", b"").decode("latin1")
    begin, end = text.find("### ASCCONV BEGIN"), text.find("### ASCCONV END")
    if begin < 0 or end < begin:
        raise ValueError(f"{path} has no ASCCONV protocol")
    return text[begin:end]


def ascconv_value(text, key, default=None):
    """A numeric ASCCONV entry; Siemens leaves zeros out, so a missing position is 0."""
    m = re.search(r"^" + re.escape(key) + r"\s*=\s*([-+0-9.eE]+)", text, re.M)
    if m is None:
        if default is None:
            raise KeyError(key)
        return default
    return float(m.group(1))


def read_protocol(path):
    """TR, flip, acquisition delay (TE), transmitter reference, field and VOI of a twix file."""
    text = read_ascconv(path)

    def vector(prefix):
        return [ascconv_value(text, f"{prefix}.d{a}", 0.0) for a in ("Sag", "Cor", "Tra")]

    return dict(
        file=str(path),
        tr_s=ascconv_value(text, "alTR[0]") * 1e-6,
        flip_deg=ascconv_value(text, "adFlipAngleDegree[0]"),
        te_s=ascconv_value(text, "alTE[0]") * 1e-6,
        reference_amplitude_v=ascconv_value(text, "sTXSPEC.asNucleusInfo[0].flReferenceAmplitude"),
        larmor_hz=ascconv_value(text, "sTXSPEC.asNucleusInfo[0].lFrequency"),
        b0_t=ascconv_value(text, "sProtConsistencyInfo.flNominalB0"),
        slab_thickness_mm=ascconv_value(text, "sSliceArray.asSlice[0].dThickness"),
        voi_thickness_mm=ascconv_value(text, "sSpecPara.sVoI.dThickness"),
        partitions=int(ascconv_value(text, "sKSpace.lPartitions")),
        slab_position=vector("sSliceArray.asSlice[0].sPosition"),
        slab_normal=vector("sSliceArray.asSlice[0].sNormal"),
        voi_position=vector("sSpecPara.sVoI.sPosition"),
    )


# Relaxation


def saturation(flip_rad, tr, t1):
    """Transverse signal of a spoiled steady state per unit M0: sin a (1 - E1) / (1 - cos a E1)."""
    e1 = np.exp(-tr / t1)
    return np.sin(flip_rad) * (1 - e1) / (1 - np.cos(flip_rad) * e1)


def attenuation(flip_rad, tr, te, t1, t2):
    """R of Near 2021 for a voxel's actual flip angle: saturation times exp(-TE / T2)."""
    return saturation(flip_rad, tr, t1) * np.exp(-te / t2)


def water_mole_fractions(volume, content):
    """f_t = v_t d_t / sum v d (Near 2021 eq. 5), volume (3, ...) as GM, WM, CSF."""
    weighted = np.asarray(volume, float) * np.reshape(
        [content[t] for t in TISSUES], (3,) + (1,) * (np.ndim(volume) - 1)
    )
    total = weighted.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, weighted / np.where(total > 0, total, 1), np.nan)


def molal_factor(volume, r_water, r_met, content):
    """[M] / c_iu per voxel (Near 2021 eq. 4).

    r_water (3, ...) is GM, WM, CSF; r_met (2, ...) is GM, WM.
    """
    f = water_mole_fractions(volume, content)
    numerator = (f * np.asarray(r_water)).sum(axis=0)
    denominator = f[0] * r_met[0] + f[1] * r_met[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        return H2O_MOLAL / IU_WATER_MM * numerator / denominator


def actual_flip(protocol, b1_relative, b1_reference_v):
    """Flip angle (rad) a voxel gets: the scan's nominal flip at the B1 map's relative field,
    scaled by the two transmitter reference amplitudes."""
    scale = protocol["reference_amplitude_v"] / b1_reference_v
    return np.deg2rad(protocol["flip_deg"]) * np.asarray(b1_relative) * scale


def correction_factors(volume, b1_relative, met, water, b1_reference_v, field=None, which=None):
    """{group: [M] / c_iu} maps. volume (3, ...) GM, WM, CSF fractions; met and water the
    read_protocol of the two scans. which = {(table, ..., tissue): "low"/"high"} picks a
    constant's bound instead of its value (uncertainty budget)."""
    field = field or field_of(met["b0_t"])
    table = CONSTANTS[field]

    def pick(path):
        node = {"water_content": WATER_CONTENT}.get(path[0]) or table[path[0]]
        for key in path[1:]:
            node = node[key]
        return node[(which or {}).get(path, "value")]

    a_met = actual_flip(met, b1_relative, b1_reference_v)
    a_water = actual_flip(water, b1_relative, b1_reference_v)
    content = {t: pick(("water_content", t)) for t in TISSUES}
    r_water = np.stack(
        [
            attenuation(
                a_water, water["tr_s"], water["te_s"], pick(("water_T1", t)), pick(("water_T2", t))
            )
            for t in TISSUES
        ]
    )
    factors = {}
    for group in GROUPS:
        r_met = np.stack(
            [
                attenuation(
                    a_met, met["tr_s"], met["te_s"], pick(("T1", group, t)), pick(("T2", group, t))
                )
                for t in ("GM", "WM")
            ]
        )
        factors[group] = molal_factor(volume, r_water, r_met, content)
    return factors


def group_of(name):
    """The GROUPS entry a fitter's map name belongs to: tNAA, NAA, NAA+NAAG, Ins, ... else None."""
    members = {
        "tNAA": {"NAA", "NAAG", "tNAA"},
        "tCr": {"Cr", "PCr", "tCr"},
        "tCho": {"GPC", "PCh", "Cho", "tCho"},
        "Glu": {"Glu"},
        "mIns": {"Ins", "mIns", "mI"},
    }
    parts = set(name.split("+"))
    return next((group for group, names in members.items() if parts <= names), None)


# Grids


def mrsi_kspace_weight(shape, circles_per_partition=None, hamming=True, acquired=None):
    """The k-space weighting of the CRT reconstruction on the Cartesian grid of `shape`
    (FFT order): the disk each partition's circles cover (circle i at radius (i - 1/2) / N,
    its annulus to i / N), the in-plane Hamming of the radius and the Hamming along kz, as
    MRSI.jl applies them in ICE mode. `acquired`: the matrix of a scan of the same FOV that
    was zero-filled onto `shape`, whose (smaller) weighting is then returned."""
    acquired = acquired or shape
    nx, ny, nz = acquired
    n_circles = nx // 2
    if circles_per_partition is None:
        part_max = nz // 2
        circles_per_partition = [
            (
                n_circles
                if part_max == 0
                else max(2, int(np.ceil(np.sqrt(n_circles**2 - (p * n_circles / part_max) ** 2))))
            )
            for p in range(-part_max, part_max + 1)
        ]
    # Radius in cycles per voxel of the acquired grid, on the FFT grid of `shape`
    kx = np.fft.fftfreq(shape[0], 1 / shape[0]) / nx
    ky = np.fft.fftfreq(shape[1], 1 / shape[1]) / ny
    radius = np.sqrt(kx[:, None] ** 2 + ky[None, :] ** 2)
    in_plane = 0.54 + 0.46 * np.cos(2 * np.pi * radius) if hamming else np.ones_like(radius)
    z_window = (
        0.54 + 0.46 * np.cos(2 * np.pi * np.linspace(-0.5, 0.5, nz))
        if hamming and nz > 1
        else np.ones(nz)
    )
    weight = np.zeros(shape)
    for index, p in enumerate(range(-(nz // 2), nz - nz // 2)):
        limit = circles_per_partition[index] / nx
        weight[:, :, p % shape[2]] = np.where(radius < limit, in_plane, 0) * z_window[index]
    return weight / weight[0, 0, 0]


def psf_average(
    volumes, affine, csi_affine, csi_shape, oversample=(4, 4, 4), weight=None, profile=None
):
    """volumes (k, ...) on the grid of `affine`, as the MRSI sees them: sampled on a grid q
    times finer than the CSI's, multiplied by `profile` (a function of the fine grid's
    partition coordinate in CSI voxels), filtered with the reconstruction's k-space weight
    and read at the CSI voxel centres."""
    from scipy.ndimage import map_coordinates

    q = np.asarray(oversample)
    fine_shape = tuple(int(n) for n in np.asarray(csi_shape) * q)
    grids = np.meshgrid(*[np.arange(n) / qq for n, qq in zip(fine_shape, q)], indexing="ij")
    csi_index = np.stack([g.ravel() for g in grids] + [np.ones(grids[0].size)])
    source_index = (np.linalg.inv(affine) @ csi_affine @ csi_index)[:3]
    weight = mrsi_kspace_weight(csi_shape) if weight is None else weight
    # The CSI's k-space weight placed in the fine grid's spectrum (both in cycles per CSI voxel)
    fine_weight = np.zeros(fine_shape)
    index = [np.round(np.fft.fftfreq(n, 1 / n)).astype(int) for n in csi_shape]
    fine_weight[np.ix_(*[i % n for i, n in zip(index, fine_shape)])] = weight
    excitation = 1.0 if profile is None else profile(grids[2])
    out = []
    for volume in volumes:
        fine = map_coordinates(
            np.asarray(volume, float), source_index, order=1, mode="constant"
        ).reshape(fine_shape)
        filtered = np.fft.ifftn(np.fft.fftn(fine * excitation) * fine_weight).real
        out.append(filtered[:: q[0], :: q[1], :: q[2]])
    return np.stack(out)


def partition_offsets_mm(protocol, n_partitions):
    """Distance of each partition centre from the VOI centre along the slab normal (mm)."""
    thickness = protocol["slab_thickness_mm"] / n_partitions
    normal = np.asarray(protocol["slab_normal"], float)
    normal = normal / np.linalg.norm(normal)
    shift = float(np.dot(np.subtract(protocol["slab_position"], protocol["voi_position"]), normal))
    return (np.arange(n_partitions) - (n_partitions - 1) / 2) * thickness + shift


def voi_profile(protocol, n_partitions):
    """The excited VOI as a box along the partition axis, a function of the CSI partition
    coordinate (voxels), for psf_average."""
    thickness = protocol["slab_thickness_mm"] / n_partitions
    first = partition_offsets_mm(protocol, n_partitions)[0]
    half = protocol["voi_thickness_mm"] / 2

    def profile(z_voxels):
        return (np.abs(first + z_voxels * thickness) <= half).astype(float)

    return profile


def inner_partitions(protocol, n_partitions, margin_mm):
    """Partitions whose centre lies at least margin_mm inside the VOI edge."""
    offsets = partition_offsets_mm(protocol, n_partitions)
    return np.abs(offsets) <= protocol["voi_thickness_mm"] / 2 - margin_mm


# B1 maps


def siemens_b1_map(path):
    """Relative B1 of a Siemens tfl_b1map flip angle map (enhanced DICOM): the pixel is ten
    times the flip angle the preparation reached, its nominal angle the FlipAngle of the
    series; also its RAS affine and transmitter reference amplitude."""
    import pydicom

    ds = pydicom.dcmread(path)
    if "FLIP ANGLE MAP" not in [str(v) for v in ds.ImageType]:
        raise ValueError(f"{path} is not a flip angle map (ImageType {list(ds.ImageType)})")
    shared = ds.SharedFunctionalGroupsSequence[0]
    nominal = float(shared.MRTimingAndRelatedParametersSequence[0].FlipAngle)
    frames = ds.PerFrameFunctionalGroupsSequence
    orientation = np.asarray(
        (frames[0].get("PlaneOrientationSequence") or shared.PlaneOrientationSequence)[
            0
        ].ImageOrientationPatient,
        float,
    )
    measures = (frames[0].get("PixelMeasuresSequence") or shared.PixelMeasuresSequence)[0]
    row_spacing, col_spacing = (float(v) for v in measures.PixelSpacing)
    positions = np.array([f.PlanePositionSequence[0].ImagePositionPatient for f in frames], float)
    pixels = ds.pixel_array.astype(float)  # frames, rows, columns
    slope = float(frames[0].PixelValueTransformationSequence[0].RescaleSlope)
    intercept = float(frames[0].PixelValueTransformationSequence[0].RescaleIntercept)
    pixels = pixels * slope + intercept
    normal = np.cross(orientation[:3], orientation[3:])
    order = np.argsort(positions @ normal)
    positions, pixels = positions[order], pixels[order]
    step = (positions[-1] - positions[0]) / max(len(positions) - 1, 1)
    # LPS affine of (column, row, frame) indices, then to RAS
    lps = np.eye(4)
    lps[:3, 0] = orientation[:3] * col_spacing
    lps[:3, 1] = orientation[3:] * row_spacing
    lps[:3, 2] = step
    lps[:3, 3] = positions[0]
    affine = np.diag([-1, -1, 1, 1]) @ lps
    with open(path, "rb") as f:
        text = f.read().decode("latin1")
    m = re.search(r"sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude\s*=\s*([-+0-9.eE]+)", text)
    if m is None:
        raise ValueError(f"{path} carries no transmitter reference amplitude")
    volume = pixels.transpose(2, 1, 0) / (10.0 * nominal)  # (column, row, frame)
    return (
        volume,
        affine,
        dict(
            file=str(path),
            nominal_flip_deg=nominal,
            reference_amplitude_v=float(m.group(1)),
            convention="relative B1 = pixel / (10 x FlipAngle of the series)",
        ),
    )


def b1_on_grid(b1, affine, csi_affine, csi_shape, valid=(0.3, 2.0), sigma_vox=1.0):
    """The relative B1 at the CSI voxel centres: values outside `valid` (noise, no signal)
    left out and filled from their neighbours by normalised Gaussian smoothing."""
    from scipy.ndimage import gaussian_filter, map_coordinates

    b1 = np.asarray(b1, float)
    good = (b1 > valid[0]) & (b1 < valid[1])
    w = gaussian_filter(good.astype(float), sigma_vox)
    smooth = gaussian_filter(np.where(good, b1, 0), sigma_vox) / np.maximum(w, 1e-6)
    filled = np.where(good, b1, np.where(w > 0.05, smooth, 0))
    grid = np.meshgrid(*[np.arange(n) for n in csi_shape], indexing="ij")
    index = np.stack([g.ravel() for g in grid] + [np.ones(grid[0].size)])
    source = (np.linalg.inv(affine) @ csi_affine @ index)[:3]
    on_grid = map_coordinates(filled, source, order=1, mode="constant", cval=0.0).reshape(csi_shape)
    support = map_coordinates(good.astype(float), source, order=1, mode="constant").reshape(
        csi_shape
    )
    return on_grid, support
