#!/usr/bin/env python3
"""C_fit of water_scaling.py: synthetic spectra of known concentration through a fitter.

The synthetic metabolite FID is sum_i c_i * basis_i, built from the fitter's own basis on
the grid of the data; its water is 2 * 55556 * a_p, a_p the basis' per-proton amplitude
(the Cr CH3 singlet at the first stored point; for a deep basis made from an LCModel basis,
that basis' a_p times the conversion's scale). With F the fitter's result and A_w the
water amplitude from the estimator used in vivo, C_fit = c * A_w / F. PHIVE's synthetic
sits at its in vivo operating point: mix, line shape and the dataset's frequency offset.
"""

import glob
import hashlib
import json
import os
import re
import subprocess

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

WATER_MM = 55556.0
CR_CH3_PPM, CR_CH3_PROTONS = 3.027, 3
# A brain-like mix; a fitter's synthetic keeps what its basis has and its control fits.
TRUTH = {
    "NAA": 10.0,
    "NAAG": 1.5,
    "Cr": 4.0,
    "PCr": 4.0,
    "GPC": 1.0,
    "PCh": 1.0,
    "Ins": 6.0,
    "Glu": 8.0,
    "Gln": 3.0,
    "GSH": 2.0,
    "GABA": 1.5,
    "Tau": 1.5,
    "Asp": 2.0,
}
GROUPS = {
    "tNAA": ("NAA", "NAAG"),
    "tCr": ("Cr", "PCr"),
    "tCho": ("GPC", "PCh", "Cho"),
    "Glu": ("Glu",),
    "mIns": ("Ins",),
}
# The deep fitters' names for the groups
DEEP_NAMES = {"tNAA": "tNAA", "tCr": "tCr", "tCho": "tCho", "Glu": "Glu", "mIns": "Ins"}
# Line widths in ppm (3 Hz metabolites, 5 Hz Lorentzian and 3.7 Hz Gaussian water at 3T),
# noise relative to the metabolite FID's first point. LCModel's result depends on the Lorentzian
# broadening beyond its basis' own (its lineshape model misses Lorentzian tails, not Gaussian
# ones); 3 Hz matches its in vivo fits on average.
METABOLITE_LW_PPM, WATER_LW_PPM, WATER_GAUSS_PPM = 0.024, 0.04, 0.03
NOISE = 2e-4
DEEP_GRID = (8, 8, 3)
# PHIVE's gain depends on where its input sits in its training distribution, so its synthetic
# carries its median in vivo mix at 7T (basis units per mM, tCr 8) and line shape (Lorentzian,
# Gaussian FWHM in ppm), shifted by the offset of the dataset's NAA and Cr singlets.
PHIVE_MIX = {
    "NAA": 11.55,
    "NAAG": 2.125,
    "Cr": 3.921,
    "PCr": 4.079,
    "GPC": 1.023,
    "PCh": 0.937,
    "Ins": 7.26,
    "Glu": 8.113,
    "Gln": 1.923,
    "GSH": 2.098,
    "GABA": 2.852,
    "Tau": 2.144,
    "Asp": 0.964,
    "Scyllo": 1.851,
    "Glc": 0.068,
    "MM_mea": 0.01,
    "Lip_c": 1.946,
}
PHIVE_LINE_PPM = (0.027, 0.0067)
OFFSET_PEAKS = {"NAA": 2.008, "Cr": 3.027}
# A deep basis derives from an LCModel basis if |deep| / |source| is one scale times an
# envelope for all metabolites, to these tolerances, over its first TRANSFER_SECONDS.
TRANSFER_SECONDS, TRANSFER_TOLERANCE = 0.15, 1e-3


def line_amplitude(fid, dt, f_expected, tmax=0.6, max_points=4000):
    """|A| at t = 0 of the line nearest f_expected, from a fit of all lines with one envelope.

    The FID is modelled as exp(-a t - b t^2) sum_j A_j exp(2 pi i f_j t), with the lines
    taken from the peaks of its spectrum.
    """
    n_fit = min(fid.size, int(tmax / dt))
    n = np.arange(0, n_fit, max(1, n_fit // max_points))
    t, y = n * dt, fid[n]
    zf = max(fid.size, 16 * n_fit)
    # A half Hann window keeps truncation ripples out of the peak list
    spec = np.abs(np.fft.fft(fid[:n_fit] * np.hanning(2 * n_fit)[n_fit:], zf))
    peaks = np.where(
        (spec > np.roll(spec, 1)) & (spec >= np.roll(spec, -1)) & (spec > 0.02 * spec.max())
    )[0]
    freqs = np.fft.fftfreq(zf, dt)[peaks]
    alias = 1.0 / dt
    target = int(np.argmin(np.abs((freqs - f_expected + alias / 2) % alias - alias / 2)))

    def design(p):
        return np.exp(-p[0] * t - p[1] * t**2)[:, None] * np.exp(2j * np.pi * np.outer(t, p[2:]))

    def residual(p):
        m = design(p)
        r = (y - m @ np.linalg.lstsq(m, y, rcond=None)[0]) / np.abs(y).max()
        return np.concatenate([r.real, r.imag])

    sol = least_squares(residual, np.concatenate([[10.0, 0.0], freqs]), diff_step=1e-10)
    return float(np.abs(np.linalg.lstsq(design(sol.x), y, rcond=None)[0][target]))


def synthetic_water(amplitude, n, dt, larmor_hz):
    """A Voigt water line of the given amplitude at t = 0."""
    t = np.arange(n) * dt
    lw, fwhm_gauss = WATER_LW_PPM * larmor_hz * 1e-6, WATER_GAUSS_PPM * larmor_hz * 1e-6
    return amplitude * np.exp(
        -np.pi * lw * t - (np.pi * fwhm_gauss * t) ** 2 / (4 * np.log(2)) + 2j * np.pi * 5.0 * t
    )


def add_noise(fid, seed):
    rng = np.random.default_rng(seed)
    sd = NOISE * np.abs(fid[..., :1]).max()
    return fid + sd * (rng.standard_normal(fid.shape) + 1j * rng.standard_normal(fid.shape))


def truth(names, mix=None):
    """The entries of a mix (TRUTH by default) a fitter's synthetic carries."""
    return {name: c for name, c in (TRUTH if mix is None else mix).items() if name in names}


def group_truth(present):
    return {group: sum(present.get(m, 0.0) for m in members) for group, members in GROUPS.items()}


def group_of(name):
    """The GROUPS entry a map belongs to, by its own or a member's name, else None."""
    return next((g for g, members in GROUPS.items() if name == g or name in members), None)


def constants(fitted, water_amplitude, present):
    """C per group (c * A_w / F) and C_fit, their mean."""
    known = group_truth(present)
    per_group = {
        group: known[group] * water_amplitude / fitted[group]
        for group in GROUPS
        if fitted.get(group, 0) > 0 and known[group] > 0
    }
    if not per_group:
        raise RuntimeError("the fitter reported none of " + ", ".join(GROUPS))
    return float(np.mean(list(per_group.values()))), per_group


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def operating_point(fitter, offset_hz=0.0):
    """Mix, Lorentzian and Gaussian FWHM (ppm) and frequency offset (Hz) of a fitter's synthetic."""
    if fitter == "dlfit":
        lorentzian, gaussian = PHIVE_LINE_PPM
        return dict(
            mix=PHIVE_MIX, lorentzian_ppm=lorentzian, gaussian_ppm=gaussian, offset_hz=offset_hz
        )
    return dict(mix=TRUTH, lorentzian_ppm=METABOLITE_LW_PPM, gaussian_ppm=0.0, offset_hz=0.0)


def synthetic_record(present, water_amplitude, point=None):
    point = point or operating_point("lcmodel")
    return {
        "concentrations_mM": present,
        "metabolite_linewidth_ppm": point["lorentzian_ppm"],
        "metabolite_gaussian_ppm": point["gaussian_ppm"],
        "frequency_offset_hz": point["offset_hz"],
        "water_lorentzian_ppm": WATER_LW_PPM,
        "water_gaussian_ppm": WATER_GAUSS_PPM,
        "noise_relative": NOISE,
        "water_amplitude": water_amplitude,
    }


# LCModel


def read_lcmodel_basis(path):
    """(HZPPPM, BADELT, NDATAB) and {name: (CONC, ISHIFT, spectrum as stored)} of a .BASIS file."""
    text = open(path, errors="ignore").read()
    number = r"\s*=\s*([0-9.eEdD+-]+)"
    header = [
        float(re.search(key + number, text).group(1).replace("D", "E"))
        for key in ("HZPPPM", "BADELT", "NDATAB")
    ]
    entries = {}
    for part in re.split(r"\$NMUSED", text)[1:]:
        m = re.search(r"\$BASIS\s(.*?)\$END(.*)", part, re.S)
        if m is None:
            continue
        head, data = m.groups()
        name = re.search(r"METABO\s*=\s*'\s*([^']+?)\s*'", head).group(1)
        values = np.array(data.replace("D", "E").split(), dtype=float)
        entries[name] = (
            float(re.search("CONC" + number, head).group(1)),
            int(re.search(r"ISHIFT\s*=\s*(-?\d+)", head).group(1)),
            values[0::2] + 1j * values[1::2],
        )
    return header, entries


def lcmodel_a_p(header, entries):
    """a_p of a read_lcmodel_basis result: its Cr CH3 singlet per proton and mM."""
    hz_basis, dt_basis, n_basis = header
    conc, ishift, spectrum = entries["Cr"]
    df = 1.0 / (int(n_basis) * dt_basis)
    a_p = line_amplitude(
        np.fft.ifft(spectrum), dt_basis, (4.65 - CR_CH3_PPM) * hz_basis + ishift * df
    )
    return a_p / CR_CH3_PROTONS / conc


def lcmodel_synthetic(basis_path, omitted, hzpppm, deltat, nunfil, linewidth_ppm):
    """Metabolite FID on the data grid (LCModel .RAW convention), its concentrations and a_p.

    Each basis FID gets its ISHIFT, is stretched to the data's HZPPPM (ppm preserved) and
    interpolated onto the data dwell.
    """
    header, basis = read_lcmodel_basis(basis_path)
    hz_basis, dt_basis, n_basis = header
    n_basis = int(n_basis)
    df = 1.0 / (n_basis * dt_basis)
    tb = np.arange(n_basis // 2) * dt_basis
    tau = np.arange(nunfil) * deltat * hzpppm / hz_basis
    last = int(tau[-1] / dt_basis) + 8

    def on_data_grid(spectrum, ishift):
        fid = np.fft.ifft(spectrum)[: n_basis // 2] * np.exp(-2j * np.pi * ishift * df * tb)
        return CubicSpline(tb[:last], fid[:last].real)(tau) + 1j * CubicSpline(
            tb[:last], fid[:last].imag
        )(tau)

    present = truth(set(basis) - set(omitted))
    met = sum(
        c / basis[name][0] * on_data_grid(basis[name][2], basis[name][1])
        for name, c in present.items()
    )
    t = np.arange(nunfil) * deltat
    return met * np.exp(-np.pi * linewidth_ppm * hzpppm * t), present, lcmodel_a_p(header, basis)


def read_lcmodel_table(path):
    """{name: (concentration, %SD)} of an LCModel table, empty for a failed fit."""
    text = open(path, errors="ignore").read()
    if "FATAL" in text:
        return {}
    rows = re.findall(r"^\s*([-\d.Ee+]+)\s+(\d+)%\s+[-\d.Ee+]+\s+(\S+)\s*$", text, re.M)
    return {name: (float(conc), int(sd)) for conc, sd, name in rows}


def lcmodel_groups(table, present):
    """The GROUPS' concentrations, summed over the members LCModel reports."""
    return {
        group: sum(table[m][0] for m in members if m in table)
        for group, members in GROUPS.items()
        if any(m in present for m in members) and all(m in table for m in members if m in present)
    }


def fit_lcmodel_synthetic(
    control_path,
    lcmodel,
    workdir,
    water_fit,
    water_amplitude,
    seed=0,
    linewidth_ppm=METABOLITE_LW_PPM,
):
    """The synthetic spectrum fitted with a voxel's control: F per group, A_w, the
    concentrations and the inputs."""
    os.makedirs(workdir, exist_ok=True)
    control = open(control_path).read()
    grid = re.search(
        r"HZPPPM\s*=\s*([-\d.eE+]+),\s*DELTAT\s*=\s*([-\d.eE+]+),\s*NUNFIL\s*=\s*(\d+)", control
    )
    hzpppm, deltat, nunfil = float(grid.group(1)), float(grid.group(2)), int(grid.group(3))
    basis_path = re.search(r"FILBAS\s*=\s*'([^']*)'", control).group(1)
    raw_path = re.search(r"FILRAW\s*=\s*'([^']*)'", control).group(1)
    omitted = re.findall(r"CHOMIT\(\d+\)\s*=\s*'([^']*)'", control)

    met, present, a_p = lcmodel_synthetic(
        basis_path, omitted, hzpppm, deltat, nunfil, linewidth_ppm
    )
    met = add_noise(water_amplitude / (2 * WATER_MM * a_p) * met, seed)
    header = open(raw_path).read().rsplit("$END", 1)[0] + "$END\n"
    raw = os.path.join(workdir, "synthetic.RAW")
    with open(raw, "w") as f:
        f.write(header + "".join(f"{v.real:14.5e}{v.imag:14.5e}\n" for v in met))
    control = re.sub(r"^\s*(DOWS|FILH2O)\s*=.*\n", "", control, flags=re.M)
    control = re.sub(
        r"(FIL(?:TAB|PS|COO|CSV|PRI))\s*=\s*'[^']*'",
        lambda m: f"{m.group(1)}='{os.path.join(workdir, 'synthetic.' + m.group(1)[3:].lower())}'",
        control,
    )
    control = re.sub(r"FILRAW\s*=\s*'[^']*'", f"FILRAW='{raw}'", control)
    control_out = os.path.join(workdir, "synthetic.control")
    with open(control_out, "w") as f:
        f.write(control)
    with open(control_out) as stdin:
        subprocess.run(
            [lcmodel],
            stdin=stdin,
            cwd=workdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    table = read_lcmodel_table(os.path.join(workdir, "synthetic.tab"))
    if not table:
        raise RuntimeError(f"LCModel did not fit the synthetic spectrum, see {workdir}")
    estimated = water_fit(synthetic_water(water_amplitude, nunfil, deltat, hzpppm * 1e6), deltat)
    return (
        lcmodel_groups(table, present),
        estimated,
        present,
        {
            "a_p": a_p,
            "basis": basis_path,
            "basis_sha256": sha256(basis_path),
            "hzpppm": hzpppm,
            "dwell_s": deltat,
            "points": nunfil,
            "control": control_path,
            "omitted": omitted,
        },
    )


def calibrate_lcmodel(control_path, lcmodel, workdir, water_fit, water_amplitude):
    """C_fit of LCModel with the control, basis and grid of the in vivo fit."""
    fitted, estimated, present, inputs = fit_lcmodel_synthetic(
        control_path, lcmodel, workdir, water_fit, water_amplitude
    )
    c_fit, per_group = constants(fitted, estimated, present)
    return c_fit, dict(
        inputs, C_per_metabolite=per_group, synthetic=synthetic_record(present, water_amplitude)
    )


# Deep fitters


def deep_basis(fitter, info):
    """The basis FIDs as the fitter fits them, their names, the dwell time, the basis'
    reference frequency (MHz) and the file."""
    from mrs_utils.basisset import Basisset

    if fitter == "gpufit":
        from forD_gpu_fit import GPUFIT

        # The run config as JSON: its parsed form allocates tensors on the fit device.
        config_path, path = GPUFIT.BASIS_SETS[info["gpufit_basis"]]
        with open(config_path) as f:
            config = json.load(f)
        io = config["io_config"]
        names, conjugate = io["metabolite_names"], config["preprocessor_config"]["conjugate_basis"]
        dwell, larmor_mhz, nucleus = io["dwelltime"], io["reference_frequency"], io["nucleus"]
    else:
        from forD import DLFIT

        io, pre = DLFIT.config.io_config, DLFIT.config.preprocessor_config
        path, names, conjugate = str(DLFIT.basis_path), pre.basis_signal_names, pre.conjugate_basis
        dwell, larmor_mhz, nucleus = io.dwelltime, io.reference_frequency, io.nucleus
    basis = Basisset.from_matlab(
        path=path,
        metabolite_names=names,
        conjugate_basis=conjugate,
        dwelltime=dwell,
        reference_frequency=larmor_mhz,
        nucleus=nucleus,
    )
    return np.asarray(basis.fids), list(basis.metabolite_names), float(dwell), larmor_mhz, path


def source_candidates(extra=()):
    """LCModel bases a deep basis may be made from: the given ones and those shipped with forD."""
    try:
        from forD import DLFIT

        shipped = sorted(glob.glob(os.path.join(os.path.dirname(str(DLFIT.basis_path)), "*.basis")))
    except ImportError:
        shipped = []
    paths = [os.path.normpath(p) for p in (*extra, *shipped) if os.path.isfile(p)]
    return list(dict.fromkeys(paths))


def transfer_scale(fids, names, dwell, header, entries):
    """s with |deep| = s exp(-a t - b t^2) |source per mM, linearly interpolated|, the same for
    every shared metabolite; None if the deep basis is not such a copy of the source."""
    hz_basis, dt_basis, n_basis = header
    t = np.arange(fids.shape[0]) * dwell
    t = t[t < TRANSFER_SECONDS]
    t_source = np.arange(int(n_basis) // 2) * dt_basis
    scales, residuals = [], []
    for name in sorted(set(names) & set(entries)):
        conc, _, spectrum = entries[name]
        source = np.fft.ifft(spectrum)[: t_source.size] / conc
        lin = np.interp(t, t_source, source.real) + 1j * np.interp(t, t_source, source.imag)
        keep = np.abs(lin) > 0.1 * np.abs(lin).max()
        log_ratio = np.log(np.abs(fids[: t.size, names.index(name)][keep]) / np.abs(lin[keep]))
        design = np.stack([np.ones(keep.sum()), -t[keep], -t[keep] ** 2], axis=1)
        coef = np.linalg.lstsq(design, log_ratio, rcond=None)[0]
        scales.append(np.exp(coef[0]))
        residuals.append(np.std(log_ratio - design @ coef))
    if "Cr" not in entries or len(scales) < 3 or max(residuals) > TRANSFER_TOLERANCE:
        return None
    return float(np.median(scales)) if max(scales) / min(scales) - 1 < TRANSFER_TOLERANCE else None


def deep_a_p(fids, names, dwell, larmor_mhz, sources=()):
    """a_p of a deep basis and how it was found: its scale to the first source it derives from
    times that source's a_p, else the Cr CH3 fit on the basis itself (biased on short records)."""
    singlet = (
        line_amplitude(fids[:, names.index("Cr")], dwell, (CR_CH3_PPM - 4.7) * larmor_mhz)
        / CR_CH3_PROTONS
    )
    for path in sources:
        header, entries = read_lcmodel_basis(path)
        scale = transfer_scale(fids, names, dwell, header, entries)
        if scale is not None:
            source_a_p = lcmodel_a_p(header, entries)
            return scale * source_a_p, {
                "a_p_method": "source basis",
                "a_p_source_basis": path,
                "a_p_source_sha256": sha256(path),
                "a_p_source": source_a_p,
                "a_p_transfer_scale": scale,
                "a_p_singlet_fit": singlet,
            }
    return singlet, {"a_p_method": "singlet fit", "a_p_sources_tried": list(sources)}


def peak_hz(magnitude, hz, k0, half=6):
    """Frequency of each row's maximum within half points of k0, parabolically interpolated."""
    first, stop = max(k0 - half, 1), k0 + half + 1
    k = first + np.argmax(magnitude[:, first:stop], axis=1)
    k = np.clip(k, 1, magnitude.shape[1] - 2)
    y0, y1, y2 = (magnitude[np.arange(k.size), k + d] for d in (-1, 0, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.interp(k + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2), np.arange(hz.size), hz)


def frequency_offset(spectra, ppm, fids, names, dwell, larmor_hz, n_points):
    """Record of the median offset "hz" of the NAA and Cr singlets of the spectra (voxels x
    points on the ppm axis, model convention) from the same peaks of PHIVE's synthetic at 0 Hz,
    over the better half by S/N."""
    ppm = np.asarray(ppm, float)
    spectra = np.asarray(spectra).reshape(-1, ppm.size)
    spectra = spectra[np.any(spectra != 0, axis=1)]
    mhz = larmor_hz * 1e-6
    hz = (ppm - 4.7) * mhz
    magnitude = np.abs(spectra)

    def nearest(f):
        return int(np.argmin(np.abs(hz - f)))

    # S/N: the NAA peak over the spread of the 40 highest-ppm points
    noise = np.std(spectra[:, np.argsort(ppm)[-40:]].real, axis=1)
    snr = magnitude[:, nearest((OFFSET_PEAKS["NAA"] - 4.7) * mhz)] / np.maximum(noise, 1e-30)
    good = snr > np.median(snr)
    if good.sum() < 10:
        raise ValueError(f"{good.sum()} spectra are too few to measure the frequency offset")
    t = np.arange(n_points) * dwell
    lorentzian, gaussian = (w * mhz for w in PHIVE_LINE_PPM)
    shape = np.exp(-np.pi * lorentzian * t - (np.pi * gaussian * t) ** 2 / (4 * np.log(2)))
    # Overlapping lines (NAAG, PCr, ...) move the peaks; the synthetic carries them too.
    mix = truth(names, PHIVE_MIX)
    fid = sum(c * fids[:n_points, names.index(name)] for name, c in mix.items()) * shape
    reference = np.abs(np.exp(-2j * np.pi * np.outer(hz, t)) @ fid)[None]
    offsets = {}
    for name, ppm0 in OFFSET_PEAKS.items():
        f_basis = peak_hz(reference, hz, nearest((ppm0 - 4.7) * mhz))[0]
        offsets[name] = peak_hz(magnitude[good], hz, nearest(f_basis)) - f_basis
    record = {name + "_hz": float(np.nanmedian(v)) for name, v in offsets.items()}
    offset = float(np.nanmedian(np.concatenate(list(offsets.values()))))
    return dict(hz=offset, voxels=int(good.sum()), **record)


def window_ppm(grid, n):
    """ppm of a deep fitter's input window, as run_deepmrsi.py writes it."""
    parent = int(grid["signal_length"])
    step = 1.0 / (grid["dwelltime_s"] * parent * grid["reference_frequency_mhz"])
    return 4.7 + (int(grid["interval_bounds"][0]) + np.arange(n) - (parent >> 1)) * step


def fit_deep_synthetic(
    fitter, info, water_fit, water_amplitude, n_points, point, seed=0, sources=()
):
    """The synthetic spectrum at an operating point through PHIVE (dlfit) or SpatialRegu
    (gpufit) on a small grid of identical voxels: F per group, A_w, the concentrations, the
    inputs and the fitter's input spectra with their ppm axis."""
    import deep_crt_mrsi.deepmrsi as fire

    info = dict(info, b0_corrected=True)
    if fitter == "dlfit":
        fire._check_field_matches_basis(info, True)
    fit = (fire._dlfit_fitter if fitter == "dlfit" else fire._gpufit_fitter)(info, False)
    fids, names, dt, larmor_mhz, basis_path = deep_basis(fitter, info)
    a_p, a_p_record = deep_a_p(fids, names, dt, larmor_mhz, source_candidates(sources))
    present = truth(names, point["mix"])
    met = np.zeros(n_points, complex)
    for name, c in present.items():
        column = fids[:n_points, names.index(name)]
        met[: column.size] += c * column
    larmor = float(info["larmor_frequency"])
    t = np.arange(n_points) * dt
    lorentzian, gaussian = (point[k] * larmor * 1e-6 for k in ("lorentzian_ppm", "gaussian_ppm"))
    met *= np.exp(
        -np.pi * lorentzian * t
        - (np.pi * gaussian * t) ** 2 / (4 * np.log(2))
        + 2j * np.pi * point["offset_hz"] * t
    )
    met *= water_amplitude / (2 * WATER_MM * a_p)
    # The fitters conjugate their input, so the data carry the conjugate of the model.
    data = add_noise(np.broadcast_to(np.conj(met), DEEP_GRID + (n_points,)), seed).astype(
        np.complex64
    )
    mask = np.ones(DEEP_GRID, dtype=bool)
    images, out_names, _, _, windowed, grid = fit(data, mask)
    out_names = list(out_names)
    fitted = {
        group: float(np.median(np.asarray(images[out_names.index(name)])[mask]))
        for group, name in DEEP_NAMES.items()
        if name in out_names
    }
    estimated = water_fit(synthetic_water(water_amplitude, n_points, dt, larmor), dt)
    windowed = np.asarray(windowed).reshape(-1, np.shape(windowed)[-1])
    inputs = {
        "a_p": a_p,
        **a_p_record,
        "basis": os.path.normpath(basis_path),
        "basis_sha256": sha256(basis_path),
        "larmor_hz": larmor,
        "dwell_s": dt,
        "points": n_points,
        "fitter": fitter,
    }
    return fitted, estimated, present, inputs, (windowed, window_ppm(grid, windowed.shape[1]))


def calibrate_deep(fitter, info, water_fit, water_amplitude, n_points, invivo=None, sources=()):
    """C_fit of PHIVE (dlfit) or SpatialRegu (gpufit) for the data's field and grid. PHIVE's is
    per group, at the frequency offset of invivo = (its input spectra, their ppm axis)."""
    offset = None
    if fitter == "dlfit":
        if invivo is None:
            raise ValueError("PHIVE's calibration needs the spectra it fitted")
        fids, names, dt, _, _ = deep_basis(fitter, info)
        offset = frequency_offset(*invivo, fids, names, dt, info["larmor_frequency"], n_points)
    point = operating_point(fitter, offset["hz"] if offset else 0.0)
    fitted, estimated, present, inputs, _ = fit_deep_synthetic(
        fitter, info, water_fit, water_amplitude, n_points, point, sources=sources
    )
    c_fit, per_group = constants(fitted, estimated, present)
    record = dict(
        inputs,
        C_per_metabolite=per_group,
        per_group=fitter == "dlfit",
        synthetic=synthetic_record(present, water_amplitude, point),
    )
    if offset:
        record["frequency_offset"] = offset
    return c_fit, record
