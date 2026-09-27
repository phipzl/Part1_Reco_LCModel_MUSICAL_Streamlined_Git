#!/usr/bin/env python3
"""C_fit of water_scaling.py: synthetic spectra of known concentration through a fitter.

The synthetic metabolite FID is sum_i c_i * basis_i, built from the fitter's own basis on
the grid of the data; its water is 2 * 55556 * a_p, a_p the basis' per-proton amplitude
(the Cr CH3 singlet at the first stored point). With F the fitter's result and A_w the
water amplitude from the estimator used in vivo, C_fit = c * A_w / F.
"""

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


def truth(names):
    """The TRUTH entries a fitter's synthetic carries."""
    return {name: c for name, c in TRUTH.items() if name in names}


def group_truth(present):
    return {group: sum(present.get(m, 0.0) for m in members) for group, members in GROUPS.items()}


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


def synthetic_record(present, water_amplitude):
    return {
        "concentrations_mM": present,
        "metabolite_linewidth_ppm": METABOLITE_LW_PPM,
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


def lcmodel_synthetic(basis_path, omitted, hzpppm, deltat, nunfil, linewidth_ppm):
    """Metabolite FID on the data grid (LCModel .RAW convention), its concentrations and a_p.

    Each basis FID gets its ISHIFT, is stretched to the data's HZPPPM (ppm preserved) and
    interpolated onto the data dwell.
    """
    (hz_basis, dt_basis, n_basis), basis = read_lcmodel_basis(basis_path)
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

    conc, ishift, spectrum = basis["Cr"]
    a_p = line_amplitude(
        np.fft.ifft(spectrum), dt_basis, (4.65 - CR_CH3_PPM) * hz_basis + ishift * df
    )
    present = truth(set(basis) - set(omitted))
    met = sum(
        c / basis[name][0] * on_data_grid(basis[name][2], basis[name][1])
        for name, c in present.items()
    )
    t = np.arange(nunfil) * deltat
    return met * np.exp(-np.pi * linewidth_ppm * hzpppm * t), present, a_p / CR_CH3_PROTONS / conc


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
    """The basis FIDs as the fitter fits them, their names, the dwell time, a_p and the file."""
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
    fids = np.asarray(basis.fids)
    names = list(basis.metabolite_names)
    a_p = (
        line_amplitude(fids[:, names.index("Cr")], dwell, (CR_CH3_PPM - 4.7) * larmor_mhz)
        / CR_CH3_PROTONS
    )
    return fids, names, float(dwell), a_p, os.path.normpath(path)


def fit_deep_synthetic(
    fitter, info, water_fit, water_amplitude, n_points, seed=0, linewidth_ppm=METABOLITE_LW_PPM
):
    """The synthetic spectrum through PHIVE (dlfit) or SpatialRegu (gpufit) on a small grid of
    identical voxels: F per group, A_w, the concentrations and the inputs."""
    import deep_crt_mrsi.deepmrsi as fire

    info = dict(info, b0_corrected=True)
    if fitter == "dlfit":
        fire._check_field_matches_basis(info, True)
    fit = (fire._dlfit_fitter if fitter == "dlfit" else fire._gpufit_fitter)(info, False)
    fids, names, dt, a_p, basis_path = deep_basis(fitter, info)
    present = truth(names)
    met = np.zeros(n_points, complex)
    for name, c in present.items():
        column = fids[:n_points, names.index(name)]
        met[: column.size] += c * column
    larmor = float(info["larmor_frequency"])
    met *= np.exp(-np.pi * linewidth_ppm * larmor * 1e-6 * np.arange(n_points) * dt)
    met *= water_amplitude / (2 * WATER_MM * a_p)
    # The fitters conjugate their input, so the data carry the conjugate of the model.
    data = add_noise(np.broadcast_to(np.conj(met), DEEP_GRID + (n_points,)), seed).astype(
        np.complex64
    )
    mask = np.ones(DEEP_GRID, dtype=bool)
    images, out_names = fit(data, mask)[:2]
    out_names = list(out_names)
    fitted = {
        group: float(np.median(np.asarray(images[out_names.index(name)])[mask]))
        for group, name in DEEP_NAMES.items()
        if name in out_names
    }
    estimated = water_fit(synthetic_water(water_amplitude, n_points, dt, larmor), dt)
    return (
        fitted,
        estimated,
        present,
        {
            "a_p": a_p,
            "basis": basis_path,
            "basis_sha256": sha256(basis_path),
            "larmor_hz": larmor,
            "dwell_s": dt,
            "points": n_points,
            "fitter": fitter,
        },
    )


def calibrate_deep(fitter, info, water_fit, water_amplitude, n_points):
    """C_fit of PHIVE (dlfit) or SpatialRegu (gpufit) for the data's field and grid."""
    fitted, estimated, present, inputs = fit_deep_synthetic(
        fitter, info, water_fit, water_amplitude, n_points
    )
    c_fit, per_group = constants(fitted, estimated, present)
    return c_fit, dict(
        inputs, C_per_metabolite=per_group, synthetic=synthetic_record(present, water_amplitude)
    )
