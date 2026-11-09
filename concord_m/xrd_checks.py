"""CONCORD-M pilot (formerly checkMATE): physics checks for a powder XRD pattern claimed to be hydroxyapatite (Cu Kα).

Each rule returns a dict: {rule, status: PASS | FLAG | N/A, evidence}. Rules use physics only:
  R1 single lattice   - all peaks must be explained by ONE hexagonal P6_3/m lattice (Bragg's law)
  R2 Cu Kα doublet    - lab Cu Kα data carry a Kα2 component (λ = 1.544414 Å, about half of Kα1)
  R3 counting noise   - detector counts follow Poisson statistics: noise grows as sqrt(intensity)
  R4 integer counts   - raw counts are whole numbers
Validated first on a physics-simulated control (see run_pilot.py) so that correct data pass.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks, savgol_filter

KA1, KA2 = 1.540593, 1.544414  # Cu Kα1, Kα2 wavelengths (Å)
HA_A, HA_C = 9.418, 6.884       # stoichiometric hydroxyapatite reference lattice (Å)


def d_hex(h, k, l, a, c):
    return 1.0 / np.sqrt(4.0 / 3.0 * (h * h + h * k + k * k) / a ** 2 + l * l / c ** 2)


def two_theta(d, lam=KA1):
    s = lam / (2.0 * d)
    return np.degrees(2.0 * np.arcsin(s)) if s < 1 else np.nan


def allowed_reflections(a=HA_A, c=HA_C, tt_max=72.0):
    """Unique d-spacings allowed in P6_3/m (only condition: 00l with l = 2n)."""
    seen = {}
    for h in range(0, 9):
        for k in range(0, h + 1):
            for l in range(0, 9):
                if h == k == 0 and (l == 0 or l % 2):
                    continue
                key = (h * h + h * k + k * k, l)
                if key in seen:
                    continue
                tt = two_theta(d_hex(h, k, l, a, c))
                if np.isfinite(tt) and tt <= tt_max:
                    seen[key] = (h, k, l)
    return list(seen.values())


def find_xrd_peaks(tt, y, min_prom_frac=0.04):
    ys = savgol_filter(y, 11, 3)
    prom = min_prom_frac * (ys.max() - np.percentile(ys, 5))
    idx, props = find_peaks(ys, prominence=prom, distance=5)
    peaks = []
    for i in idx:
        lo, hi = max(i - 3, 0), min(i + 4, len(tt))
        p = np.polyfit(tt[lo:hi] - tt[i], ys[lo:hi], 2)
        x0 = tt[i] - p[1] / (2 * p[0]) if p[0] < 0 else tt[i]
        peaks.append(float(x0))
    return np.array(peaks), ys


def peak_fwhm(tt, ys, x0):
    """FWHM of a peak in the smoothed curve, measured against the local minimum within ±1°."""
    m = (tt > x0 - 1.0) & (tt < x0 + 1.0)
    base = ys[m].min()
    i0 = int(np.argmin(np.abs(tt - x0)))
    half = base + (ys[i0] - base) / 2
    lo = i0
    while lo > 0 and ys[lo] > half:
        lo -= 1
    hi = i0
    while hi < len(ys) - 1 and ys[hi] > half:
        hi += 1
    return float(tt[hi] - tt[lo])


def rule_single_lattice(tt, y):
    """R1: index every resolved peak with one (a, c); flag peaks one lattice cannot explain.

    Scope: only peaks whose nearest other allowed reflection is more than 1.5 FWHM away (overlapping broad
    peaks shift apparent maxima). Nanocrystalline patterns with too few resolved peaks return N/A; they need
    whole-pattern fitting instead.
    """
    peaks, ys = find_xrd_peaks(tt, y)
    peaks = peaks[(peaks > 20) & (peaks < 70)]
    widths = np.array([peak_fwhm(tt, ys, p) for p in peaks])
    refl = allowed_reflections()
    ref0 = np.sort(np.array([two_theta(d_hex(*hkl, HA_A, HA_C)) for hkl in refl]))
    keep = []
    for p, w in zip(peaks, widths):
        near = np.sort(np.abs(ref0 - p))
        keep.append(len(near) > 1 and near[1] > 1.5 * w)  # near[0] is the peak's own reflection
    keep = np.array(keep, dtype=bool)
    all_peaks = peaks
    peaks, widths = peaks[keep], widths[keep]
    if len(peaks) < 4:
        return {"rule": "R1 single lattice (Bragg)", "status": "N/A",
                "evidence": f"only {len(peaks)} resolved peaks (broad, overlapping pattern)", "bad": [],
                "peaks": all_peaks}

    def assign(a, c):
        ref = np.array([two_theta(d_hex(*hkl, a, c)) for hkl in refl])
        j = np.abs(peaks[:, None] - ref[None, :]).argmin(axis=1)
        return j, ref[j]

    j, _ = assign(HA_A, HA_C)

    def resid(p):
        return peaks - np.array([two_theta(d_hex(*refl[m], p[0], p[1])) for m in j])

    fit = least_squares(resid, [HA_A, HA_C], bounds=([9.2, 6.7], [9.6, 7.1]))
    a, c = fit.x
    j, ref = assign(a, c)
    dev = peaks - ref
    tols = np.maximum(0.06, 0.25 * widths)
    bad = [(round(float(p), 2), round(float(dv), 2), refl[m]) for p, dv, m, t in zip(peaks, dev, j, tols)
           if abs(dv) > t]
    status = "FLAG" if len(bad) >= 2 else "PASS"
    ev = (f"{len(peaks)} resolved peaks; best single lattice a = {a:.3f} Å, c = {c:.3f} Å; "
          f"rms misfit {np.sqrt(np.mean(dev ** 2)):.3f}°; {len(bad)} peaks off beyond tolerance "
          f"(worst {max((abs(b[1]) for b in bad), default=0):.2f}°)")
    return {"rule": "R1 single lattice (Bragg)", "status": status, "evidence": ev, "bad": bad,
            "peaks": all_peaks, "a": a, "c": c}


def _pv(x, x0, w, eta=0.5):
    g = np.exp(-4 * np.log(2) * ((x - x0) / w) ** 2)
    lz = 1.0 / (1.0 + 4 * ((x - x0) / w) ** 2)
    return eta * lz + (1 - eta) * g


def rule_kalpha2(tt, y, peaks):
    """R2: for isolated peaks above 40°, compare a single peak with a fixed Kα1/Kα2 doublet (same parameter count)."""
    wins, tested = 0, 0
    details = []
    for x0 in peaks[peaks > 40]:
        if np.min(np.abs(peaks[peaks != x0] - x0), initial=9) < 0.8:
            continue
        m = (tt > x0 - 0.6) & (tt < x0 + 0.6)
        x, yy = tt[m], y[m]
        if len(x) < 20:
            continue
        th = np.radians(x0 / 2)
        sep = np.degrees(2 * np.tan(th) * (KA2 - KA1) / KA1)

        def single(p):
            return p[0] * _pv(x, p[1], p[2]) + p[3] + p[4] * (x - x0) - yy

        def doublet(p):
            return p[0] * (_pv(x, p[1], p[2]) + 0.5 * _pv(x, p[1] + sep, p[2])) + p[3] + p[4] * (x - x0) - yy

        base = np.percentile(yy, 10)
        p0 = [yy.max() - base, x0, 0.15, base, 0.0]
        lb, ub = [0, x0 - 0.3, 0.02, -np.inf, -np.inf], [np.inf, x0 + 0.3, 1.0, np.inf, np.inf]
        rs = least_squares(single, p0, bounds=(lb, ub))
        rd = least_squares(doublet, p0, bounds=(lb, ub))
        fwhm = rs.x[2]
        if sep < 0.25 * fwhm:  # doublet not resolvable for very broad peaks
            continue
        tested += 1
        better = np.sum(rd.fun ** 2) < np.sum(rs.fun ** 2)
        wins += int(better)
        details.append((round(float(x0), 2), round(float(sep), 3), round(float(fwhm), 3), bool(better)))
    if tested == 0:
        return {"rule": "R2 Cu Kα1/Kα2 doublet", "status": "N/A", "evidence": "no isolated, resolvable peak above 40°",
                "details": details}
    status = "FLAG" if wins <= tested / 3 else "PASS"
    ev = f"doublet model fits better in {wins} of {tested} isolated peaks above 40° (expected: most)"
    return {"rule": "R2 Cu Kα1/Kα2 doublet", "status": status, "evidence": ev, "details": details}


def rule_poisson(tt, y):
    """R3: noise must scale as sqrt(intensity). Slope of log(noise) vs log(level): Poisson ≈ 0.5."""
    sm = savgol_filter(y, 9, 3)
    r = y - sm
    lvl = np.clip(sm, 1e-9, None)
    edges = np.quantile(lvl, np.linspace(0, 1, 11))
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (lvl >= lo) & (lvl < hi)
        if m.sum() > 30:
            xs.append(np.log(np.median(lvl[m])))
            ys.append(np.log(np.std(r[m]) + 1e-12))
    slope = np.polyfit(xs, ys, 1)[0] if len(xs) >= 4 else np.nan
    span = np.exp(max(xs) - min(xs)) if xs else 1
    if not np.isfinite(slope) or span < 3:
        return {"rule": "R3 counting (Poisson) noise", "status": "N/A", "evidence": "intensity range too narrow"}
    status = "PASS" if 0.3 <= slope <= 0.7 else "FLAG"
    ev = f"noise grows with intensity^{slope:.2f} (Poisson counting: ^0.50)"
    return {"rule": "R3 counting (Poisson) noise", "status": status, "evidence": ev, "slope": slope}


def rule_integer_counts(y):
    """R4: raw detector counts are whole numbers."""
    frac = float(np.mean(np.abs(y - np.round(y)) > 1e-6))
    status = "FLAG" if frac > 0.05 else "PASS"
    return {"rule": "R4 whole-number counts", "status": status,
            "evidence": f"{frac * 100:.0f}% of points are not whole numbers"}


def audit(tt, y):
    r1 = rule_single_lattice(tt, y)
    return [r1, rule_kalpha2(tt, y, r1["peaks"]), rule_poisson(tt, y), rule_integer_counts(y)]
