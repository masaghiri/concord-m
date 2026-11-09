"""Physics-simulated hydroxyapatite control pattern: one lattice, Cu Kα1 + Kα2, Scherrer broadening, Poisson counts.

Correct data must PASS every rule (or return N/A); this is how the rules are checked for false flags.
"""
import numpy as np

from . import xrd_checks as xc

# Approximate relative intensities of the main hydroxyapatite reflections (calculated reference pattern)
REL_I = {(0, 0, 2): 49.0, (1, 0, 2): 12.2, (2, 1, 0): 39.8, (2, 1, 1): 100.0, (1, 1, 2): 92.9, (3, 0, 0): 90.4,
         (2, 0, 2): 22.5, (3, 1, 0): 42.2, (1, 1, 3): 13.6, (2, 0, 3): 10.5, (2, 2, 2): 45.0, (3, 1, 2): 18.4,
         (2, 1, 3): 60.1, (3, 2, 1): 28.3, (4, 1, 0): 17.1, (4, 0, 2): 17.6, (0, 0, 4): 18.8, (1, 0, 0): 20.8,
         (1, 1, 0): 15.8, (2, 0, 0): 22.3, (1, 1, 1): 8.0, (3, 0, 2): 2.6, (3, 1, 3): 8.7, (5, 0, 2): 13.1,
         (3, 0, 4): 10.1, (3, 2, 3): 14.1, (2, 4, 0): 10.2}


def control_pattern(a=9.416, c=6.893, size_nm=50.0, seed=20261109):
    """Return (two_theta, counts) for a physically consistent hydroxyapatite powder pattern."""
    rng = np.random.default_rng(seed)
    tt = np.round(np.arange(10, 70.0001, 0.02), 2)
    y = 400 * np.exp(-(tt - 10) / 12) + 120  # background (counts)
    for hkl, rel in REL_I.items():
        d = xc.d_hex(*hkl, a, c)
        for lam, frac in ((xc.KA1, 1.0), (xc.KA2, 0.5)):
            t0 = xc.two_theta(d, lam)
            th = np.radians(t0 / 2)
            fwhm = np.degrees(0.9 * 0.15406 / (size_nm * np.cos(th))) + 0.06
            y = y + 25 * rel * frac * xc._pv(tt, t0, fwhm)
    return tt, rng.poisson(y).astype(float)
