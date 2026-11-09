"""Rules must not flag physically correct data, and must flag simple naive fakes."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from concord_m import xrd_checks as xc  # noqa: E402
from concord_m.control import REL_I, control_pattern  # noqa: E402


@pytest.mark.parametrize("size_nm", [8, 15, 25, 50, 100])
def test_no_false_flags_on_physics_control(size_nm):
    tt, y = control_pattern(size_nm=size_nm)
    statuses = [r["status"] for r in xc.audit(tt, y)]
    assert "FLAG" not in statuses


def _naive_fake(jitter_deg=0.0, seed=1):
    """A naive fake: Gaussian peaks near hydroxyapatite positions, no Kα2, constant noise, decimal values."""
    rng = np.random.default_rng(seed)
    tt = np.round(np.arange(10, 70.0001, 0.02), 2)
    y = 300 * np.exp(-(tt - 10) / 12) + 150
    for hkl, rel in REL_I.items():
        t0 = xc.two_theta(xc.d_hex(*hkl, xc.HA_A, xc.HA_C)) + rng.uniform(-jitter_deg, jitter_deg)
        y = y + 25 * rel * np.exp(-4 * np.log(2) * ((tt - t0) / 0.18) ** 2)
    return tt, y + rng.normal(0, 12, tt.size)


def test_naive_fake_flags_noise_and_counts():
    tt, y = _naive_fake()
    res = {r["rule"][:2]: r["status"] for r in xc.audit(tt, y)}
    assert res["R3"] == "FLAG"  # constant noise is not counting noise
    assert res["R4"] == "FLAG"  # counts are not whole numbers


def test_jittered_peaks_fail_single_lattice():
    tt, y = _naive_fake(jitter_deg=0.25, seed=7)
    res = {r["rule"][:2]: r["status"] for r in xc.audit(tt, np.round(np.clip(y, 0, None)))}
    assert res["R1"] == "FLAG"  # peaks no single lattice can explain
