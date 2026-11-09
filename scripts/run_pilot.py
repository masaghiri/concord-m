"""Reproduce the CONCORD-M pilot: rules on a physics control (should pass) and on a blind AI forgery.

Usage (from the repository root):  python scripts/run_pilot.py
Writes: data/physics_control_hydroxyapatite_xrd.csv and results/pilot_results.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from concord_m import xrd_checks as xc  # noqa: E402
from concord_m.control import control_pattern  # noqa: E402


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return np.array([float(r["two_theta"]) for r in rows]), np.array([float(r["intensity"]) for r in rows])


def summary(res):
    return [{"rule": r["rule"], "status": r["status"], "evidence": r["evidence"]} for r in res]


if __name__ == "__main__":
    out = {}
    tt, y = control_pattern()
    np.savetxt(ROOT / "data" / "physics_control_hydroxyapatite_xrd.csv", np.c_[tt, y], delimiter=",",
               header="two_theta,intensity", comments="", fmt=["%.2f", "%d"])
    out["physics_control"] = summary(xc.audit(tt, y))
    tt, y = load(ROOT / "data" / "ai_forgery_hydroxyapatite_xrd.csv")
    out["ai_forgery"] = summary(xc.audit(tt, y))
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "pilot_results.json").write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                         encoding="utf-8")
    for name, res in out.items():
        print(f"\n== {name}")
        for r in res:
            print(f"  {r['status']:4s}  {r['rule']}: {r['evidence']}")
