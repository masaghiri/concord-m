"""Public pilot figure: (A) the blind AI forgery with rule results, (B) the planned cross-technique test.

Usage (from the repository root, after run_pilot.py):  python scripts/make_figure.py
Writes: figures/pilot_ai_forgery.png
"""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
RED, GREEN, NAVY, GREY = "#e34948", "#1baf7a", "#1F3A5F", "#7a7a7a"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})

res = json.loads((ROOT / "results" / "pilot_results.json").read_text(encoding="utf-8"))
rows = list(csv.DictReader(open(ROOT / "data" / "ai_forgery_hydroxyapatite_xrd.csv", encoding="utf-8")))
tt = np.array([float(r["two_theta"]) for r in rows])
y = np.array([float(r["intensity"]) for r in rows])

fig = plt.figure(figsize=(8.0, 3.0), dpi=200)
ax = fig.add_axes([0.07, 0.15, 0.62, 0.72])
axc = fig.add_axes([0.73, 0.04, 0.26, 0.84])
ax.plot(tt, y, color=NAVY, lw=0.6)
ax.set_xlim(tt.min(), tt.max())
ax.set_xlabel("2θ (degrees)")
ax.set_ylabel("Intensity (counts)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_title("A  AI-made “hydroxyapatite” XRD (blind agent)", loc="left", fontweight="bold", color=NAVY)
names = {"R1": "One crystal lattice (Bragg)", "R2": "Cu Kα1/Kα2 doublet", "R3": "Counting (Poisson) noise",
         "R4": "Whole-number counts"}
yt = 0.97
for r in res["ai_forgery"]:
    col = {"PASS": GREEN, "FLAG": RED}.get(r["status"], GREY)
    mark = {"PASS": "✓", "FLAG": "✗"}.get(r["status"], "–")
    ax.text(0.99, yt, f"{names[r['rule'][:2]]} {mark}", transform=ax.transAxes, ha="right", va="top", color=col,
            fontsize=7.5)
    yt -= 0.09
ax.text(0.99, yt - 0.01, "Passes every XRD-only check:\nthe agent simulated the crystal", transform=ax.transAxes,
        ha="right", va="top", color=RED, fontsize=7)

axc.set_xlim(0, 1)
axc.set_ylim(0, 1)
axc.axis("off")
axc.set_title("B  Next test", loc="left", fontweight="bold", color=NAVY)
for yb, col, t1, t2 in ((0.68, RED, "AI forger", "XRD + EDS + SEM bundle"), (0.25, GREEN, "CONCORD-M", "coupled checks")):
    axc.add_patch(FancyBboxPatch((0.06, yb), 0.88, 0.22, boxstyle="round,pad=0.01,rounding_size=0.04", fc=col,
                                 ec="none"))
    axc.text(0.5, yb + 0.145, t1, ha="center", va="center", color="white", fontweight="bold")
    axc.text(0.5, yb + 0.06, t2, ha="center", va="center", color="white", fontsize=7)
axc.add_patch(FancyArrowPatch((0.33, 0.67), (0.33, 0.48), arrowstyle="-|>", mutation_scale=8, color="#444"))
axc.add_patch(FancyArrowPatch((0.67, 0.48), (0.67, 0.67), arrowstyle="-|>", mutation_scale=8, color="#444"))
axc.text(0.30, 0.575, "fake", ha="right", va="center", fontsize=7, color="#444")
axc.text(0.70, 0.575, "flag", ha="left", va="center", fontsize=7, color="#444")
axc.text(0.5, 0.02, "Does requiring agreement\nraise the cost of forgery?", ha="center", va="bottom", fontsize=7,
         color=GREY)

(ROOT / "figures").mkdir(exist_ok=True)
out = ROOT / "figures" / "pilot_ai_forgery.png"
fig.savefig(out, facecolor="white")
print("wrote", out)
