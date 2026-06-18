# C_D sensitivity sweep -- wake centerline deficit vs following distance
# anchored drag coefficients, provenance: vault research/cd-drag-coefficient-sourcing.md
#   closed DRS  C_D = 0.90   (whole-car, frontal area A ~ 1.4 m^2; band 0.70-1.10)
#   open DRS    C_D = 0.72   (closed * (1 - 0.20), ~20% whole-car drag reduction)
# since no single public whole-car C_D exists, the result is reported as a band
# from a sweep, not a single curve (see decision 2026-06-16 in vault).

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wake import WakeParameters, centerline_deficit

FIGURES = Path(__file__).resolve().parents[1] / "figures"
FIGURES.mkdir(exist_ok=True)

U_INF = 60.0
D_REF = 1.0
DELTA_D = D_REF

# anchored baselines - from research in vault, see research/cd-drag-coefficient-sourcing.md
CD_CLOSED = 0.90
DRS_REDUCTION = 0.20
CD_OPEN = CD_CLOSED * (1 - DRS_REDUCTION)   # 0.72

# sweep ranges
CD_CLOSED_RANGE = [0.70, 0.90, 1.10]
DRS_RANGE = [0.15, 0.20, 0.25]

x_over_d = np.array([1, 2, 3, 5, 7, 10, 15, 20], dtype=float)
x_vals = x_over_d * D_REF


def deficit_frac(cd, x):
    # centerline speed loss as a fraction of freestream
    p = WakeParameters(U_inf=U_INF, C_D=cd, d=D_REF, delta_d=DELTA_D)
    return centerline_deficit(x, p) / U_INF


closed_frac = np.array([deficit_frac(CD_CLOSED, x) for x in x_vals])
open_frac = np.array([deficit_frac(CD_OPEN, x) for x in x_vals])

# band: every closed C_D in the range, plus every closed*DRS open value
all_cd = set(CD_CLOSED_RANGE)
for cd in CD_CLOSED_RANGE:
    for r in DRS_RANGE:
        all_cd.add(round(cd * (1 - r), 4))
all_cd = sorted(all_cd)
band = np.array([[deficit_frac(cd, x) for cd in all_cd] for x in x_vals])
band_lo = band.min(axis=1)
band_hi = band.max(axis=1)

print("wake centerline velocity deficit -- how much slower the air is at the wake center")
print("(fraction of clean-air speed, so 0.21 = 21% slower)")
print()
print(f"{'x/d':>5} {'closed 0.90':>12} {'open 0.72':>11} {'band lo':>9} {'band hi':>9}")
for i, xod in enumerate(x_over_d):
    print(f"{xod:5.0f} {closed_frac[i]*100:11.1f}% {open_frac[i]*100:10.1f}% "
          f"{band_lo[i]*100:8.1f}% {band_hi[i]*100:8.1f}%")

print()
cut = (1 - open_frac[0] / closed_frac[0]) * 100
print(f"DRS effect at x/d = 1:  closed {closed_frac[0]*100:.1f}%  ->  open {open_frac[0]*100:.1f}%  "
      f"(wake deficit cut by {cut:.0f}%)")
print(f"sweep spans C_D = {all_cd[0]:.3f} to {all_cd[-1]:.3f}")

# plot
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.fill_between(x_over_d, band_lo * 100, band_hi * 100,
                color="0.8", alpha=0.7, label="C_D / DRS sensitivity band")
ax.plot(x_over_d, closed_frac * 100, "o-", color="C3", label="DRS closed (C_D = 0.90)")
ax.plot(x_over_d, open_frac * 100, "s-", color="C0", label="DRS open (C_D = 0.72)")
ax.set_xlabel("following distance   x / d")
ax.set_ylabel("centerline velocity deficit   (%)")
ax.set_title("wake deficit vs following distance (anchored C_D)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "03_cd_sensitivity.png", dpi=150)
plt.close(fig)
print(f"\nfigure written to {FIGURES / '03_cd_sensitivity.png'}")
