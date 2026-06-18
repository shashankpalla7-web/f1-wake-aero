# Stage 2 -- effective inflow U_eff(x) and downforce loss vs following distance
# derivation: derivations/02_effective_inflow_derivation.md
#
# reduces the wake field u(x,y) to the single inflow speed the following car's wing
# sees, then to the downforce it loses to dirty air. headline reduction is the RMS
# (dynamic-pressure-equivalent) inflow; centerline and span-mean are reported as the
# sensitivity range, same posture as the C_D sourcing.
#
# anchored C_D, provenance: vault research/cd-drag-coefficient-sourcing.md
#   closed DRS C_D = 0.90 (band 0.70-1.10);  open DRS C_D = 0.72 (closed * (1 - 0.20))

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wake import (
    WakeParameters,
    centerline_inflow, span_mean_inflow, rms_inflow,
    effective_inflow, downforce_loss_fraction,
)

FIGURES = Path(__file__).resolve().parents[1] / "figures"
FIGURES.mkdir(exist_ok=True)

U_INF = 60.0      # m/s, ~215 km/h
D_REF = 1.0
DELTA_D = D_REF
SPAN = 2.0        # following front-wing span, m (FIA reg max 2000 mm; derivation 2 sec 1)

CD_CLOSED = 0.90
DRS_REDUCTION = 0.20
CD_OPEN = CD_CLOSED * (1 - DRS_REDUCTION)   # 0.72

# sweep for the sensitivity band (mirrors scripts/cd_sensitivity.py)
CD_CLOSED_RANGE = [0.70, 0.90, 1.10]
DRS_RANGE = [0.15, 0.20, 0.25]

x_over_d = np.array([1, 2, 3, 5, 7, 10, 15, 20], dtype=float)
x_vals = x_over_d * D_REF

closed = WakeParameters(U_inf=U_INF, C_D=CD_CLOSED, d=D_REF, delta_d=DELTA_D)
open_ = WakeParameters(U_inf=U_INF, C_D=CD_OPEN, d=D_REF, delta_d=DELTA_D)


def loss(params, x):
    return downforce_loss_fraction(x, params, SPAN, method="rms")


# ----------------------------------------------------------------------------
# main table: U_eff/U_inf (three reductions) + downforce loss, DRS closed vs open
# ----------------------------------------------------------------------------
print("Stage 2 -- effective inflow and downforce loss vs following distance")
print(f"U_inf = {U_INF:.0f} m/s,  wing span b = {SPAN:.1f} m,  reductions of the wake profile\n")
print("U_eff / U_inf  (1.00 = clean air):")
print(f"{'x/d':>5} | {'cl_center':>9} {'cl_mean':>8} {'cl_RMS':>7} | {'op_center':>9} {'op_mean':>8} {'op_RMS':>7}")
for xod, x in zip(x_over_d, x_vals):
    cc = centerline_inflow(x, closed) / U_INF
    cm = span_mean_inflow(x, closed, SPAN) / U_INF
    cr = rms_inflow(x, closed, SPAN) / U_INF
    oc = centerline_inflow(x, open_) / U_INF
    om = span_mean_inflow(x, open_, SPAN) / U_INF
    orr = rms_inflow(x, open_, SPAN) / U_INF
    print(f"{xod:5.0f} | {cc:9.3f} {cm:8.3f} {cr:7.3f} | {oc:9.3f} {om:8.3f} {orr:7.3f}")

print("\ndownforce loss  1 - (U_eff/U_inf)^2  using RMS inflow (the headline):")
print(f"{'x/d':>5} | {'closed':>8} {'open':>8} | {'band lo':>8} {'band hi':>8}")

# sensitivity band over the full C_D x DRS sweep
all_cd = set(CD_CLOSED_RANGE)
for cd in CD_CLOSED_RANGE:
    for r in DRS_RANGE:
        all_cd.add(round(cd * (1 - r), 4))
all_cd = sorted(all_cd)

closed_loss = np.array([loss(closed, x) for x in x_vals])
open_loss = np.array([loss(open_, x) for x in x_vals])
band = np.array([[loss(WakeParameters(U_INF, cd, D_REF, DELTA_D), x) for cd in all_cd] for x in x_vals])
band_lo, band_hi = band.min(axis=1), band.max(axis=1)

for i, xod in enumerate(x_over_d):
    print(f"{xod:5.0f} | {closed_loss[i]*100:7.1f}% {open_loss[i]*100:7.1f}% | "
          f"{band_lo[i]*100:7.1f}% {band_hi[i]*100:7.1f}%")

print()
print(f"HEADLINE  at x/d = 1 (one car length back): following car loses "
      f"{closed_loss[0]*100:.0f}% of downforce (DRS closed) / {open_loss[0]*100:.0f}% (DRS open).")
print(f"          by x/d = 20 the penalty falls to {closed_loss[-1]*100:.0f}% / {open_loss[-1]*100:.0f}%.")
drs_gain = (closed_loss[0] - open_loss[0]) * 100
print(f"          DRS recovers ~{drs_gain:.0f} downforce-points at one car length back.")

# ----------------------------------------------------------------------------
# validation checks (derivation 2 sec 5)
# ----------------------------------------------------------------------------
print("\nvalidation checks:")

# 1. limits: b -> 0 collapses to centerline; x -> inf recovers U_inf
tiny = 1e-4
near_center = span_mean_inflow(5.0, closed, tiny)
print(f"  [1a] b->0  span-mean({tiny}) = {near_center:.4f}  vs centerline {centerline_inflow(5.0, closed):.4f}  "
      f"-> {'OK' if abs(near_center - centerline_inflow(5.0, closed)) < 1e-2 else 'FAIL'}")
# deficit decays only as 1/sqrt(x), so the limit has to be taken far out
far = rms_inflow(1e8, closed, SPAN) / U_INF
print(f"  [1b] x->inf (x/d=1e8)  U_eff/U_inf = {far:.6f}  -> {'OK' if far > 0.9999 else 'FAIL'}")

# 2. ordering U_c <= U_mean <= U_rms <= U_inf at every x
ok_order = all(
    centerline_inflow(x, closed) <= span_mean_inflow(x, closed, SPAN) + 1e-9
    <= 1e18 and span_mean_inflow(x, closed, SPAN) <= rms_inflow(x, closed, SPAN) + 1e-9
    <= U_INF + 1e-9
    for x in x_vals
)
print(f"  [2]  ordering U_c <= U_mean <= U_rms <= U_inf  -> {'OK' if ok_order else 'FAIL'}")

# 3. numeric RMS^2 matches the erf closed form (derivation 2 sec 1c)
def rms_sq_closed_form(x, p, b):
    from wake import half_width, centerline_deficit
    delta = half_width(x, p); du = centerline_deficit(x, p)
    ln2 = np.log(2.0)
    t1 = p.U_inf ** 2
    t2 = (2 * p.U_inf * du * delta / b) * np.sqrt(np.pi / ln2) * erf((b / (2 * delta)) * np.sqrt(ln2))
    t3 = (du ** 2 * delta / b) * np.sqrt(np.pi / (2 * ln2)) * erf((b / (2 * delta)) * np.sqrt(2 * ln2))
    return t1 - t2 + t3
cf = np.sqrt(rms_sq_closed_form(3.0, closed, SPAN))
num = rms_inflow(3.0, closed, SPAN)
print(f"  [3]  RMS closed-form {cf:.5f}  vs numeric {num:.5f}  -> {'OK' if abs(cf - num) < 1e-4 else 'FAIL'}")

# 4. 2x rule: at large x, dL/L0 ~ 2*(1 - U_eff/U_inf)
eps = 1 - rms_inflow(20.0, closed, SPAN) / U_INF
approx = 2 * eps
exact = loss(closed, 20.0)
print(f"  [4]  2x rule at x/d=20: 2*eps = {approx*100:.2f}%  vs exact {exact*100:.2f}%  "
      f"-> {'OK' if abs(approx - exact) < 0.01 else 'FAIL'}")

# ----------------------------------------------------------------------------
# figure: two panels -- U_eff/U_inf, and downforce loss with sensitivity band
# ----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# panel 1: effective inflow, three reductions (closed DRS)
center = np.array([centerline_inflow(x, closed) / U_INF for x in x_vals])
mean_ = np.array([span_mean_inflow(x, closed, SPAN) / U_INF for x in x_vals])
rms = np.array([rms_inflow(x, closed, SPAN) / U_INF for x in x_vals])
ax1.axhline(1.0, color="0.6", lw=1, ls=":", label="clean air")
ax1.plot(x_over_d, center, "v--", color="0.5", label="centerline (worst case)")
ax1.plot(x_over_d, mean_, "^--", color="0.4", label="span-mean")
ax1.plot(x_over_d, rms, "o-", color="C3", label="RMS  =  U_eff (headline)")
ax1.set_xlabel("following distance   x / d")
ax1.set_ylabel("effective inflow   U_eff / U_inf")
ax1.set_title("what the following wing sees (DRS closed)")
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# panel 2: downforce loss, closed vs open + band
ax2.fill_between(x_over_d, band_lo * 100, band_hi * 100,
                 color="0.8", alpha=0.7, label="C_D / DRS sensitivity band")
ax2.plot(x_over_d, closed_loss * 100, "o-", color="C3", label="DRS closed (C_D = 0.90)")
ax2.plot(x_over_d, open_loss * 100, "s-", color="C0", label="DRS open (C_D = 0.72)")
ax2.set_xlabel("following distance   x / d")
ax2.set_ylabel("downforce loss   1 - (U_eff/U_inf)²   (%)")
ax2.set_title("downforce lost to dirty air")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(FIGURES / "04_effective_inflow.png", dpi=150)
plt.close(fig)
print(f"\nfigure written to {FIGURES / '04_effective_inflow.png'}")
