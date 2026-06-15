"""Validation checks for the wake model -- see derivations/01_wake_profile_derivation.md §5.

1. Profile self-similarity: u'(x,y)/Delta_u(x) vs y/delta(x) collapses to one Gaussian.
2. Scaling: Delta_u/U_inf ~ (x/d)^-1/2, delta/d ~ (x/d)^+1/2 on log-log axes.
3. Momentum conservation: integral of u'(x,y) dy is constant across x.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wake import WakeParameters, half_width, centerline_deficit, velocity_profile, momentum_deficit

FIGURES = Path(__file__).resolve().parents[1] / "figures"
FIGURES.mkdir(exist_ok=True)

# Placeholder C_D values -- DRS closed (high downforce) vs open (low drag).
# TODO: replace with literature-sourced values for the leading car's whole-car
# drag coefficient in each DRS state (see research/ in the vault project).
U_INF = 60.0  # m/s, ~215 km/h
D_REF = 1.0   # reference length [m] -- normalized; x/d is the quantity of interest
DELTA_D = D_REF  # spreading param: half-width = d at x = d (derivation §3 starting value)

closed = WakeParameters(U_inf=U_INF, C_D=0.90, d=D_REF, delta_d=DELTA_D)
open_ = WakeParameters(U_inf=U_INF, C_D=0.80, d=D_REF, delta_d=DELTA_D)

x_over_d = np.array([1, 2, 5, 10, 20], dtype=float)
x_vals = x_over_d * D_REF


# ---------------------------------------------------------------------------
# Check 1: self-similarity -- profiles collapse onto exp(-ln2 * eta^2)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4.5))
eta = np.linspace(-3, 3, 400)
ax.plot(eta, np.exp(-np.log(2) * eta**2), "k--", lw=2, label="exp(-ln2 * eta^2)")

for xod, x in zip(x_over_d, x_vals):
    delta = half_width(x, closed)
    du = centerline_deficit(x, closed)
    y = eta * delta
    deficit = closed.U_inf - velocity_profile(x, y, closed)
    ax.plot(eta, deficit / du, alpha=0.7, label=f"x/d = {xod:g}")

ax.set_xlabel("eta = y / delta(x)")
ax.set_ylabel("u'(x,y) / Delta_u(x)")
ax.set_title("Self-similarity check: all profiles collapse to one Gaussian")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / "01_self_similarity.png", dpi=150)
plt.close(fig)


# ---------------------------------------------------------------------------
# Check 2: scaling exponents on log-log axes
# ---------------------------------------------------------------------------
du_over_uinf = np.array([centerline_deficit(x, closed) / U_INF for x in x_vals])
delta_over_d = np.array([half_width(x, closed) / D_REF for x in x_vals])

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.loglog(x_over_d, du_over_uinf, "o-", label="Delta_u / U_inf (data)")
ax.loglog(x_over_d, delta_over_d, "s-", label="delta / d (data)")

# reference slopes anchored at x/d = 1
ax.loglog(x_over_d, du_over_uinf[0] * x_over_d ** -0.5, "k--", lw=1, label="slope -1/2")
ax.loglog(x_over_d, delta_over_d[0] * x_over_d ** 0.5, "k:", lw=1, label="slope +1/2")

ax.set_xlabel("x / d")
ax.set_ylabel("normalized quantity")
ax.set_title("Scaling check (log-log)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / "02_scaling.png", dpi=150)
plt.close(fig)


# ---------------------------------------------------------------------------
# Check 3: momentum conservation -- integral of u'(x,y) dy constant across x
# ---------------------------------------------------------------------------
print("Momentum check (integral of u'(x,y) dy, should be ~constant across x):")
for xod, x in zip(x_over_d, x_vals):
    print(f"  x/d = {xod:5.1f}  ->  {momentum_deficit(x, closed):.6f}")

print()
print("Scaling check (Delta_u/U_inf and delta/d):")
for xod, du, delta in zip(x_over_d, du_over_uinf, delta_over_d):
    print(f"  x/d = {xod:5.1f}  ->  Delta_u/U_inf = {du:.4f}   delta/d = {delta:.4f}")

print()
print(f"Figures written to {FIGURES}")
