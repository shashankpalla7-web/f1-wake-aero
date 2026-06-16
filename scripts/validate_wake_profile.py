# sanity checks for the wake model (derivation pt 5)
#   1. profiles at different x/d collapse onto one gaussian (self-similarity)
#   2. Delta_u/U_inf ~ (x/d)^-1/2, delta/d ~ (x/d)^+1/2 on log-log
#   3. integral of u'(x,y) dy is constant across x (momentum conservation)

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wake import WakeParameters, half_width, centerline_deficit, velocity_profile, momentum_deficit

FIGURES = Path(__file__).resolve().parents[1] / "figures"
FIGURES.mkdir(exist_ok=True)

# anchored C_D values, DRS closed vs open
# provenance: vault research/cd-drag-coefficient-sourcing.md (whole-car, frontal area A ~ 1.4 m^2)
U_INF = 60.0     # m/s, ~215 km/h
D_REF = 1.0      # reference length, normalized so x/d is what matters
DELTA_D = D_REF  # half-width = d at x = d, starting guess (pt 3)

closed = WakeParameters(U_inf=U_INF, C_D=0.90, d=D_REF, delta_d=DELTA_D)  # baseline, band 0.70-1.10
open_ = WakeParameters(U_inf=U_INF, C_D=0.72, d=D_REF, delta_d=DELTA_D)   # closed * (1 - 0.20 DRS)

x_over_d = np.array([1, 2, 5, 10, 20], dtype=float)
x_vals = x_over_d * D_REF


# --- check 1: self similarity ---
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
ax.set_title("self-similarity check")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / "01_self_similarity.png", dpi=150)
plt.close(fig)


# --- check 2: scaling exponents ---
du_over_uinf = np.array([centerline_deficit(x, closed) / U_INF for x in x_vals])
delta_over_d = np.array([half_width(x, closed) / D_REF for x in x_vals])

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.loglog(x_over_d, du_over_uinf, "o-", label="Delta_u / U_inf")
ax.loglog(x_over_d, delta_over_d, "s-", label="delta / d")
ax.loglog(x_over_d, du_over_uinf[0] * x_over_d ** -0.5, "k--", lw=1, label="slope -1/2")
ax.loglog(x_over_d, delta_over_d[0] * x_over_d ** 0.5, "k:", lw=1, label="slope +1/2")

ax.set_xlabel("x / d")
ax.set_ylabel("normalized quantity")
ax.set_title("scaling check (log-log)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / "02_scaling.png", dpi=150)
plt.close(fig)


# --- check 3: momentum conservation ---
print("momentum check (should be ~constant across x):")
for xod, x in zip(x_over_d, x_vals):
    print(f"  x/d = {xod:5.1f}  ->  {momentum_deficit(x, closed):.6f}")

print()
print("scaling check:")
for xod, du, delta in zip(x_over_d, du_over_uinf, delta_over_d):
    print(f"  x/d = {xod:5.1f}  ->  Delta_u/U_inf = {du:.4f}   delta/d = {delta:.4f}")

print()
print(f"figures written to {FIGURES}")
