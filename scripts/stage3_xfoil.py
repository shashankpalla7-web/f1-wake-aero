# stage 3 -- XFOIL closes the loop from wake inflow to real downforce loss
#
# stage 2 gave the downforce loss as a pure dynamic-pressure proxy:
#     dL/L0 = 1 - (U_eff/U_inf)^2          (assumes C_l is unchanged in the wake)
#
# but the wake also drops the Reynolds number the following wing sees
# (Re_wake = Re_clean * U_eff/U_inf), and C_l depends on Re. stage 3 measures that
# dependence with XFOIL and folds it in exactly:
#
#     L_wake/L_clean = (U_eff/U_inf)^2 * (C_l,wake / C_l,clean)
#     dL/L0          = 1 - (U_eff/U_inf)^2 * (C_l,wake / C_l,clean)
#
# C_l,clean is XFOIL at Re_clean; C_l,wake is XFOIL at the reduced Re_wake, same
# section and same angle of attack. the ratio r = C_l,wake/C_l,clean is the second-
# order correction stage 2 dropped.
#
# test matrix: a controlled NACA camber ladder (2412/4412/6412, same 12% thickness,
# only camber varies) at a fixed high-load AoA. answers: does the wake's downforce
# penalty depend on how hard the wing is already working? single-element is a proxy
# for the loaded front-wing element -- XFOIL can't do the real multi-element slot
# flow (named in limitations, derivation 03). sections run upright; downforce is
# |lift| and the C_l ratio is sign-independent.

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wake import WakeParameters, rms_inflow, cl_at

FIGURES = Path(__file__).resolve().parents[1] / "figures"
FIGURES.mkdir(exist_ok=True)

# --- wake / flow conditions (same anchors as stage 2) -----------------------
U_INF = 60.0          # m/s
D_REF = 1.0
DELTA_D = D_REF
SPAN = 2.0            # following front-wing span, m (FIA max)
CD_CLOSED = 0.90
CD_OPEN = 0.72        # DRS open (20% whole-car drag reduction)

# --- airfoil / XFOIL conditions ---------------------------------------------
CHORD = 0.25          # m, representative front-wing element chord
NU_AIR = 1.5e-5       # m^2/s, kinematic viscosity of air ~20 C
RE_CLEAN = U_INF * CHORD / NU_AIR        # = 1.0e6 in clean air
ALPHA = 5.0           # deg, fixed high-load operating point across the ladder
RAMP_STEP = 1.0       # deg, alpha ramp increment for convergence

AIRFOILS = ["NACA 2412", "NACA 4412", "NACA 6412"]   # camber ladder: 2% / 4% / 6%

x_over_d = np.array([1, 2, 3, 5, 7, 10, 15, 20], dtype=float)
x_vals = x_over_d * D_REF


def wake_speed_ratios(C_D):
    # U_eff(x)/U_inf along the distance grid, RMS (dynamic-pressure-equivalent) inflow
    p = WakeParameters(U_inf=U_INF, C_D=C_D, d=D_REF, delta_d=DELTA_D)
    return np.array([rms_inflow(x, p, SPAN) / U_INF for x in x_vals])


def cl_curve(airfoil, speed_ratios):
    # XFOIL C_l at each wake-reduced Reynolds number, fixed AoA
    out = []
    for s in speed_ratios:
        cl = cl_at(airfoil, RE_CLEAN * s, ALPHA, step=RAMP_STEP)
        if cl is None:
            raise RuntimeError(f"{airfoil}: XFOIL failed to converge at Re={RE_CLEAN*s:.0f}")
        out.append(cl)
    return np.array(out)


print("Stage 3 -- XFOIL-corrected downforce loss vs following distance")
print(f"U_inf={U_INF:.0f} m/s  chord={CHORD:.2f} m  Re_clean={RE_CLEAN:.2e}  "
      f"AoA={ALPHA:.1f} deg  span={SPAN:.1f} m  (Ncrit=9, free transition)\n")

s_closed = wake_speed_ratios(CD_CLOSED)
s_open = wake_speed_ratios(CD_OPEN)
proxy_closed = 1.0 - s_closed ** 2      # stage-2 loss (airfoil-independent)
proxy_open = 1.0 - s_open ** 2

results = {}   # airfoil -> dict of arrays
for foil in AIRFOILS:
    print(f"  running {foil} ...", flush=True)
    cl_clean = cl_at(foil, RE_CLEAN, ALPHA, step=RAMP_STEP)
    cl_wake_closed = cl_curve(foil, s_closed)
    cl_wake_open = cl_curve(foil, s_open)
    r_closed = cl_wake_closed / cl_clean
    r_open = cl_wake_open / cl_clean
    results[foil] = dict(
        cl_clean=cl_clean,
        r_closed=r_closed, r_open=r_open,
        loss3_closed=1.0 - s_closed ** 2 * r_closed,
        loss3_open=1.0 - s_open ** 2 * r_open,
    )

# --- table: the C_l correction and the corrected loss, DRS closed -----------
print("\nDRS closed -- C_l ratio (wake/clean) and downforce loss:")
hdr = f"{'x/d':>4} | {'U_eff/U':>7} {'proxy%':>7} |"
for foil in AIRFOILS:
    tag = foil.split()[-1]
    hdr += f" {tag+' r':>9} {tag+' L%':>7} |"
print(hdr)
for i, xod in enumerate(x_over_d):
    row = f"{xod:4.0f} | {s_closed[i]:7.3f} {proxy_closed[i]*100:6.1f}% |"
    for foil in AIRFOILS:
        R = results[foil]
        row += f" {R['r_closed'][i]:9.4f} {R['loss3_closed'][i]*100:6.1f}% |"
    print(row)

print("\nheadline at x/d = 1 (one car length back), DRS closed:")
for foil in AIRFOILS:
    R = results[foil]
    d2 = proxy_closed[0] * 100
    d3 = R['loss3_closed'][0] * 100
    print(f"  {foil}: C_l,clean={R['cl_clean']:.3f}  C_l ratio={R['r_closed'][0]:.4f}  "
          f"-> stage-2 proxy {d2:.1f}%  |  stage-3 corrected {d3:.1f}%  "
          f"(Re effect adds {d3-d2:+.1f} pts)")

mean_r1 = np.mean([results[f]['r_closed'][0] for f in AIRFOILS])
print(f"\n  the C_l ratio at x/d=1 is ~{mean_r1:.3f} across the ladder: the wake's lower Re")
print(f"  costs only ~{(1-mean_r1)*100:.1f}% of C_l, so the dynamic-pressure proxy was the")
print(f"  dominant term -- stage 3 confirms it and nudges the loss slightly deeper.")

# --- figure: two panels -----------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
colors = {"NACA 2412": "C0", "NACA 4412": "C1", "NACA 6412": "C3"}

# panel 1: the C_l correction factor r(x) = C_l,wake / C_l,clean (DRS closed)
ax1.axhline(1.0, color="0.6", lw=1, ls=":", label="no Re effect (stage-2 assumption)")
for foil in AIRFOILS:
    ax1.plot(x_over_d, results[foil]['r_closed'], "o-", color=colors[foil],
             label=f"{foil}  (clean $C_l$={results[foil]['cl_clean']:.2f})")
ax1.set_xlabel("following distance   x / d")
ax1.set_ylabel(r"$C_{l,\mathrm{wake}} / C_{l,\mathrm{clean}}$")
ax1.set_title("XFOIL correction: wake Reynolds-number effect on $C_l$")
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# panel 2: downforce loss -- stage-2 proxy vs stage-3 corrected (DRS closed)
ax2.plot(x_over_d, proxy_closed * 100, "k--", lw=1.5,
         label="stage-2 proxy  $1-(U_{eff}/U)^2$")
for foil in AIRFOILS:
    ax2.plot(x_over_d, results[foil]['loss3_closed'] * 100, "o-", color=colors[foil],
             label=f"stage-3  {foil}")
ax2.set_xlabel("following distance   x / d")
ax2.set_ylabel("downforce loss   (%)")
ax2.set_title("downforce lost to dirty air (DRS closed)")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
out = FIGURES / "05_xfoil_downforce.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"\nfigure written to {out}")
