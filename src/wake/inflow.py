# effective inflow seen by the following car's wing + downforce-loss proxy
# derivation: derivations/02_effective_inflow_derivation.md

import numpy as np
from scipy.integrate import quad

from .wake_model import WakeParameters, half_width, centerline_deficit, velocity_profile

LN2 = np.log(2.0)


def centerline_inflow(x, params: WakeParameters):
    # U_c(x) = u(x,0): deepest point of the deficit (derivation 2, 1a)
    return params.U_inf - centerline_deficit(x, params)


def span_mean_inflow(x, params: WakeParameters, span):
    # arithmetic span-average of u(x,y) over y in [-b/2, b/2] (derivation 2, 1b)
    # closed form: U_inf - Delta_u (delta/b) sqrt(pi/ln2) erf((b/2delta) sqrt(ln2))
    from scipy.special import erf
    delta = half_width(x, params)
    du = centerline_deficit(x, params)
    a = (span / (2.0 * delta)) * np.sqrt(LN2)
    mean_deficit = du * (delta / span) * np.sqrt(np.pi / LN2) * erf(a)
    return params.U_inf - mean_deficit


def rms_inflow(x, params: WakeParameters, span):
    # dynamic-pressure-equivalent inflow: the uniform speed delivering the same
    # span-integrated q = 1/2 rho u^2 as the real profile (derivation 2, 1c).
    # this is the headline U_eff fed to XFOIL in stage 3.
    u_sq = lambda y: velocity_profile(x, y, params) ** 2
    mean_u_sq, _ = quad(u_sq, -span / 2.0, span / 2.0)
    mean_u_sq /= span
    return np.sqrt(mean_u_sq)


def effective_inflow(x, params: WakeParameters, span, method="rms"):
    # representative inflow speed U_eff(x). method: "rms" (headline), "mean", "centerline"
    if method == "rms":
        return rms_inflow(x, params, span)
    if method == "mean":
        return span_mean_inflow(x, params, span)
    if method == "centerline":
        return centerline_inflow(x, params)
    raise ValueError(f"unknown method: {method!r}")


def downforce_loss_fraction(x, params: WakeParameters, span, method="rms"):
    # fraction of clean-air downforce lost to the wake at distance x.
    # downforce ~ u^2 at fixed C_l, so dL/L0 = 1 - (U_eff/U_inf)^2 (derivation 2, 2)
    u_eff = effective_inflow(x, params, span, method)
    return 1.0 - (u_eff / params.U_inf) ** 2
