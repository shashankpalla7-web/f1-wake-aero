# effective inflow, basically the speed seen by the car
# allowing us to calculate X% velocity loss into y% downforce loss

import numpy as np
from scipy.integrate import quad

from .wake_model import WakeParameters, half_width, centerline_deficit, velocity_profile

LN2 = np.log(2.0)

# function that efectively inverts the wake model to give an inflow speed, which we can then use to calculate downforce loss.
def centerline_inflow(x, params: WakeParameters):
    return params.U_inf - centerline_deficit(x, params)


def span_mean_inflow(x, params: WakeParameters, span):
    from scipy.special import erf
    delta = half_width(x, params)
    du = centerline_deficit(x, params)
    a = (span / (2.0 * delta)) * np.sqrt(LN2)
    mean_deficit = du * (delta / span) * np.sqrt(np.pi / LN2) * erf(a)
    return params.U_inf - mean_deficit


def rms_inflow(x, params: WakeParameters, span):
    # dynamic-pressure-equivalent inflow: the uniform speed delivering the same
    # this is the headline U_eff fed to XFOIL for simulations
    u_sq = lambda y: velocity_profile(x, y, params) ** 2
    mean_u_sq, _ = quad(u_sq, -span / 2.0, span / 2.0)
    mean_u_sq /= span
    return np.sqrt(mean_u_sq)


def effective_inflow(x, params: WakeParameters, span, method="rms"):
    if method == "rms":
        return rms_inflow(x, params, span)
    if method == "mean":
        return span_mean_inflow(x, params, span)
    if method == "centerline":
        return centerline_inflow(x, params)
    raise ValueError(f"unknown method: {method!r}")

#how much downforce we're losing ( the actual data collection function)
def downforce_loss_fraction(x, params: WakeParameters, span, method="rms"):
    u_eff = effective_inflow(x, params, span, method)
    return 1.0 - (u_eff / params.U_inf) ** 2
