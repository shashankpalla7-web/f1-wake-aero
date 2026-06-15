# far wake velocity profile model
# derivation: derivations/01_wake_profile_derivation.md

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

LN2 = np.log(2.0)


@dataclass(frozen=True)
class WakeParameters:
    U_inf: float    # freestream velocity, m/s
    C_D: float      # leading car's drag coeff (DRS open or closed)
    d: float        # reference length, m
    delta_d: float  # wake half-width at x = d (spreading param, derivation pt 3)


def half_width(x, params: WakeParameters):
    # delta(x) = delta_d * sqrt(x/d)
    return params.delta_d * np.sqrt(x / params.d)


def centerline_deficit(x, params: WakeParameters):
    # from momentum conservation, derivation pt 2
    delta = half_width(x, params)
    return params.U_inf * (params.C_D * params.d / 2.0) * np.sqrt(LN2 / np.pi) / delta


def velocity_profile(x, y, params: WakeParameters):
    delta = half_width(x, params)
    du = centerline_deficit(x, params)
    return params.U_inf - du * np.exp(-LN2 * (y / delta) ** 2)


def momentum_deficit(x, params: WakeParameters):
    # integral of u'(x,y) dy -- should be constant across x.
    # this is the conservation law from pt 2, used here as a check on the
    # implementation rather than as an input
    deficit = lambda y: params.U_inf - velocity_profile(x, y, params)
    width = 10.0 * half_width(x, params)  # gaussian tails are ~0 past this
    value, _ = quad(deficit, -width, width)
    return value
