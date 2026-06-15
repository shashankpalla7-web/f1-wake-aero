"""Far-wake velocity profile model. See derivations/01_wake_profile_derivation.md."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

LN2 = np.log(2.0)


@dataclass(frozen=True)
class WakeParameters:
    """
    U_inf: freestream velocity [m/s]
    C_D:   leading car's drag coefficient (DRS open or closed)
    d:     reference length [m] (wake width scale of the leading body)
    delta_d: wake half-width at x = d [m] -- the spreading-rate model
             parameter from derivation §3. delta(x) = delta_d * sqrt(x/d).
    """

    U_inf: float
    C_D: float
    d: float
    delta_d: float


def half_width(x, params: WakeParameters):
    """delta(x) = delta_d * sqrt(x/d)  -- derivation §3."""
    return params.delta_d * np.sqrt(x / params.d)


def centerline_deficit(x, params: WakeParameters):
    """Delta_u(x), from momentum conservation -- derivation §2:

        Delta_u(x)/U_inf = (C_D*d/2) * sqrt(ln2/pi) / delta(x)
    """
    delta = half_width(x, params)
    return params.U_inf * (params.C_D * params.d / 2.0) * np.sqrt(LN2 / np.pi) / delta


def velocity_profile(x, y, params: WakeParameters):
    """u(x, y) = U_inf - Delta_u(x) * exp(-ln2 * (y/delta(x))^2) -- derivation §4."""
    delta = half_width(x, params)
    du = centerline_deficit(x, params)
    return params.U_inf - du * np.exp(-LN2 * (y / delta) ** 2)


def momentum_deficit(x, params: WakeParameters):
    """Numerically integrate u'(x,y) = U_inf - u(x,y) over all y.

    Should be constant across x -- this is the conservation law that
    fixed Delta_u(x)*delta(x) in derivation §2, used here as a
    validation check on the implementation rather than as an input.
    """
    deficit = lambda y: params.U_inf - velocity_profile(x, y, params)
    width = 10.0 * half_width(x, params)  # Gaussian tails are negligible beyond this
    value, _ = quad(deficit, -width, width)
    return value
