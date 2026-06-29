from .wake_model import WakeParameters, half_width, centerline_deficit, velocity_profile, momentum_deficit
from .inflow import (
    centerline_inflow,
    span_mean_inflow,
    rms_inflow,
    effective_inflow,
    downforce_loss_fraction,
)
from .xfoil_runner import run_alpha_ramp, cl_at, XfoilError

__all__ = [
    "WakeParameters",
    "half_width",
    "centerline_deficit",
    "velocity_profile",
    "momentum_deficit",
    "centerline_inflow",
    "span_mean_inflow",
    "rms_inflow",
    "effective_inflow",
    "downforce_loss_fraction",
    "run_alpha_ramp",
    "cl_at",
    "XfoilError",
]


# initiallizng all parameters and functions at the package level for easy import in scripts -  from wake import WakeParameters, centerline_deficit
