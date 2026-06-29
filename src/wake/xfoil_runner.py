# headless XFOIL driver for stage 3
#
# drives the XFOIL 6.99 binary as a subprocess to get the viscous lift
# coefficient of a single-element section at a given Reynolds number. stage 3
# feeds it Re = Re_clean * (U_eff/U_inf) so the wake's effect on C_l (not just
# on dynamic pressure) is captured -- the part stage 2's constant-C_l proxy drops.
#
# two implementation notes that cost real time to find (see derivation 03):
#  1. the prebuilt binary was compiled with -ffpe-trap=invalid,zero -finit-real=inf,
#     which crashes (SIGFPE) on the exact float ops XFOIL's BL solver relies on.
#     fix: rebuild without those flags (keep -fallow-argument-mismatch).
#  2. XFOIL's PACC polar-save file clobbers stdin in this build -- after the first
#     accumulated point the process hits EOF. so we DON'T use PACC; we step alpha
#     up from 0 (each ALFA seeds from the last converged BL solution) and parse the
#     converged C_l / C_d straight from stdout.

import os
import re
import subprocess
from pathlib import Path

XFOIL_BIN = os.environ.get("XFOIL_BIN", str(Path.home() / "Xfoil" / "bin" / "xfoil"))

# one Newton iteration prints across three consecutive lines:
#     <n>   rms:  <rms>   max: <max>   <flag> at ...
#        a =  <alpha>      CL =  <cl>
#       Cm = <cm>     CD =  <cd>   =>   ...
# so alpha+CL share a line and CD is on the next -- parse them separately and pair.
_RMS = re.compile(r"rms:\s*([-\d.E+]+)")
_ACL = re.compile(r"a =\s*(-?\d+\.\d+)\s+CL =\s*(-?\d+\.\d+)")
_CD = re.compile(r"CD =\s*(-?\d+\.\d+)")
_FAILED = re.compile(r"Convergence failed", re.I)


class XfoilError(RuntimeError):
    pass


def _airfoil_cmds(airfoil):
    # accept either a NACA 4/5-digit designation ("NACA 2412" or "2412") or a
    # path to a coordinate .dat file (loaded with LOAD).
    s = str(airfoil).strip()
    if s.upper().startswith("NACA") or (s.isdigit() and len(s) in (4, 5)):
        digits = s.upper().replace("NACA", "").strip()
        return [f"NACA {digits}"]
    path = Path(s)
    if not path.exists():
        raise XfoilError(f"airfoil not a NACA code or existing .dat file: {airfoil!r}")
    return [f"LOAD {path.resolve()}", ""]  # blank line accepts the parsed-name prompt


def run_alpha_ramp(airfoil, Re, alpha_target, *, step=1.0,
                   n_iter=300, timeout=120):
    """ramp angle of attack from 0 up to alpha_target and return per-alpha results.

    returns {alpha_rounded: {"CL": float, "CD": float, "converged": bool}}.
    stepping from 0 lets each viscous point initialise from the previous converged
    boundary layer, which is what keeps the solver off the cold-start cliff.
    """
    if alpha_target < 0:
        raise XfoilError("ramp expects alpha_target >= 0 (run sections upright; "
                         "downforce is |lift|, and the C_l ratio is sign-independent)")

    n = max(1, round(alpha_target / step))
    alphas = [round(i * alpha_target / n, 4) for i in range(n + 1)]  # 0 .. target

    cmds = ["PLOP", "G F", ""]                 # disable graphics
    cmds += _airfoil_cmds(airfoil)
    # Ncrit defaults to 9.0 (free transition), which is the value we want and state;
    # we don't touch VPAR -- its inline-arg prompts are fragile in this build.
    cmds += ["OPER", f"VISC {Re:.0f}", f"ITER {n_iter}"]
    cmds += [f"ALFA {a:.4f}" for a in alphas]
    cmds += ["", "QUIT"]
    stdin = "\n".join(cmds) + "\n"

    try:
        proc = subprocess.run([XFOIL_BIN], input=stdin, capture_output=True,
                              text=True, timeout=timeout, cwd="/tmp")
    except subprocess.TimeoutExpired:
        raise XfoilError(f"XFOIL timed out after {timeout}s (airfoil={airfoil}, Re={Re:.0f})")
    log = proc.stdout + proc.stderr
    if "SIGFPE" in log:
        raise XfoilError("XFOIL crashed with SIGFPE -- binary still has the FPE trap; "
                         "rebuild without -ffpe-trap/-finit-real (see derivation 03).")

    # walk the log once, splitting iterations by the commanded-alpha boundaries.
    # the converged value for an alpha is its last iteration; convergence is taken
    # from that iteration's rms residual (XFOIL stops iterating once it's tiny).
    out = {}
    last_rms = {}
    cur_rms = None
    pending = None  # (alpha, cl) awaiting its CD line
    for line in log.splitlines():
        m_rms = _RMS.search(line)
        if m_rms:
            try:
                cur_rms = abs(float(m_rms.group(1)))
            except ValueError:
                cur_rms = None
            continue
        m_acl = _ACL.search(line)
        if m_acl:
            pending = (round(float(m_acl.group(1)), 3), float(m_acl.group(2)))
            continue
        m_cd = _CD.search(line)
        if m_cd and pending is not None:
            al, cl = pending
            out[al] = {"CL": cl, "CD": float(m_cd.group(1))}
            if cur_rms is not None:
                last_rms[al] = cur_rms
            pending = None
    for al, rec in out.items():
        rec["converged"] = last_rms.get(al, 1.0) < 5e-3

    if not out:
        raise XfoilError(f"XFOIL produced no converged points (airfoil={airfoil}, Re={Re:.0f})")
    return out


def cl_at(airfoil, Re, alpha, **kw):
    """converged viscous C_l at one (airfoil, Re, alpha); None if it didn't converge."""
    res = run_alpha_ramp(airfoil, Re, alpha, **kw)
    target = round(alpha, 3)
    rec = res.get(target) or res.get(min(res, key=lambda a: abs(a - alpha)))
    if rec is None or not rec["converged"]:
        return None
    return rec["CL"]
