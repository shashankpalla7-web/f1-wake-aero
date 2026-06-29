# authentic XFOIL plot dumps -- proof-of-modelling screenshots for the repo
#
# for each airfoil, runs XFOIL to the operating point and dumps its native
# pressure-distribution plot (Cp vs x/c with the section and the full force
# header: Re, alpha, C_L, C_M, C_D, L/D, Ncrit) straight from the solver via the
# HARD command. no manual screenshotting, no interactive window -- one pass per
# airfoil, identical every time.
#
# pipeline per airfoil:
#   XFOIL  -- ramp alpha to the operating point, HARD -> plot.ps  (graphics on so
#             the Cp plot buffer is populated; runs headless against XQuartz)
#   gs     -- PostScript -> PNG
#   PIL    -- rotate to landscape (XFOIL renders sideways on the page) + autocrop
#
# output: figures/xfoil/cp_naca<dddd>.png  (tracked in git, unlike the regenerable
# figures/*.png).  needs ghostscript on PATH (brew install ghostscript).

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageChops

XFOIL = os.environ.get("XFOIL_BIN", str(Path.home() / "Xfoil" / "bin" / "xfoil"))
GS = shutil.which("gs")
OUT = Path(__file__).resolve().parents[1] / "figures" / "xfoil"

RE_CLEAN = 1_000_000     # clean-air operating Reynolds number (matches stage 3)
ALPHA = 5.0              # operating angle of attack, deg
AIRFOILS = ["2412", "4412", "6412"]   # the stage-3 camber ladder


def dump_plot_ps(digits, workdir):
    # ramp alpha from 0 to the operating point (each point seeds the next), then
    # HARD dumps the current Cp plot to plot.ps in the working directory.
    alphas = range(0, int(ALPHA) + 1)
    cmds = [f"NACA {digits}", "OPER", f"VISC {RE_CLEAN}", "ITER 250"]
    cmds += [f"ALFA {a}" for a in alphas]
    cmds += ["HARD", "QUIT"]
    proc = subprocess.run([XFOIL], input="\n".join(cmds) + "\n", text=True,
                          capture_output=True, timeout=120, cwd=workdir)
    if "SIGFPE" in proc.stdout + proc.stderr:
        raise RuntimeError("XFOIL SIGFPE -- rebuild without the FPE trap (derivation 03)")
    ps = Path(workdir) / "plot.ps"
    if not ps.exists():
        raise RuntimeError(f"NACA {digits}: XFOIL did not write plot.ps")
    return ps


def ps_to_png(ps, png, dpi=200, pad=24):
    raw = png.with_suffix(".raw.png")
    subprocess.run([GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
                    f"-r{dpi}", f"-sOutputFile={raw}", str(ps)], check=True)
    im = Image.open(raw).convert("RGB").rotate(-90, expand=True)  # page -> landscape
    bbox = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
    if bbox:
        bbox = (max(0, bbox[0] - pad), max(0, bbox[1] - pad),
                min(im.width, bbox[2] + pad), min(im.height, bbox[3] + pad))
        im = im.crop(bbox)
    im.save(png)
    raw.unlink(missing_ok=True)


def main():
    if GS is None:
        raise SystemExit("ghostscript not found on PATH -- run: brew install ghostscript")
    OUT.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    print(f"XFOIL Cp plots at Re={RE_CLEAN:.0e}, alpha={ALPHA:.0f} deg -> {OUT.relative_to(root)}/")
    for d in AIRFOILS:
        with tempfile.TemporaryDirectory() as wd:
            ps = dump_plot_ps(d, wd)
            out = OUT / f"cp_naca{d}.png"
            ps_to_png(ps, out)
            print(f"  NACA {d}  ->  {out.relative_to(root)}")


if __name__ == "__main__":
    main()
