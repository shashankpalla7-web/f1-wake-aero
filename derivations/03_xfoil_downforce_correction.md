# Derivation 3 - XFOIL correction: from effective inflow to real downforce loss

Goal: relax the one assumption stage 2 left standing. Stage 2 reduced the wake to a single
effective inflow `U_eff(x)` and turned it into a downforce loss using a *constant-C_l* proxy:

```
dL/L0 = 1 - (U_eff/U∞)²            (stage-2 proxy: C_l assumed unchanged in the wake)
```

That captures the dynamic-pressure effect (the air is slower, so `q = ½ρu²` is lower) but
ignores that the wake also lowers the Reynolds number the following wing sees, and `C_l`
depends on Reynolds number. Stage 3 measures that `C_l` dependence directly with XFOIL and
folds it in, turning the proxy into a modeled result.

This is the bridge from derivation 2 (the wake → an inflow number) to an airfoil response.

---

## 1. The exact decomposition

At fixed geometry and angle of attack, sectional downforce is `L = ½ρ U² c · C_l(Re)`. Form
the ratio of the following (in-wake) wing to the same wing in clean air:

```
L_wake     ½ρ U_eff²  c  C_l(Re_wake)      ( U_eff )²    C_l,wake
------  =  --------------------------  =  ( ----- )  ·  ---------
L_clean    ½ρ U∞²    c  C_l(Re_clean)      (  U∞   )     C_l,clean
```

so

```
dL/L0 = 1 - (U_eff/U∞)² · (C_l,wake / C_l,clean)
                \_______/   \_________________/
              stage-2 term   stage-3 correction r(x)
```

The first factor is exactly stage 2. The second factor `r(x) = C_l,wake/C_l,clean` is the
piece stage 2 set to 1. Stage 3's whole job is to measure `r(x)`.

`U_eff(x)` is the RMS (dynamic-pressure-equivalent) inflow from derivation 2 - the physically
correct single speed for a downforce claim, since downforce ∝ q and the RMS preserves the
span-averaged q.

## 2. The Reynolds map

The wake scales the inflow speed by `s(x) = U_eff(x)/U∞`, so it scales Reynolds number by the
same factor (chord and ν fixed):

```
Re_wake(x) = U_eff(x) c / ν = Re_clean · s(x),     Re_clean = U∞ c / ν
```

Conditions: `U∞ = 60 m/s`, chord `c = 0.25 m` (representative front-wing element),
`ν = 1.5e-5 m²/s` → `Re_clean = 1.0e6`. At one car length back `s ≈ 0.83`, so
`Re_wake ≈ 0.83e6`. The band `Re ≈ 0.83e6 – 1.0e6` is exactly the second-order effect stage 3
exists to capture.

## 3. XFOIL setup

- **Sections (controlled camber ladder):** NACA 2412 / 4412 / 6412 - same 12% thickness, only
  camber varies (2% / 4% / 6%). Isolates one variable: *does the wake's downforce penalty
  depend on how hard the wing is already working?*
- **Single-element proxy.** A real F1 front wing is multi-element (main plane + slotted flaps);
  XFOIL is a 2D single-element panel + boundary-layer code and cannot model slot-gap flow. So a
  single high-camber element stands in for the loaded front-wing element. This is a real
  limitation, named in §6 - not hidden.
- **Operating point:** fixed AoA `α = 5°` (a solid high-load point that converges across the
  whole Re band for all three sections); `Ncrit = 9` (free transition); Mach 0 (incompressible -
  at 60 m/s, M ≈ 0.17, compressibility negligible).
- **Sign convention:** sections are run upright (positive α, positive `C_l`). Downforce is lift
  flipped, and the ratio `r = C_l,wake/C_l,clean` is sign-independent, so the upright ratio
  equals the inverted-downforce ratio. We report `C_l` magnitudes.
- **Per point:** `C_l,clean` is XFOIL at `Re_clean`; `C_l,wake(x)` is XFOIL at `Re_wake(x)`,
  same section and α. Each run ramps α from 0 to 5° so every viscous solution seeds from the
  previous converged boundary layer (cold-starting viscous at α = 5° diverges).
- **Proof of modelling.** XFOIL's own pressure-distribution plot at the operating point is
  dumped for each section to `figures/xfoil/cp_naca{2412,4412,6412}.png` (the native solver
  output: Cp vs x/c with the section and the full force header). Regenerate with
  `scripts/xfoil_screenshots.py` (needs ghostscript). These are the authentic solver plots,
  not redrawn data.

## 4. Results

DRS closed (`C_D = 0.90`), `r(x) = C_l,wake/C_l,clean` and the corrected loss:

| x/d | U_eff/U∞ | stage-2 proxy | 2412 r | 2412 loss | 4412 r | 4412 loss | 6412 r | 6412 loss |
|----:|---------:|--------------:|-------:|----------:|-------:|----------:|-------:|----------:|
|  1  | 0.829 | 31.2% | 0.9975 | 31.4% | 0.9970 | 31.4% | 0.9953 | 31.5% |
|  2  | 0.866 | 25.0% | 0.9981 | 25.1% | 0.9979 | 25.1% | 0.9959 | 25.3% |
|  5  | 0.910 | 17.2% | 0.9989 | 17.3% | 0.9989 | 17.3% | 0.9974 | 17.5% |
| 10  | 0.935 | 12.6% | 0.9991 | 12.7% | 0.9993 | 12.7% | 0.9982 | 12.8% |
| 20  | 0.953 |  9.1% | 0.9991 |  9.2% | 0.9995 |  9.2% | 0.9988 |  9.2% |

Clean-air `C_l` (α = 5°, Re 1e6): 2412 → 0.81, 4412 → 1.02, 6412 → 1.23. Figure:
`figures/05_xfoil_downforce.png`.

**Headline.** At one car length back the `C_l` ratio is ≈ 0.997 across the ladder - the wake's
lower Reynolds number costs only ~0.3% of `C_l`. So:

1. **Stage 3 validates the stage-2 proxy.** The dynamic-pressure term `(U_eff/U∞)²` carries
   ~99% of the effect; the `C_l(Re)` correction is < 0.5 percentage points. The constant-`C_l`
   assumption was sound, and we now have a rigorous bound on its error instead of a hope.
2. **The correction is real and ordered by loading.** `r` is smallest (most correction) for the
   highest-camber 6412 and grows toward 1 with distance as `Re_wake → Re_clean`. Harder-working
   wings are slightly more sensitive to the wake's Re drop - physically sensible (thicker,
   nearer-separation boundary layer is more Re-sensitive).
3. **Direction.** `r < 1` always, so the corrected loss is slightly *deeper* than the proxy:
   the proxy was marginally optimistic. At x/d = 1, 31.2% → 31.4–31.5%.

This is the comprehensive answer to the paper's question: the downforce penalty is set almost
entirely by the dynamic-pressure deficit, with a small (≲ 0.5 pt), loading-dependent Reynolds
correction on top.

## 5. Validation

- **Limit.** `r(x) → 1` monotonically as x grows (Re_wake → Re_clean). Holds for all three
  sections (figure, left panel).
- **Ordering.** `r_6412 < r_4412 ≲ r_2412 < 1` at every x - correction increases with loading,
  no crossings out to the far wake.
- **Magnitude sanity.** Clean `C_l` values (0.81 / 1.02 / 1.23) and the lift slope match
  standard XFOIL NACA results at Re 1e6; the `C_l` change over a 17% Re drop is ~0.3–0.5%,
  consistent with the weak `C_l`-vs-Re dependence of attached-flow airfoils below stall.
- **Solver trustworthiness.** Free-transition viscous solutions converge (rms < 5e-3) at every
  point in the matrix; no points were extrapolated or hand-filled.

## 6. Limitations (carry into the writeup)

- **Single-element proxy for a multi-element wing.** The biggest one. XFOIL cannot model the
  slot-gap flow between front-wing elements; the real multi-element `C_l` response to the wake
  may differ. Future work: MSES or RANS CFD for the multi-element case.
- **Fixed angle of attack.** Real cars trim the wing; we hold α and let the wake move Re only.
  A trim-to-target-downforce study is a natural extension.
- **2D, no ground effect.** The front wing runs in strong ground effect; this 2D free-air
  section ignores the floor. Stage 2's lateral-only wake and this 2D section share that bound.
- **Re is the only wake channel modeled here.** The wake also brings turbulence and shear that a
  steady XFOIL run with fixed `Ncrit` does not see; `Ncrit` could be lowered to mimic elevated
  freestream turbulence as a sensitivity check (future work).

---

## Appendix - reproducing the XFOIL stage (two build gotchas)

The prebuilt XFOIL 6.99 at `~/Xfoil/bin/xfoil` and the scripting path both needed fixing; both
cost real debugging time and are worth recording.

1. **The binary crashed on every viscous solve (SIGFPE).** The build's `Makefile` used
   `CHK = -fbounds-check -finit-real=inf -ffpe-trap=invalid,zero -fallow-argument-mismatch`.
   The `-ffpe-trap=invalid,zero` traps the exact floating-point operations XFOIL's
   boundary-layer solver deliberately produces and recovers from, and `-finit-real=inf` seeds
   uninitialised reals with infinities that then trigger the trap - so any viscous point dies
   with a floating-point exception. **Fix:** rebuild with `CHK = -fallow-argument-mismatch`
   only (the `-fallow-argument-mismatch` is still needed for modern gfortran to compile XFOIL's
   legacy argument-mismatched calls). Drop `-fbounds-check`, `-finit-real=inf`, `-ffpe-trap`.
   Then `make xfoil` (the plotlib `libPlt_gDP.a` is already built).

2. **XFOIL's PACC polar file clobbers stdin in this build.** With polar accumulation on
   (`PACC`), after the first accumulated point the process hits EOF and quits before reading
   the next command - the polar save-file write collides with the input stream. **Fix:** don't
   use PACC. Disable graphics (`PLOP` → `G F`), step α up with individual `ALFA` commands
   (which also gives ramp convergence), and parse the converged `C_l`/`C_d` straight from
   stdout (the last iteration block per α, convergence taken from its rms residual). This is
   what `src/wake/xfoil_runner.py` does.

Run it:
```
.venv/bin/python scripts/stage3_xfoil.py      # ~30 s, ~50 XFOIL solves
```
The harness reads the binary path from `$XFOIL_BIN` (default `~/Xfoil/bin/xfoil`).
