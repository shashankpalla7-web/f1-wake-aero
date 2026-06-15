# Derivation 1 — The Far-Wake Velocity Profile

**Goal:** derive a self-consistent, dimensionally sound model for the streamwise velocity
deficit `u'(x, y) = U∞ − u(x, y)` behind a bluff body (the leading car), as a function of
downstream distance `x` and lateral position `y`.

This is the inflow condition fed to XFOIL for the following car's wing (Step 2+).

---

## 1. Shape: why Gaussian

Start from the 2D boundary-layer form of the streamwise momentum equation in the wake,
where pressure gradients vanish far downstream and the deficit is small (`u' ≪ U∞`):

```
U∞ · ∂u'/∂x  =  ν_t · ∂²u'/∂y²
```

`ν_t` here is an **eddy (turbulent) diffusivity**, not the molecular viscosity. This is the
standard closure for turbulent free-shear flows (mixing-length / eddy-viscosity models):
replace molecular ν with an empirical ν_t and the laminar equations carry over formally
unchanged. This equation is the heat equation with `x` playing the role of time.

Seeking a self-similar solution `u'(x,y) = Δu(x) · f(η)`, `η = y/δ(x)`, the heat-equation
analogy gives the point-source solution shape:

```
f(η) = exp(−ln2 · η²)
```

The `ln2` is a normalization choice, not a free parameter — it makes `δ(x)` the
**half-width**: `f(±1) = 1/2`, i.e. `u'(x, ±δ) = Δu(x)/2`. This is the standard
"Gaussian wake" form used in 2D turbulent wake theory (Schlichting, *Boundary-Layer
Theory*, wake chapter).

So:

```
u'(x, y) = Δu(x) · exp( −ln2 · (y/δ(x))² )
```

Two unknown functions remain: the **centerline deficit** `Δu(x)` and the **half-width**
`δ(x)`. Both are fixed below.

---

## 2. Amplitude: momentum conservation fixes `Δu(x)·δ(x)`

Drag on the leading body equals the rate of momentum deficit carried downstream
(far-wake, `u' ≪ U∞`, so `u(U∞−u) ≈ U∞u'`):

```
D' = ρ ∫ u(U∞ − u) dy  ≈  ρ U∞ ∫ u'(x,y) dy
```

This integral is the same at *every* downstream station `x` — drag is constant,
momentum deficit is conserved as the wake spreads and decays. Substituting the
Gaussian profile:

```
∫ u'(x,y) dy = Δu(x) · δ(x) · √(π / ln2)
```

so:

```
D' = ρ U∞ · √(π/ln2) · Δu(x)·δ(x)   =  const
  ⟹  Δu(x) · δ(x) = const  (independent of x)
```

Writing drag per unit span in terms of a drag coefficient and reference length
(`D' = ½ ρ U∞² C_D d`):

```
Δu(x)/U∞  =  (C_D · d / 2) · √(ln2/π)  /  δ(x)
```

**This is the key relation.** It says: once you know how the wake *spreads*
(`δ(x)`), the centerline deficit is fully determined by `C_D`, `d`, and that
spreading — no extra free amplitude parameter.

---

## 3. Spreading: `δ(x) ∝ √x` (empirical, stated as a model assumption)

The heat-equation analogy (§1) predicts `δ(x) ∝ √x` — this exponent is robust and
well-supported for turbulent wakes. The *prefactor* is not universal; it depends on
the body's turbulence characteristics and isn't something a 2D analytical model can
derive from first principles alone.

**Model choice (stated explicitly, this is a limiting assumption of the paper):**

```
δ(x) = δ_d · √(x / d)
```

where `δ_d` is the wake half-width at one reference length downstream (`x = d`),
treated as a **model parameter**. A reasonable starting value is `δ_d ≈ d`
(half-width comparable to the body size just behind it), refined later against
published CFD/wind-tunnel wake data if available.

---

## 4. Complete model

```
δ(x)     = δ_d · √(x/d)
Δu(x)/U∞ = (C_D·d/2)·√(ln2/π) / δ(x)
u'(x,y)  = Δu(x) · exp(−ln2·(y/δ(x))²)
u(x,y)   = U∞ − u'(x,y)
```

Free parameters: `C_D` (leading car's drag coefficient — two values, DRS open/closed),
`d` (reference length), `δ_d` (wake spreading parameter).

## 5. Validation checks (implemented in `scripts/validate_wake_profile.py`)

1. **Scaling check** — `Δu(x)/U∞` plotted vs `x/d` on log-log axes should be a straight
   line of slope **−1/2**, and `δ(x)` of slope **+1/2**.
2. **Momentum conservation** — `∫u'(x,y) dy` (computed numerically by integrating the
   profile) should be **constant across x**, independent of the spreading model. This
   checks the implementation is self-consistent with §2, not just curve-fitted.
3. **Profile shapes** — plot `u'(x,y)/Δu(x)` vs `y/δ(x)` at several `x/d`; all curves
   should collapse onto the single Gaussian `exp(−ln2·η²)` (self-similarity).

## Note on the original concept doc

The concept doc's §1.4 "complete wake solution" used molecular `ν` directly in the
final formula. That's not dimensionally appropriate for a turbulent wake — using
molecular viscosity would predict a wake far too narrow to be physical. The model
above keeps the doc's correct results (Gaussian shape, `x^±1/2` scaling, the role of
`C_D` and `d`) but replaces the molecular-viscosity amplitude with the momentum-integral
relation (§2) and makes the spreading rate an explicit, stated model parameter (§3)
rather than something derived from `ν`. This is a strictly more defensible position for
the paper's "Limitations and Simplifying Assumptions" section.
