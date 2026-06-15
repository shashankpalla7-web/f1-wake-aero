# Derivation 1 - far wake velocity profile

Goal: get a velocity deficit profile u'(x,y) = U∞ - u(x,y) behind the leading car that's
actually physically defensible, as a function of downstream distance x and lateral
position y. This is what gets fed into XFOIL as the inflow for the following car's wing
later on.

---

## 1. Why the profile is Gaussian

Start from the boundary layer form of the streamwise momentum eq in the wake (far
downstream, pressure gradient ~0, deficit small so u' << U∞):

```
U∞ du'/dx = ν_t d²u'/dy²
```

ν_t here is an eddy (turbulent) diffusivity, not the molecular viscosity. This is the
normal closure for turbulent free shear flows - swap molecular ν for an empirical ν_t and
the laminar equations carry over basically unchanged. This is just the heat equation with
x playing the role of time.

Looking for a self similar solution u'(x,y) = Δu(x) f(η), η = y/δ(x), the heat equation
analogy gives the point source shape:

```
f(η) = exp(-ln2 · η²)
```

(the ln2 isn't really a free parameter - it's there so that δ(x) ends up being the
half-width, f(±1) = 1/2, so u'(x,±δ) = Δu(x)/2). This Gaussian wake form is all over 2D
turbulent wake theory.

So:

```
u'(x,y) = Δu(x) exp(-ln2 (y/δ(x))²)
```

Two things still unknown - Δu(x) and δ(x). Both get pinned down below.

---

## 2. Amplitude - momentum conservation gives Δu(x)δ(x)

Drag on the leading car = rate of momentum deficit carried downstream (far wake,
u' << U∞, so u(U∞-u) ≈ U∞u'):

```
D' = ρ ∫ u(U∞-u) dy ≈ ρU∞ ∫ u'(x,y) dy
```

This integral has to be the same at every x - drag doesn't change with x, so the momentum
deficit is conserved as the wake spreads out and weakens. Plug in the Gaussian:

```
∫ u'(x,y) dy = Δu(x) δ(x) √(π/ln2)
```

so D' = ρU∞ √(π/ln2) Δu(x)δ(x) = const, which means

```
Δu(x) δ(x) = const   (no x dependence)
```

Writing drag per unit span as D' = ½ρU∞²C_D d:

```
Δu(x)/U∞ = (C_D d / 2) √(ln2/π) / δ(x)
```

So once δ(x) is known, Δu(x) is fully determined by C_D and d - no extra amplitude
parameter floating around.

---

## 3. Spreading - δ(x) ~ √x

The heat equation analogy says δ(x) ~ √x. The exponent is solid and shows up in basically
every turbulent wake paper, but the prefactor isn't universal, it depends on the body's
turbulence and a 2D analytical model can't really get it from scratch.

So this is just stated as a model parameter:

```
δ(x) = δ_d √(x/d)
```

where δ_d is the half-width at x = d. Starting guess δ_d ≈ d (half-width about the size of
the car right behind it) - can come back and tune this against real wake data if I find
some.

---

## 4. Putting it together

```
δ(x)      = δ_d √(x/d)
Δu(x)/U∞  = (C_D d / 2) √(ln2/π) / δ(x)
u'(x,y)   = Δu(x) exp(-ln2 (y/δ(x))²)
u(x,y)    = U∞ - u'(x,y)
```

Free params: C_D (drag coefficient, two values for DRS open/closed), d (reference
length), δ_d (spreading parameter).

## 5. Checks (scripts/validate_wake_profile.py)

1. **Scaling** - Δu/U∞ vs x/d on log-log should be a straight line of slope -1/2, and
   δ(x) should be +1/2.
2. **Momentum conservation** - ∫u'(x,y)dy should come out constant across x. This checks
   the code actually matches part 2 and isn't just curve fit to look right.
3. **Profile shape** - u'(x,y)/Δu(x) vs y/δ(x) at a few different x/d should all land on
   the same Gaussian curve.

## Note on the concept doc

The original concept doc's wake formula used molecular ν directly in the final result,
which doesn't really work for a turbulent wake - you'd get a wake way too narrow to be
real. Kept the parts that are right (Gaussian shape, x^±1/2 scaling, C_D and d showing
up) but swapped the amplitude for the momentum integral relation from part 2, and made
the spreading rate an explicit stated parameter instead of something derived from ν.
Worth a line about this in the limitations section of the paper.
