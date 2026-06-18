# Derivation 2 - effective inflow U_eff(x) and the downforce-loss proxy

Goal: take the far-wake velocity field u(x,y) from derivation 1 and reduce it to a single
effective inflow speed U_eff(x) that the *following* car's wing actually sees at a following
distance x. Then convert that into the quantity the paper is really about: how much downforce
the following car loses to the leading car's dirty air, as a function of x.

This is the bridge between the wake physics (derivation 1) and the XFOIL stage (stage 3).
Stage 2 isolates the *dynamic-pressure* effect of the wake (how much slower the air is);
stage 3 will refine the *C_l* effect (how the wing's lift coefficient responds to that inflow).

---

## 0. What we're starting from

From derivation 1, the wake behind the leading car is

```
u(x,y)  = U∞ - u'(x,y)
u'(x,y) = Δu(x) exp(-ln2 (y/δ(x))²)
δ(x)    = δ_d √(x/d)
Δu(x)/U∞ = (C_D d / 2) √(ln2/π) / δ(x)
```

x = downstream (following) distance, y = lateral position across the track, y = 0 the wake
centerline.

---

## 1. The reduction problem - a field has to become a number

A wing doesn't see one speed. Across its span it sees the whole profile u(x,y): slowest on
the centerline, recovering toward U∞ at the tips. XFOIL (stage 3) takes a single freestream
speed, so the profile has to be collapsed to one representative number U_eff(x).

Model the following car's front wing as a flat span of width b centered on the wake centerline,
occupying y ∈ [-b/2, b/2]. Baseline b = 2.0 m (FIA front-wing regulation max width, 2000 mm);
centered-on-wake is the straight-line, in-line-following assumption. Both are stated limitations
(§5).

There are three defensible ways to reduce the profile, and - same posture as the C_D sourcing -
all three are reported so the sensitivity to the choice is explicit rather than hidden.

### 1a. Centerline
```
U_c(x) = u(x,0) = U∞ - Δu(x)
```
The deepest point of the deficit. Pessimistic - ignores that the wing extends out to where
the wake is weaker. A lower bound on U_eff.

### 1b. Arithmetic span-mean
```
Ū(x) = (1/b) ∫_{-b/2}^{b/2} u(x,y) dy
```
With the Gaussian, the integral is a clean error function:
```
∫_{-b/2}^{b/2} exp(-ln2 (y/δ)²) dy = δ √(π/ln2) · erf( (b/2δ)√ln2 )
```
so
```
Ū(x) = U∞ - Δu(x) (δ/b) √(π/ln2) erf( (b/2δ)√ln2 )
```
The "average speed the span sees." Natural, but not the right thing for a downforce claim -
see 1c.

### 1c. RMS / dynamic-pressure-equivalent  ← headline
The wing makes downforce, and downforce scales with dynamic pressure q = ½ρu², not with u.
Under strip theory (§3) the wing's total downforce ∝ ∫ u(y)² dy across the span. The single
uniform speed that delivers the *same span-integrated dynamic pressure* is the root-mean-square
of the profile:
```
U_eff(x) ≡ U_rms(x) = sqrt( (1/b) ∫_{-b/2}^{b/2} u(x,y)² dy )
```
This is the physically correct reduction for a downforce statement, so U_eff ≡ U_rms is the
headline definition fed to XFOIL in stage 3. Closed form (expanding (U∞ - u')²):
```
(1/b)∫ u² dy = U∞²
              - (2 U∞ Δu δ / b) √(π/ln2)   erf( (b/2δ)√ln2 )
              + (Δu² δ / b)   √(π/2ln2)    erf( (b/2δ)√2ln2 )
```
(the code integrates numerically and the closed form is the cross-check.)

Ordering: U_c ≤ Ū ≤ U_rms ≤ U∞. Centerline gives the biggest loss; Ū and U_rms sit close
together when the inflow is fairly uniform across the span (low variance), and the gap between
them *is* the spanwise non-uniformity of the inflow - itself a real, separate degrader of wing
performance (a wing in sheared inflow underperforms even at matched mean q). Noted, not modeled
here.

---

## 2. From U_eff to the downforce loss

At fixed wing geometry and fixed lift coefficient C_l, sectional downforce ∝ u². So the
following car's downforce relative to the same wing in clean air is

```
L(x) / L_0 = (U_eff(x) / U∞)²        with U_eff = U_rms
```

and the **downforce-loss fraction** - the paper's headline result in proxy form - is

```
ΔL/L_0 (x) = 1 - (U_eff(x)/U∞)²
```

Equivalent, and the more intuitive statement, using U_rms² = mean of u² across the span:

```
L(x)/L_0 = (1/b) ∫_{-b/2}^{b/2} (u(x,y)/U∞)² dy   =  ⟨(u/U∞)²⟩
```

i.e. the downforce ratio is just the span-average of the squared speed ratio.

### The 2× rule of thumb
For a small uniform-ish deficit ε = 1 - U_eff/U∞,
```
ΔL/L_0 = 1 - (1-ε)² = 2ε - ε² ≈ 2ε.
```
So **the downforce penalty is roughly twice the speed deficit.** A 17% effective-speed deficit
costs ~31% of downforce, not 17%. This quadratic amplification is the core reason dirty air
hurts following cars so much more than the raw "air is X% slower" number suggests - a clean,
reviewer-facing takeaway.

---

## 3. Strip-theory assumption behind L ∝ ∫u²

Total downforce on a span-b wing in a non-uniform inflow u(y):
```
L = ∫_{-b/2}^{b/2} ½ ρ u(y)² c C_l(y) dy
```
Assume chord c and sectional lift coefficient C_l are constant across the span - rectangular
wing, and the mild spanwise speed gradient doesn't change the local C_l. Then
```
L = ½ ρ c C_l ∫ u(y)² dy,   L_0 = ½ ρ c C_l U∞² b,
L/L_0 = (1/b)∫(u/U∞)² dy,
```
which is exactly §2. The constant-C_l assumption is what stage 3 (XFOIL) relaxes: it recomputes
C_l at the actual inflow so the *shape*/angle-of-attack response is captured, not just the
dynamic-pressure scaling isolated here.

---

## 4. Putting it together

```
δ(x)        = δ_d √(x/d)
Δu(x)/U∞    = (C_D d / 2) √(ln2/π) / δ(x)
U_c(x)      = U∞ - Δu(x)
Ū(x)        = U∞ - Δu(x)(δ/b)√(π/ln2) erf((b/2δ)√ln2)
U_eff(x)    = U_rms(x) = sqrt( (1/b)∫_{-b/2}^{b/2} u(x,y)² dy )      ← headline
ΔL/L_0 (x)  = 1 - (U_eff(x)/U∞)²
```
New free parameter introduced this stage: b (following wing span), baseline 2.0 m. Everything
else (C_D, d, δ_d) carries over from derivation 1.

## 5. Checks (scripts/effective_inflow.py)

1. **Limits.** b → 0 ⇒ all three reductions → centerline u(x,0). x → ∞ ⇒ all → U∞ and ΔL/L_0 → 0.
2. **Ordering.** U_c ≤ Ū ≤ U_rms ≤ U∞ at every x.
3. **Closed form.** Numerically integrated U_rms² matches the erf closed form in §1c.
4. **2× rule.** ΔL/L_0 ≈ 2(1 - U_eff/U∞) at large x where the deficit is small.

## 6. Limitations to carry into the writeup

- **Wing centered on the wake.** Real following cars run an offset racing line to escape the
  wake; off-center would see a shallower, asymmetric profile. Baseline is the worst-case in-line
  follow.
- **Constant C_l (strip theory).** Spanwise pressure equalization and 3D/induced effects ignored;
  stage 3 addresses the C_l response, not the spanwise coupling.
- **Single transverse coordinate.** Derivation 1's wake is 2D in (x,y); averaging over the lateral
  span only. The vertical extent of the real wake and ground effect are out of scope here.
- **Span value b = 2.0 m** is the regulation max; effective aerodynamic span is smaller. b enters
  only weakly (through erf) once the wing is inside the wake, but it's a stated input, sweepable.
```