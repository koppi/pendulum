# Double Inverted Pendulum — Control in the Browser

This simulation brings a double inverted pendulum on a cart from the hanging
rest position into the upright **up-up position** and holds it there. The
controller has three parts:

* **LQR state feedback** about the upper equilibrium, as a genuine solution of
  the algebraic Riccati equation (matrix sign function),
* **swing-up by trajectory optimisation** (iLQR on the full nonlinear model),
  tracked by a **time-varying LQR**,
* **settling** through energy extraction, which puts the system into the start
  state of the planned trajectory before each attempt.

The hand-over point between swing-up and balancing is not a heuristic — it is a
sub-level set of the Riccati matrix.

The entire computation — Riccati solver and trajectory optimisation included —
runs in JavaScript in the browser, without libraries, in real time.

---

## 1. Physical model

### 1.1 Coordinates and state space

The state vector has six components:

```
x = [x, ẋ, θ₁, θ̇₁, θ₂, θ̇₂]ᵀ
```

| Quantity | Meaning | Unit |
|----------|---------|------|
| `x`   | cart position | m |
| `ẋ`   | cart velocity | m/s |
| `θ₁`  | angle of the first link (0 = upright) | rad |
| `θ̇₁`  | angular velocity of the first link | rad/s |
| `θ₂`  | angle of the second link (0 = upright) | rad |
| `θ̇₂`  | angular velocity of the second link | rad/s |
| `θ̈₁`  | angular acceleration of the first link | rad/s² |
| `θ̈₂`  | angular acceleration of the second link | rad/s² |

**Control input:** `u` – horizontal force on the cart (N).

### 1.2 Parameters (default values)

| Symbol | Value | Meaning |
|--------|-------|---------|
| `M`    | 1.5 kg | cart mass |
| `m₁`   | 0.50 kg | mass of the first link |
| `m₂`   | 0.50 kg | mass of the second link |
| `L₁`   | 0.40 m  | total length of the first link (`L₁ = 2·l₁`) |
| `L₂`   | 0.40 m  | total length of the second link (`L₂ = 2·l₂`) |
| `l₁`   | 0.20 m | centre-of-mass distance of the first link |
| `l₂`   | 0.20 m | centre-of-mass distance of the second link |
| `I₁`   | `¹⁄₁₂·m₁·(0.02 + 4·l₁²)` kg·m² | moment of inertia of the first link |
| `I₂`   | `¹⁄₁₂·m₂·(0.02 + 4·l₂²)` kg·m² | moment of inertia of the second link |
| `g`    | 9.81 m/s² | gravitational acceleration |
| `c`    | 0.00 Nms/rad | joint damping |

For the default values (`m₁ = m₂ = 0.50 kg`, `l₁ = l₂ = 0.20 m`) this gives
`I₁ = I₂ ≈ 0.0075 kg·m²`.

*Note: the moments of inertia are recomputed automatically whenever a parameter
changes. The slider-controlled parameters are `l₁` and `l₂`; the total length
follows as `L = 2·l`.*

### 1.3 Equations of motion (Lagrangian formalism)

The dynamics follow from the Euler-Lagrange equation:

```
d/dt (∂L/∂q̇) − ∂L/∂q = Q
```

with the Lagrangian `L = T − V` (kinetic minus potential energy). The
generalised coordinates are `q = [x, θ₁, θ₂]ᵀ`.

#### Inertia matrix M(q)

```
M₁₁ = M + m₁ + m₂
M₁₂ = (m₁·l₁ + m₂·L₁)·cos(θ₁)
M₁₃ = m₂·l₂·cos(θ₂)
M₂₂ = I₁ + m₁·l₁² + m₂·L₁²
M₂₃ = m₂·L₁·l₂·cos(θ₁−θ₂)
M₃₃ = I₂ + m₂·l₂²
```

Symmetric: `M₂₁ = M₁₂, M₃₁ = M₁₃, M₃₂ = M₂₃`.

#### Coriolis matrix C(q, q̇)

```
C = [[0, −(m₁·l₁+m₂·L₁)·θ̇₁·sin(θ₁), −m₂·l₂·θ̇₂·sin(θ₂)],
     [0, 0, m₂·L₁·l₂·θ̇₂·sin(θ₁−θ₂)],
     [0, −m₂·L₁·l₂·θ̇₁·sin(θ₁−θ₂), 0]]
```

#### Gravity vector G(q)

```
G(q) = [0, −(m₁·l₁+m₂·L₁)·g·sin(θ₁), −m₂·g·l₂·sin(θ₂)]ᵀ
```

#### Equation of motion (implicit form)

```
M(q)·q̈ + C(q,q̇)·q̇ + G(q) = H(u)
```

with `H = [u, 0, 0]ᵀ`. The accelerations are solved for q̈:

```
q̈ = M⁻¹(q)·(H(u) − C(q,q̇)·q̇ − G(q))
```

Integration uses **classical Runge-Kutta 4 (RK4)** with a fixed step size of
`dt = 6.25·10⁻⁵ s` (16 kHz).

#### Energies

Potential energy (relative to the upright rest position). Link 1 carries the
mass of link 2 at its tip, which is why `m₂·L₁` enters:

```
Epot = (m₁·l₁ + m₂·L₁)·g·(cos(θ₁)−1) + m₂·l₂·g·(cos(θ₂)−1)
```

Kinetic energy of the links — including the coupling term `M₂₃`:

```
Ekin = ½·(M₂₂·θ̇₁² + 2·M₂₃·θ̇₁·θ̇₂ + M₃₃·θ̇₂²)
```

---

## 2. Partial Feedback Linearization (PFL)

PFL decouples the desired cart acceleration `v` from the pendulum dynamics.
Instead of choosing a force `u` directly, a virtual cart acceleration `v` is
specified and the force `u` needed to realise it is computed.

### 2.1 Derivation

From the first row of `M(q)·q̈ + C·q̇ + G = H(u)`:

```
M₁₁·ẍ + M₁₂·θ̈₁ + M₁₃·θ̈₂ + C₁ + G₁ = u
```

The link accelerations `θ̈₁, θ̈₂` depend affinely on `v = ẍ`:

```
θ̈₁ = θ̈₁ᵖ + θ̈₁ᵛ·v
θ̈₂ = θ̈₂ᵖ + θ̈₂ᵛ·v
```

where `θ̈ᵖ` is the link acceleration without cart motion (Coriolis + gravity +
damping only) and `θ̈ᵛ` the acceleration per unit of `v`:

```
θ̈₁ᵖ = (M₃₃·(−C₂−G₂−c·θ̇₁) − M₂₃·(−C₃−G₃−c·θ̇₂)) / (M₂₂·M₃₃ − M₂₃²)
θ̈₂ᵖ = (−M₃₂·(−C₂−G₂−c·θ̇₁) + M₂₂·(−C₃−G₃−c·θ̇₂)) / (M₂₂·M₃₃ − M₂₃²)
θ̈₁ᵛ = (−M₃₃·M₂₁ + M₂₃·M₃₁) / (M₂₂·M₃₃ − M₂₃²)
θ̈₂ᵛ = (M₃₂·M₂₁ − M₂₂·M₃₁) / (M₂₂·M₃₃ − M₂₃²)
```

Here `C₁, C₂, C₃` are the components of `C(q,q̇)·q̇` and therefore
**quadratic** in the joint rates:

```
C₁ = −(m₁·l₁+m₂·L₁)·sin(θ₁)·θ̇₁²  −  m₂·l₂·sin(θ₂)·θ̇₂²
C₂ =  m₂·L₁·l₂·sin(θ₁−θ₂)·θ̇₂²
C₃ = −m₂·L₁·l₂·sin(θ₁−θ₂)·θ̇₁²
```

### 2.2 PFL control law

```
u = α·v + φ

α = M₁₁ + M₁₂·θ̈₁ᵛ + M₁₃·θ̈₂ᵛ
φ = M₁₂·θ̈₁ᵖ + M₁₃·θ̈₂ᵖ + C₁ + G₁
```

Any desired value of `v` is thus turned into a force `u` exactly.

### 2.3 Energy balance

With `v` as the input, the energy of the link pair

```
E   = ½·θ̇ᵀ·M_p·θ̇ + a₁·g·cos(θ₁) + a₂·g·cos(θ₂)
a₁  = m₁·l₁ + m₂·L₁,   a₂ = m₂·l₂
M_p = [[M₂₂, M₂₃], [M₂₃, M₃₃]]
```

obeys exactly

```
Ė = −v·W − c·(θ̇₁² + θ̇₂²),   W = a₁·cos(θ₁)·θ̇₁ + a₂·cos(θ₂)·θ̇₂
```

(the Coriolis terms cancel against `½·θ̇ᵀ·Ṁ_p·θ̇`). `W` is the only channel
through which the cart can add or remove energy — this identity carries both
the settling phase (section 6) and the energy-pumping fallback.

---

## 3. State machine

The controller runs in three states:

| `recoveryState` | Mode | Description |
|------------------|------|-------------|
| 0 | `LQR Balance` | state feedback about the upper equilibrium |
| 1 | `Settling` | drain energy down to the hanging rest position, centre the cart |
| 2 | `Swing-up` | execute the optimised swing-up trajectory |

### Transitions

```
Power ON / Reset → 0 if the state is inside the catch region, otherwise → 1
State 1: hanging rest reached ∧ trajectory available → 2
State 2: catch region reached → 0
State 2: trajectory ran out without a catch → 1 (new attempt)
State 0: |θ₁|>0.8 or |θ₂|>0.8 for >20 cycles → 1
Rail stop (±2.5 m): → 1
```

The catch region is **not** heuristic; it is derived from the LQR's Riccati
matrix (section 4.3).

---

## 4. LQR state feedback (RecoveryState = 0)

### 4.1 Linearisation

At `θ₁ = θ₂ = 0, q̇ = 0` the matrix `M(q)` is constant, `C` vanishes and `G` is
linear. With `N = M₀⁻¹`:

```
q̈ = N·([u,0,0]ᵀ − ∂G/∂q·q − D·q̇)
```

and from this the state model `ẋ = A·x + B·u` with
`x = [x, ẋ, θ₁, θ̇₁, θ₂, θ̇₂]ᵀ`. For the default parameters the open-loop plant
has the eigenvalues

```
0, 0, ±4.96, ±11.44   (1/s)
```

— two unstable poles, as expected for a double inverted pendulum.

### 4.2 Solving the Riccati equation

The gain is computed as a genuine LQR solution, not estimated. What is solved is

```
AᵀP + P·A − P·B·R⁻¹·Bᵀ·P + Q = 0
```

via the **matrix sign function** of the Hamiltonian

```
H = [[A, −B·R⁻¹·Bᵀ], [−Q, −Aᵀ]]
```

Newton's iteration `Z ← ½·(c·Z + (c·Z)⁻¹)` converges quadratically to
`sign(H)`; the stable invariant subspace is spanned by the columns of
`I − sign(H)` and has the form `[u; P·u]`, from which `P` follows by least
squares. A few Newton/Kleinman steps (the Lyapunov equation as a 36×36 system)
polish `P` to machine precision. The residual is `≈10⁻¹⁰` and the runtime a few
milliseconds — so the synthesis is rerun on every parameter change.

Weights:

```
Q = diag(5.0 ; 0.1 ; 10 ; 50 ; 10 ; 50)
R = 0.5
```

The light weighting of cart position and the heavy weighting of angular rates
maximise the basin of attraction reachable within the ±50 N actuator limit —
and that basin is exactly what the swing-up has to hit.

For the default parameters this gives

```
K ≈ [3.16 ; 6.49 ; −303.4 ; −8.91 ; 343.9 ; 41.4]
u  = −K·(x − set-point vector),   u = clamp(u, −50, 50)
```

with closed-loop poles `−24.60`, `−8.10 ± 4.31i`, `−1.70`, `−0.71 ± 0.81i`.

Measured basin of attraction (nonlinear simulation, ±50 N):

| Deflection | Limit |
|------------|-------|
| `θ₁ = θ₂`  | ±0.59 rad |
| `θ₁` only  | ±0.36 rad |
| `θ₂` only  | ±0.41 rad |
| `θ̇₁ = θ̇₂`  | ±2.1 rad/s |

### 4.3 Catch region

Control is handed to the LQR as soon as

```
eᵀ·P·e < catchLevel   ∧   |θ₁|,|θ₂| < 0.6   ∧   |θ̇₁|,|θ̇₂| < 4.0   ∧   |x−x_soll| < 2.0
```

The explicit bounds are only a sanity guard (no catch on the far side of an
angle wrap, or right at a rail stop) — the sub-level set is what decides.

Cart position is **deliberately** left out of the quadratic form: at the moment
of hand-over the LQR set point is moved to the current cart position and then
returned to `x_soll` with a time constant of 1.5 s. That way the entire force
budget is available for catching the links instead of for a position error that
is about to disappear anyway.

The threshold is the largest sub-level set of `eᵀPe` on which the unsaturated
LQR force still respects the actuator limit:

```
max |K·e| over eᵀPe ≤ c  =  √(c · K·P⁻¹·Kᵀ)
⟹ catchLevel = U_max² / (K·P⁻¹·Kᵀ)     (≈ 28.5 for the default parameters)
```

Monte-Carlo sampling (~100 accepted states each, then simulated forward for
12 s) confirms the set completely for most plants. Two outliers remained: at
`l₂ = 0.5 m` 1 of 106 accepted states was lost, at `g = 1 m/s²` 4 of 135 — so
there the bound is a few per cent optimistic, because it comes from the
*linearised* model.

Two things absorb this:

* The **planned** swing-up trajectory does not depend on this bound. Before it
  is adopted it is simulated in full — swing-up *and* the subsequent balancing
  — and accepted only if the LQR really brings the pendulum to rest upright
  (section 5.3).
* If the test does misfire after a **disturbance**, the divergence monitor
  detects it within 20 computation cycles; the controller settles and swings up
  again.

---

## 5. Swing-up by trajectory optimisation (RecoveryState = 2)

### 5.1 Why not pure energy pumping

Energy pumping (`v ∝ (E−E_d)·W`) reliably drives the **total energy** to its
value at the upper equilibrium, but says nothing about **how** that energy is
split between the two links. The set `E = E_d` is five-dimensional and the
catch region within it is vanishingly small. Measured on the default plant: in
120 s of energy pumping both angles are simultaneously below 0.5 rad for only
**0.06 s** in total, and the best value of `eᵀPe` reached was 178 against a
threshold of 28.5. Swing-up succeeded in fewer than a third of the attempts,
with median times around 24 s.

### 5.2 iLQR trajectory optimisation

Instead, the swing-up is solved as an optimal control problem — **iLQR** on the
full nonlinear model, starting from the hanging rest position:

```
min  ½·w_f·e_Nᵀ·P·e_N  +  Σ_k [ ½·r·w_k² + ½·q_x·x_k²
                               + ρ(k)·(a₁·(1−cosθ₁) + a₂·(1−cosθ₂))
                               + ½·ρ_v(k)·(θ̇₁² + θ̇₂²)  + rail barrier ]
```

| Element | Value / purpose |
|---------|-----------------|
| knots `N` | 125 (coarse) or up to 220 (fine), depending on the restart |
| horizon `T` | `2.5 s · 4.22 / ω_hang`, clamped to 1.5–9.0 s, plus ×1.5 or ×0.75 depending on the restart |
| force | `u = 50·tanh(w)` — the actuator limit holds by construction, the problem stays unconstrained |
| terminal cost | `eᵀPe` = LQR cost-to-go, i.e. exactly the measure of the catch region |
| `ρ(k), ρ_v(k)` | "get upright" shaping terms with weight `∝ (k/N)³` |
| rail barrier | soft, from `|x| > 1.2 m` |

The horizon scales with the slow eigenfrequency `ω_hang` of the hanging double
pendulum (from `det(K_g − ω²M_p) = 0`), so that long or weakly accelerated
pendulums get proportionally more time.

The terminal cost alone has too small a basin for the optimiser; the smooth
"get upright" terms lead it there. Measured over five random starts each with
identical start noise: without the shaping terms one start reached a terminal
cost below 1, with them three. Over 14 starts with the shaping terms, nine
reached a usable solution.

**Time slicing:** the optimisation runs in slices of ≤ 8 ms per animation
frame, while the controller is still settling — the interface stays smooth. A
restart still far off after 30 iterations (more than 22× the catch region) is
discarded rather than nursed along; the last six attempts of a round always run
to completion. Restarts cycle through four combinations of horizon and
resolution — 125 or up to 220 knots, horizon ×1, ×1.5 or ×0.75. No choice wins
everywhere: the coarse grid converges more often on the default parameters,
very long links only produce usable trajectories on the fine one, and at low
`g` it only succeeds with a stretched horizon. Since every candidate is
simulated through anyway (section 5.3), extra diversity here can only help.

### 5.3 Acceptance by closed-loop simulation

A swing-up trajectory is **open-loop unstable**. The terminal cost of the
optimisation alone therefore says little: the same forces re-integrated on a
different time grid run away — one candidate scoring `eᵀPe = 0.1` on its own
grid was measured to end at `eᵀPe > 12 000` when replayed open loop.

What matters is the **controlled** behaviour. Every candidate is therefore
simulated through completely before adoption — step size 0.5 ms, actuator limit
and rail clearance active:

1. swing-up along the trajectory with the time-varying LQR, until the catch
   test fires;
2. then four times the planning horizon with the balancing LQR, including the
   return of the cart set point. The window length has to follow the plant: at
   `g = 1 m/s²` the closed loop is several times slower than at 9.81, and a
   fixed window would discard usable trajectories merely for not having settled
   yet.

Only what passes both is adopted — the pendulum must end upright and at rest.
The delivered trajectory therefore no longer depends on whether the catch
threshold is exactly right in a marginal case. The displayed `eᵀPe` is the
value reached at the catch. A rejected candidate is discarded and the search
continues.

If a swing-up attempt fails anyway (trajectory ran out, or a rail stop was
hit), another search round of 20 restarts starts in the background while the
controller settles again — up to ten rounds.

### 5.4 Time-varying LQR for trajectory tracking

The planned trajectory is not run open loop but **tracked**. Backwards along
the trajectory the discrete Riccati recursion is solved:

```
S_k = Q_t + AₖᵀS_{k+1}Aₖ − AₖᵀS_{k+1}Bₖ·(R_t + BₖᵀS_{k+1}Bₖ)⁻¹·BₖᵀS_{k+1}Aₖ
K_k = (R_t + BₖᵀS_{k+1}Bₖ)⁻¹·BₖᵀS_{k+1}Aₖ
```

with `S_N = P` (the LQR matrix), `Q_t = diag(1,1,20,2,20,2)`, `R_t = 0.05`.
`Aₖ, Bₖ` follow from finite differences of one RK4 step. The applied force is

```
u = u_nom[k] − K_k·(x − x_nom[k]),   u = clamp(u, −50, 50)
```

This feedback is what makes the method robust: the plan is computed with step
sizes around 20 ms while the simulation runs at 62.5 µs. Measured tolerance to
deviations in the start state: ±0.4 m, ±0.4 m/s, ±0.4 rad and ±0.4 rad/s are
still caught (15/15 attempts).

---

## 6. Settling (RecoveryState = 1)

The trajectory is planned from the hanging rest position. The controller gets
there using the energy balance from section 2.3: `v = +k_w·W` gives
`Ė = −k_w·W² ≤ 0` and therefore removes energy unconditionally.

That alone stalls, however. The hanging double pendulum has two eigenmodes —
for the default parameters at 4.22 rad/s with shape `[1; 1.448]` (in phase) and
10.94 rad/s with `[1; −2.073]` (out of phase). Their input coupling `M_pxᵀ·u`
is −0.4448 and −0.0927 respectively: the out-of-phase mode is roughly **five
times more weakly** coupled to the cart. It therefore barely decays while
`W ≈ 0`.

A second term damps exactly that mode directly, by inverting the coupling
`∂(θ̈₁−θ̈₂)/∂v`:

```
v = k_w·W
    − clamp( (k_rel·(θ̇₁−θ̇₂)) / (θ̈₁ᵛ − θ̈₂ᵛ), ±8 )
    − k_x·(x − x_soll) − k_d·ẋ
v = clamp(v, −12, 12)

k_w = 6.0   k_rel = 3.0   k_x = 1.5   k_d = 1.5
```

With this the system reaches the rest position from every tested initial state
in 1–4 s, without touching the rail stops.

---

## 7. Rail limits

The rail has hard stops at `x = ±2.5 m`. On contact the cart is set to the
limit (`ẋ = 0`) and the controller switches to settling mode (state 1). The
test runs on every integration step, not once per frame, so the cart cannot
overshoot the limit.

The planned trajectory keeps its distance by itself: the soft rail barrier in
the optimisation already acts from `|x| > 1.2 m`, and in every case tested the
trajectory stayed below 1.25 m.

In addition the force is reduced proportionally near the limits (±2.0 m to
±2.5 m) when the cart is moving towards the stop:

```
margin = 0.5 m
scale = max(0, (2.5 − |x|) / margin)
u *= scale  (if ẋ points towards the limit)
```

---

## 8. Simulation parameters

| Parameter | Control | Range | Default |
|-----------|---------|-------|---------|
| cart mass `M` | slider | 0.5–5.0 kg | 1.50 kg |
| link 1 mass `m₁` | slider | 0.05–2.00 kg | 0.50 kg |
| link 2 mass `m₂` | slider | 0.05–3.50 kg | 0.50 kg |
| centre of mass 1 `l₁` | slider | 0.10–1.00 m | 0.20 m |
| centre of mass 2 `l₂` | slider | 0.10–1.00 m | 0.20 m |
| total length `L₁` | automatic | — | 0.40 m (`L₁ = 2·l₁`) |
| total length `L₂` | automatic | — | 0.40 m (`L₂ = 2·l₂`) |
| joint damping `c` | slider | 0.00–0.50 Nms/rad | 0.00 Nms/rad |
| gravity `g` | slider | 1.0–35.0 m/s² | 9.81 m/s² |
| adaptation | toggle | ON/OFF | ON |
| rail limits | toggle | ON/OFF | ON |

*Note: with "Adaptation: ON" the LQR gains are recomputed on every parameter
change. With "OFF" the gain stays frozen at the value it had when adaptation
was switched off — which lets you watch how the controller reacts to a
mis-modelled plant. The Riccati matrix `P` (and with it the catch test and the
planner's terminal cost) follows the actual plant in both cases.*

*Every parameter change also discards the swing-up trajectory and restarts
planning — which continues in the background while the controller balances or
settles.*

### Verified operating range

Every run starts from the hanging rest position and counts as passed if the
pendulum is subsequently held in the up-up position (`|θ| < 0.03 rad`, cart at
`x_soll`):

| Parameter | tested | Result |
|-----------|--------|--------|
| `M` | 0.5 – 5.0 kg | complete |
| `m₁` | 0.05 – 2.0 kg | complete |
| `m₂` | 0.05 – 3.5 kg | complete |
| `l₁` | 0.10 – 1.00 m | complete |
| `l₂` | 0.10 – 0.50 m | complete |
| `g` | 1.0 – 20.0 m/s² | complete |
| `c` | 0.00 – 0.20 Nms/rad | complete |

Also tested: all three scenarios, the four target-angle buttons, disturbance
impulses (up to `Δẋ = −2 m/s`, `Δθ̇₁ = Δθ̇₂ = −3 rad/s` from the balanced
position) and moving `x_soll` while balancing.

**Known limits.** At three ends of the sliders the optimisation finds no
trajectory within its attempt budget that passes acceptance:

* `l₂ = 1.0 m` — a 2 m long upper link,
* `g = 35 m/s²` — the largest selectable gravity,
* `c = 0.50 Nms/rad` — the largest joint damping.

Behaviour stays well mannered there: the controller settles the pendulum, keeps
planning in the background and tries again; the display reads
"No plan — energy pumping". **Balancing** itself — scenario A and the
disturbance tests — works at those settings too.

---

## 9. Scenarios

| Scenario | Initial condition | Expectation |
|----------|-------------------|-------------|
| **A** – stabilisation | `x=0, θ₁=0.05, θ₂=−0.03` | LQR holds directly; upright after ≈3 s |
| **B** – free fall | `x=0, θ₁=1.57, θ₂=0` | pendulum falls, is settled, swings up; upright after ≈8 s |
| **C** – swing-up | `x=0, θ₁=π, θ₂=π` | swing-up from rest; upright after ≈3 s (plan ready) |

In all three cases the up-up position is held to better than `10⁻³ rad` and the
cart returns to `x_soll`.

---

## 10. Interaction

| Action | Description |
|--------|-------------|
| **Power ON/OFF** | switch the controller on/off |
| **Reset State** | reset the state to the scenario's initial value |
| **Shock Test** | random impulse on cart and pendulum |
| **Download CSV** | export data as CSV |
| **Adaptation ON/OFF** | enable/disable adaptive LQR gains |
| **Limits ON/OFF** | enable/disable the rail limits |
| **Mouse wheel** | zoom the animation (10 % – 300 %) |
| **Click on the rail** | set the position set point `x_soll` |
| **Target-angle buttons** | set `(0,0)`, `(0,π)`, `(π,0)`, `(π,π)` |

---

## 11. Data output (CSV)

Clicking "Download CSV" exports up to 10 000 rows of:

```
time, x, dx, theta1, dtheta1, theta2, dtheta2, force_N, Epot_J, Ekin_J
```

Time is the frame index. All angles are wrapped to `(−π, π]`. The energies are
relative to the upright rest position.

---

## 12. Implementation details

- **Simulation step:** `dt = 62.5 µs` (16 kHz), RK4
- **Actuator limit:** `|u| ≤ 50 N`, rail `|x| ≤ 2.5 m`
- **LQR:** exact CARE solution via the matrix sign function of the Hamiltonian,
  polished with Newton/Kleinman; residual `≈10⁻¹⁰`, runtime a few milliseconds
- **Catch test:** `eᵀPe < U_max²/(K·P⁻¹·Kᵀ)` — the largest sub-level set without
  actuator saturation
- **Swing-up:** iLQR, 125 knots, horizon tied to the slow eigenfrequency of the
  hanging pendulum, force bounded via `tanh`
- **Trajectory tracking:** time-varying LQR from the discrete Riccati recursion
  along the trajectory, terminal weight `S_N = P`
- **Background planning:** ≤ 8 ms per animation frame, at least one optimisation
  step per frame; replanning on every parameter change
- **PFL:** turns a desired cart acceleration into a force exactly; Coriolis
  terms quadratic in the joint rates
- **State machine:** 3 states (balance, settling, swing-up)
- **Display:** 6 individual plots with a 150-point sliding window
- **Language:** plain JavaScript, no libraries

### Numerics in detail

| Building block | Method | Size |
|----------------|--------|------|
| matrix inversion | Gauss-Jordan with column pivoting | up to 12×12 |
| `sign(H)` | Newton iteration with norm scaling | 12×12 |
| Lyapunov equation | direct linear system | 36×36 |
| iLQR Jacobians | central differences of one RK4 step | 6×6, 6×1 |
| line search | backtracking over 9 step sizes | — |
