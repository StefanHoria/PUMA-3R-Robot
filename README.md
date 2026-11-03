***English** · [Română](README.ro.md)*

# PUMA-3R-Robot

A **PUMA-type positioning structure with 3 revolute joints**: modelling, simulation and PI control in MATLAB / Simulink.

Semester project – Robotics, 4th year, Transilvania University of Brașov
Supervisor: Prof. Dr.-Ing. Claudiu Pozna · Author: Ștefan Horia-Eusebiu (4LF832)

---

## What it contains

- **Geometric model**, forward and inverse (homogeneous matrices 0_1T · 1_2T · 2_3T · 3_eT)
- **Kinematic model**: geometric Jacobian, velocities and accelerations, tracking of a trajectory **r(t)**
- **Dynamic model** (Lagrange): M(q), C(q, q̇), G(q) – inverse and forward dynamics, validated through
  energy conservation and the skew-symmetry of dM/dt − 2C
- **Actuators**: maxon RE 40 DC motor + GP 52 C / GP 42 C planetary gearheads
- **PI control** on each joint, tuned by pole placement
- **Simulink model** `contrPUMA.slx`, generated automatically, with the structure Controller – Motor/Gearhead –
  inverse dynamic model of the structure
- **3D visualisation** (3D view + XY and YZ projections) with animations: one tab per stage, Play / Pause
  buttons, a time slider and playback speed

## Files

| File | Role |
|---|---|
| `Date_PUMA.m` | project data: dimensions, material, motions, motors, settings – **the only file you edit** |
| `PUMA_simulare.m` | all the code: models, actuators, controllers, building the Simulink model, plots, animations |
| `Proiect_PUMA.mlx` | the project report (live script) |
| `contrPUMA.slx` | Simulink model of the control system (regenerated on every run) |
| `Poz_dorita.mat` | desired joint trajectory, read by the From File block (generated) |
| `Brat3Cuple.mlx` | reference example – another 3-joint arm (symbolic Lagrange / Newton-Euler derivation) |

## Requirements

- MATLAB **R2026b** + **Simulink** (Control System Toolbox – optional)
- `contrPUMA.slx` is saved in R2026b. On an older release, run `Date_PUMA.m`, which rebuilds the model.

## Running it

1. Open the project folder in MATLAB.
2. Run **`Date_PUMA.m`** (F5).
3. The following open:
   - the **„Braț Robotic – Vizualizare Complexă”** window – one tab per stage, with the animations
     (Play, Pause, time slider, speed);
   - the **„PUMA – Grafice”** window – tabs per stage and sub-tabs with the plots;
   - the Simulink model **`contrPUMA.slx`**.

The Simulink model can also be run directly (Run), with the project folder as the current folder; the results
end up in the `out` variable: `out.miscari` (q, q̇, q̈), `out.moMo` (joint torques), `out.erori`,
`out.comanda` (voltages), `out.turatii` (speeds).

## Stages (from `Date_PUMA.m`)

| Stage | Data | What it does |
|---|---|---|
| 1. Static pose | `q_static` | drawing of the robot, 0_eT, gravity torques |
| 2. Joint space | `q_puncte` | polynomial motion (3rd / 5th order) through several configurations |
| 3. Tracking r(t) | `r_puncte` | Cartesian trajectory → inverse kinematic model → torques |
| 4. Free motion | `miscare_libera` | falling under gravity (τ = 0), energy conservation |
| 5. PI control | `motor`, `control` | actuators, PI controllers, the Simulink model, checks |

A stage is disabled by leaving its data empty (`[]`) or setting its flag to `0`.

## Actuators

| Joint | Motor | Gearhead | Available torque | Required torque |
|---|---|---|---|---|
| 1 – base | RE 40, 24 V (148867) | GP 52 C 230:1 (223101) | 27.7 N·m | 4.1 N·m |
| 2 – shoulder | RE 40, 24 V (148867) | GP 52 C 230:1 (223101) | 27.7 N·m | 21.8 N·m (static) |
| 3 – elbow | RE 40, 24 V (148867) | GP 42 C 74:1 (203123)  | 9.4 N·m  | 4.8 N·m (static) |

Catalogue data: [RE 40](https://www.maxongroup.com/maxon/view/product/148867) ·
[GP 52 C 230:1](https://www.maxongroup.com/maxon/view/product/223101) ·
[GP 42 C 74:1](https://www.maxongroup.com/maxon/view/product/203123)

## Control

Model of the motor with gearhead (for each joint):

```
q̈ = k2·v − k1·q̇ − k3·τ,   k1 = (bm + km²/R)/Ja,   k2 = km/(Ja·ir·R),   k3 = 1/(Ja·η·ir²)
```

The PI controller is tuned by pole placement on G(s) = k2 / (s(s + k1)):

```
P(s) = s³ + k1·s² + k2·Kp·s + k2·Ki = (s + p1)(s + p2)(s + p3)
p3 = k1 − p1 − p2,   Kp = (p1p2 + p1p3 + p2p3)/k2,   Ki = p1p2p3/k2
```

- k1 and k2 use the equivalent motor + arm inertia, J_ech = Ja + M̄_ii/(η·ir²), because the PUMA arm as
  seen by the motor is not negligible (k3·M_ii up to 4.2).
- The motor ↔ structure algebraic loop is broken with a first-order filter on τ (Tf = 1 ms).
- The Simulink simulation is checked against the same loop simulated directly in MATLAB (difference ~1e-5 rad).

## Results (tracking r(t), 2 × 3 s)

| Joint | Kp [V/rad] | Ki [V/(rad·s)] | max error [rad] | max voltage [V] | RMS current [A] | max speed [rpm] |
|---|---|---|---|---|---|---|
| 1 | 108.8 | 567.5 | 1.2e-2 | 6.9 | 0.22 | 2185 |
| 2 | 216.6 | 2248 | 5.8e-3 | 8.0 | 2.77 | 2763 |
| 3 | 51.0 | 387.9 | 1.4e-2 | 2.7 | 2.17 | 681 |

Maximum end-effector error relative to r(t): **8.8 mm**; in steady state the error is zero.
