# PUMA-3R-Robot

Structura de poziționare **PUMA cu 3 cuple de rotație** – modelare, simulare și control PI în MATLAB / Simulink.

Proiect de semestru – Robotică, anul IV, Universitatea Transilvania din Brașov
Coordonator: Prof. Dr.-Ing. Claudiu Pozna · Autor: Ștefan Horia-Eusebiu (4LF832)

---

## Ce conține

- **Model geometric** direct și invers (matrici omogene 0_1T · 1_2T · 2_3T · 3_eT)
- **Model cinematic**: Jacobian geometric, viteze și accelerații, urmărirea traiectoriei **r(t)**
- **Model dinamic** (Lagrange): M(q), C(q, q̇), G(q) – dinamică inversă și directă, validare prin
  conservarea energiei și antisimetria dM/dt − 2C
- **Acționări**: motor de curent continuu maxon RE 40 + reductoare planetare GP 52 C / GP 42 C
- **Control PI** pe fiecare cuplă, acordat prin alocarea polilor
- **Model Simulink** `contrPUMA.slx`, generat automat, cu structura Controler – MotoReductor – Structura
  model dinamic invers
- **Vizualizare 3D** (vedere 3D + proiecțiile XY și YZ) cu animații: câte un tab pe etapă, butoane
  Redă / Pauză, bară de timp și viteză de redare

## Fișiere

| Fișier | Rol |
|---|---|
| `Date_PUMA.m` | datele proiectului: dimensiuni, material, mișcări, motoare, setări – **singurul fișier care se editează** |
| `PUMA_simulare.m` | tot codul: modele, acționări, regulatoare, construirea modelului Simulink, grafice, animații |
| `Proiect_PUMA.mlx` | memoriul proiectului (live script) |
| `contrPUMA.slx` | modelul Simulink al sistemului de reglare (se regenerează la fiecare rulare) |
| `Poz_dorita.mat` | traiectoria dorită a cuplelor, citită de blocul From File (generată) |
| `Brat3Cuple.mlx` | exemplu de referință – alt braț cu 3 cuple (calcul simbolic Lagrange / Newton-Euler) |

## Cerințe

- MATLAB **R2026b** + **Simulink** (Control System Toolbox – opțional)
- `contrPUMA.slx` este salvat în R2026b. Pe o versiune mai veche se rulează `Date_PUMA.m`, care
  reconstruiește modelul.

## Rulare

1. Se deschide folderul proiectului în MATLAB.
2. Se rulează **`Date_PUMA.m`** (F5).
3. Se deschid:
   - fereastra **„Braț Robotic – Vizualizare Complexă”** – câte un tab pe etapă, cu animațiile
     (Redă, Pauză, bară de timp, viteză);
   - fereastra **„PUMA – Grafice”** – tab-uri pe etape și sub-tab-uri cu graficele;
   - modelul Simulink **`contrPUMA.slx`**.

Modelul Simulink se poate rula și direct (Run), cu folderul proiectului ca folder curent; rezultatele
apar în variabila `out`: `out.miscari` (q, q̇, q̈), `out.moMo` (momentele din cuple), `out.erori`,
`out.comanda` (tensiunile), `out.turatii`.

## Etapele (din `Date_PUMA.m`)

| Etapa | Date | Ce face |
|---|---|---|
| 1. Poziție statică | `q_static` | desenul robotului, 0_eT, cuplurile gravitaționale |
| 2. Spațiul cuplelor | `q_puncte` | mișcare polinomială (grad 3 / 5) prin mai multe configurații |
| 3. Urmărire r(t) | `r_puncte` | traiectorie carteziană → model cinematic invers → cupluri |
| 4. Mișcare liberă | `miscare_libera` | căderea sub greutate (τ = 0), conservarea energiei |
| 5. Control PI | `motor`, `control` | acționări, regulatoare PI, modelul Simulink, verificări |

O etapă se dezactivează lăsând datele ei goale (`[]`) sau punând setarea pe `0`.

## Acționări

| Cupla | Motor | Reductor | Cuplu disponibil | Cuplu necesar |
|---|---|---|---|---|
| 1 – baza | RE 40, 24 V (148867) | GP 52 C 230:1 (223101) | 27.7 N·m | 4.1 N·m |
| 2 – umăr | RE 40, 24 V (148867) | GP 52 C 230:1 (223101) | 27.7 N·m | 21.8 N·m (static) |
| 3 – cot  | RE 40, 24 V (148867) | GP 42 C 74:1 (203123)  | 9.4 N·m  | 4.8 N·m (static) |

Datele de catalog: [RE 40](https://www.maxongroup.com/maxon/view/product/148867) ·
[GP 52 C 230:1](https://www.maxongroup.com/maxon/view/product/223101) ·
[GP 42 C 74:1](https://www.maxongroup.com/maxon/view/product/203123)

## Control

Modelul motorului cu reductor (pentru fiecare cuplă):

```
q̈ = k2·v − k1·q̇ − k3·τ,   k1 = (bm + km²/R)/Ja,   k2 = km/(Ja·ir·R),   k3 = 1/(Ja·η·ir²)
```

Regulatorul PI se acordă prin alocarea polilor pe G(s) = k2 / (s(s + k1)):

```
P(s) = s³ + k1·s² + k2·Kp·s + k2·Ki = (s + p1)(s + p2)(s + p3)
p3 = k1 − p1 − p2,   Kp = (p1p2 + p1p3 + p2p3)/k2,   Ki = p1p2p3/k2
```

- În k1 și k2 se folosește inerția echivalentă motor + braț, J_ech = Ja + M̄_ii/(η·ir²), deoarece
  brațul PUMA văzut de motor nu este neglijabil (k3·M_ii până la 4.2).
- Bucla algebrică motor ↔ structură este ruptă cu un filtru de ordinul 1 pe τ (Tf = 1 ms).
- Simularea Simulink este verificată cu aceeași buclă simulată direct în MATLAB (diferență ~1e-5 rad).

## Rezultate (urmărirea r(t), 2 × 3 s)

| Cupla | Kp [V/rad] | Ki [V/(rad·s)] | eroare max [rad] | tensiune max [V] | I efectiv [A] | turație max [rpm] |
|---|---|---|---|---|---|---|
| 1 | 108.8 | 567.5 | 1.2e-2 | 6.9 | 0.22 | 2185 |
| 2 | 216.6 | 2248 | 5.8e-3 | 8.0 | 2.77 | 2763 |
| 3 | 51.0 | 387.9 | 1.4e-2 | 2.7 | 2.17 | 681 |

Eroarea maximă a efectorului față de r(t): **8.8 mm**; în regim staționar eroarea este nulă.
