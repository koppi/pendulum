# Doppel-Invertpendulum — Regelung im Browser

Diese Simulation zeigt die Regelung eines doppelten invertierten Pendels auf einem
fahrbaren Wagen mittels **LQR-Zustandsrückführung**, **Partial-Feedback-Linearisierungs
(PFL) Einschwingregelung (Swing-Up)** und **PFL-Stabilisierung**.
Die gesamte Berechnung läuft in JavaScript im Browser – in Echtzeit.

---

## 1. Physikalisches Modell

### 1.1 Koordinaten und Zustandsraum

Der Zustandsvektor hat sechs Komponenten:

```
x = [x, ẋ, θ₁, θ̇₁, θ₂, θ̇₂]ᵀ
```

| Größe | Bedeutung | Einheit |
|-------|-----------|---------|
| `x`   | Wagenposition | m |
| `ẋ`   | Wagengeschwindigkeit | m/s |
| `θ₁`  | Winkel des ersten Pendels (0 = aufrecht) | rad |
| `θ̇₁`  | Winkelgeschwindigkeit des ersten Pendels | rad/s |
| `θ₂`  | Winkel des zweiten Pendels (0 = aufrecht) | rad |
| `θ̇₂`  | Winkelgeschwindigkeit des zweiten Pendels | rad/s |
| `θ̈₁`  | Winkelbeschleunigung des ersten Pendels | rad/s² |
| `θ̈₂`  | Winkelbeschleunigung des zweiten Pendels | rad/s² |

**Steuereingang:** `u` – horizontale Kraft auf den Wagen (N).

### 1.2 Parameter (Standardwerte)

| Symbol | Wert | Bedeutung |
|--------|------|-----------|
| `M`    | 1,5 kg | Masse des Wagens |
| `m₁`   | 0,50 kg | Masse des ersten Pendels |
| `m₂`   | 0,50 kg | Masse des zweiten Pendels |
| `L₁`   | 0,40 m  | Gesamtlänge des ersten Pendels (`L₁ = 2·l₁`) |
| `L₂`   | 0,40 m  | Gesamtlänge des zweiten Pendels (`L₂ = 2·l₂`) |
| `l₁`   | 0,20 m | Schwerpunktabstand des ersten Pendels |
| `l₂`   | 0,20 m  | Schwerpunktabstand des zweiten Pendels |
| `I₁`   | `¹⁄₁₂·m₁·(0,02 + 4·l₁²)` kg·m² | Trägheitsmoment des ersten Pendels |
| `I₂`   | `¹⁄₁₂·m₂·(0,02 + 4·l₂²)` kg·m² | Trägheitsmoment des zweiten Pendels |
| `g`    | 9,81 m/s² | Erdbeschleunigung |
| `c`    | 0,00 Nms/rad | Gelenkdämpfung |

Für die Standardwerte (`m₁ = m₂ = 0,50 kg`, `l₁ = l₂ = 0,20 m`) ergeben sich:
`I₁ = I₂ ≈ 0,0075 kg·m²`.

*Hinweis: Die Trägheitsmomente werden bei Parameteränderungen automatisch
neu berechnet. Die Slider-steuerbaren Parameter sind `l₁` und `l₂`, wobei
sich die Gesamtlänge als `L = 2·l` ergibt.*

### 1.3 Bewegungsgleichungen (Lagrange-Formalismus)

Die Dynamik folgt aus der Euler-Lagrange-Gleichung zweiter Art:

```
d/dt (∂L/∂q̇) − ∂L/∂q = Q
```

mit der Lagrangefunktion `L = T − V` (kinetische minus potenzielle Energie).
Die generalisierten Koordinaten sind `q = [x, θ₁, θ₂]ᵀ`.

#### Trägheitsmatrix M(q)

```
M₁₁ = M + m₁ + m₂
M₁₂ = (m₁·l₁ + m₂·L₁)·cos(θ₁)
M₁₃ = m₂·l₂·cos(θ₂)
M₂₂ = I₁ + m₁·l₁² + m₂·L₁²
M₂₃ = m₂·L₁·l₂·cos(θ₁−θ₂)
M₃₃ = I₂ + m₂·l₂²
```

Symmetrisch: `M₂₁ = M₁₂, M₃₁ = M₁₃, M₃₂ = M₂₃`.

#### Coriolis-Matrix C(q, q̇)

```
C = [[0, −(m₁·l₁+m₂·L₁)·θ̇₁·sin(θ₁), −m₂·l₂·θ̇₂·sin(θ₂)],
     [0, 0, m₂·L₁·l₂·θ̇₂·sin(θ₁−θ₂)],
     [0, −m₂·L₁·l₂·θ̇₁·sin(θ₁−θ₂), 0]]
```

#### Gravitationsvektor G(q)

```
G(q) = [0, −(m₁·l₁+m₂·L₁)·g·sin(θ₁), −m₂·g·l₂·sin(θ₂)]ᵀ
```

#### Aufruf der Bewegung (implizite Form)

```
M(q)·q̈ + C(q,q̇)·q̇ + G(q) = H(u)
```

mit `H = [u, 0, 0]ᵀ`. Die Beschleunigungen werden nach q̈ aufgelöst:

```
q̈ = M⁻¹(q)·(H(u) − C(q,q̇)·q̇ − G(q))
```

Die Integration erfolgt per **Klassischem Runge-Kutta 4 (RK4)** mit
einer festen Schrittweite von `dt = 6,25·10⁻⁵ s` (16 kHz).

#### Energien

Potenzielle Energie (relativ zur aufrechten Ruhelage):

```
Epot = m₁·g·l₁·(cos(θ₁)−1) + m₂·g·l₂·(cos(θ₂)−1)
```

Kinetische Energie:

```
Ekin = ½·(I₁+m₁·l₁²)·θ̇₁² + ½·(I₂+m₂·l₂²)·θ̇₂²
```

---

## 2. Partial Feedback Linearization (PFL)

Die PFL entkoppelt die gewünschte Wagenbeschleunigung `v` von der
Pendeldynamik. Statt direkt eine Kraft `u` zu wählen, wird eine
virtuelle Wagenbeschleunigung `v` spezifiziert und die dazu nötige
Kraft `u` berechnet.

### 2.1 Herleitung

Aus der ersten Zeile von `M(q)·q̈ + C·q̇ + G = H(u)` folgt:

```
M₁₁·ẍ + M₁₂·θ̈₁ + M₁₃·θ̈₂ + C₁ + G₁ = u
```

Die Pendelbeschleunigungen `θ̈₁, θ̈₂` hängen von `v = ẍ` ab:

```
θ̈₁ = θ̈₁ᵖ + θ̈₁ᵛ·v
θ̈₂ = θ̈₂ᵖ + θ̈₂ᵛ·v
```

wobei `θ̈ᵖ` die Pendelbeschleunigung ohne Wagenbewegung ist (nur
Coriolis + Gravitation + Dämpfung) und `θ̈ᵛ` die Beschleunigung pro
Einheit `v`:

```
θ̈₁ᵖ = (M₃₃·(−C₂−G₂−c·θ̇₁) − M₂₃·(−C₃−G₃−c·θ̇₂)) / (M₂₂·M₃₃ − M₂₃²)
θ̈₂ᵖ = (−M₃₂·(−C₂−G₂−c·θ̇₁) + M₂₂·(−C₃−G₃−c·θ̇₂)) / (M₂₂·M₃₃ − M₂₃²)
θ̈₁ᵛ = (−M₃₃·M₂₁ + M₂₃·M₃₁) / (M₂₂·M₃₃ − M₂₃²)
θ̈₂ᵛ = (M₃₂·M₂₁ − M₂₂·M₃₁) / (M₂₂·M₃₃ − M₂₃²)
```

### 2.2 PFL-Stellgesetz

```
u = α·v + φ

α = M₁₁ + M₁₂·θ̈₁ᵛ + M₁₃·θ̈₂ᵛ
φ = M₁₂·θ̈₁ᵖ + M₁₃·θ̈₂ᵖ + C₁ + G₁
```

Damit kann jeder gewünschte `v`-Wert direkt in eine Kraft `u`
umgesetzt werden, die diese Wagenbeschleunigung exakt realisiert.

---

## 3. Zustandsautomat

Der Regler arbeitet in drei Zuständen:

| `recoveryState` | Modus | Beschreibung |
|------------------|-------|--------------|
| 0 | LQR | Zustandsrückführung nahe der aufrechten Lage |
| 1 | PFL-Stabilisierung | Wagenzentrierung und Pendeldämpfung |
| 2 | PFL-Swing-Up | Energiepumpen zum Einschwingen |

### Zustandsübergänge

```
Reset/Power ON → 0 (falls |θ₁|<0,35 und |θ₂|<0,35) oder 2 (sonst)
Zustand 0: |θ₁|>0,8 oder |θ₂|>0,8 für >20 Zyklen → Zustand 2
Zustand 1 (nach 5 s Timeout):
    falls |θ₁|<0,35 und |θ₂|<0,35 → Zustand 0
    sonst → Zustand 2 (30 s Timeout)
Zustand 2 (Einfang): → Zustand 1 (5 s Stabilisierung)
Schienen-Anschlag (±2,5 m): → Zustand 1 (10 s Timeout)
```

---

## 4. LQR-Zustandsrückführung (RecoveryState = 0)

Nahe der aufrechten Gleichgewichtslage regelt ein **LQR-ähnlicher
Zustandsregler** das System. Die Verstärkungsmatrix wird adaptiv aus
der Massenmatrix berechnet.

### 4.1 Regelgesetz

```
u = −(K₀·(x − x_soll) + K₁·ẋ + K₂·θ₁ + K₃·θ̇₁ + K₄·θ₂ + K₅·θ̇₂)
u = clamp(u, −50, 50)
```

Wirksamkeitsbereich: `|θ₁| < 0,8 rad` und `|θ₂| < 0,8 rad`.

Divergenzerkennung: Wenn `|θ₁| > 0,8` oder `|θ₂| > 0,8` für mehr als
20 aufeinanderfolgende Zyklen, wird auf Swing-Up (Zustand 2) geschaltet.

### 4.2 Adaptive Verstärkungssynthese

Die Verstärkungen werden aus der Massenmatrix des aktuellen
Parameterzustands berechnet:

```
scale = M₁₁ / 2,25

K₀ =  5·scale    (Position)
K₁ = 15·scale    (Geschwindigkeit)
K₂ = −120·scale  (Winkel 1)
K₃ = −30·scale   (Winkelgeschw. 1)
K₄ = −80·scale   (Winkel 2)
K₅ = −20·scale   (Winkelgeschw. 2)
```

Jede Verstärkung ist auf `[-800, 800]` begrenzt. **Adaptation:**
Bei aktivierter Adaptation werden die Verstärkungen bei
Parameteränderungen automatisch neu berechnet. Bei deaktivierter
Adaptation bleiben die letzten Verstärkungen erhalten.

Für die Standardparameter (`M = 1,5 kg`, `m₁ = m₂ = 0,50 kg`)
ergibt sich `scale ≈ 1,11` und damit näherungsweise:

```
K ≈ [5,6, 16,7, −133,3, −33,3, −88,9, −22,2]
```

---

## 5. PFL-Einschwingregelung (Swing-Up, RecoveryState = 2)

Statt einer direkten Kraft `u` wird eine virtuelle Beschleunigung `v`
aus vier Komponenten zusammengesetzt:

```
v = v_pump + v_sin + v_cart + v_anti
```

#### Energiepumpe (v_pump)

Pumpt Energie in das System, solange die gewichtete Gesamtenergie
unter null liegt:

```
E₁ = ½·(I₁+m₁·l₁²)·θ̇₁² + m₁·g·l₁·(cos(θ₁)−1)
E₂ = ½·(I₂+m₂·l₂²)·θ̇₂² + m₂·g·l₂·(cos(θ₂)−1)
E_w = max(2·E₁ + 0,5·E₂, −10)
v_pump = pg · min(E_w, 0) · θ̇₁ · cos(θ₁)
  pg = min(15, 3 + t·0,3)
```

Die Gewichtung bevorzugt Pendel 1 (schwerer, langsamer) und bremst
Pendel 2.

#### Sinusförmiger parametrischer Antrieb (v_sin)

Treibt das System bei der Eigenfrequenz des ersten Pendels an:

```
ω = √(g / (I₁/(m₁·l₁) + l₁))
v_sin = min(10, 3 + t·0,2) · sin(ω·t)
```

#### Wagenzentrierung (v_cart)

Zieht den Wagen zur Zielposition `x_soll` zurück:

```
v_cart = −cg · (x − x_soll) − cg · 0,3 · ẋ
  cg = max(1, 3 − t·0,03)
```

Die Zentrierungsverstärkung sinkt mit der Zeit, um dem Pump-Term
mehr Freiheit zu geben.

#### Anti-Spin (v_anti)

Begrenzt übermäßige Drehzahlen direkt:

```
v_anti = 0
wenn |θ̇₂| > 8: v_anti += −1,5·(θ̇₂ − sign(θ̇₂)·8)
wenn |θ̇₁| > 7: v_anti += −1,0·(θ̇₁ − sign(θ̇₁)·7)
```

### PFL-Kraftumsetzung

```
u = pflForce(v_desired, x, ẋ, θ₁, θ̇₁, θ₂, θ̇₂)
v_desired = clamp(v_pump + v_sin + v_cart + v_anti, −15, 15)
u = clamp(u, −35, 35)
```

### Einfangkriterium

Der Swing-Up schaltet auf PFL-Stabilisierung, wenn alle Bedingungen
erfüllt sind:

```
|θ₁| < 0,55 rad  ∧  |θ₂| < 0,75 rad  ∧  E₁+E₂ > −1,5 J
|θ̇₁| < 5,0 rad/s  ∧  |θ̇₂| < 8,0 rad/s
|x| < 2,0 m  ∧  |ẋ| < 2,0 m/s
```

### Timer

- Swing-Up-Timeout: 30 s (danach neuer Versuch)
- Stabilisierungs-Timeout: 5 s (nach Swing-Up-Einfang) bzw. 10 s (nach Schienen-Anschlag)

---

## 6. PFL-Stabilisierung (RecoveryState = 1)

Wird der Schienenbegrenzer (±2,5 m) erreicht oder das System nach
dem Swing-Up eingefangen, schaltet der Regler in die PFL-Stabilisierung:

```
v_center = −30·(x − x_soll) − 10·ẋ
v_angle = 100·θ₁ + 80·θ₂
v_damp = −25·θ̇₁ − 20·θ̇₂
v_anti = −1,5·(θ̇₂ − sign(θ̇₂)·8)  wenn |θ̇₂| > 8
         −1,0·(θ̇₁ − sign(θ̇₁)·7)  wenn |θ̇₁| > 7
v_desired = clamp(v_center + v_angle + v_damp + v_anti, −15, 15)
u = pflForce(v_desired, ...)
u = clamp(u, −50, 50)
```

Zentriert den Wagen auf `x_soll`, reguliert die Pendelwinkel und
dämpft die Geschwindigkeiten.

---

## 7. Schienenbegrenzung

Die Schiene hat harte Anschläge bei `x = ±2,5 m`. Bei Erreichen
wird der Wagen auf die Grenze gesetzt (`ẋ = 0`) und in den
Stabilisierungsmodus (Zustand 1) geschaltet.

Zusätzlich wird die Kraft nahe den Grenzen (±2,0 m bis ±2,5 m)
proportional reduziert, wenn der Wagen sich auf die Grenze zubewegt:

```
margin = 0,5 m
scale = max(0, (2,5 − |x|) / margin)
u *= scale  (falls ẋ in Richtung der Grenze)
```

---

## 8. Simulationsparameter

| Parameter | Steuerung | Bereich | Standard |
|-----------|-----------|---------|----------|
| Wagenmasse `M` | Slider | 0,5–5,0 kg | 1,50 kg |
| Pendel-1-Masse `m₁` | Slider | 0,05–2,00 kg | 0,50 kg |
| Pendel-2-Masse `m₂` | Slider | 0,05–3,50 kg | 0,50 kg |
| Schwerpunkt 1 `l₁` | Slider | 0,10–1,00 m | 0,20 m |
| Schwerpunkt 2 `l₂` | Slider | 0,10–1,00 m | 0,20 m |
| Gesamtlänge `L₁` | automatisch | — | 0,40 m (`L₁ = 2·l₁`) |
| Gesamtlänge `L₂` | automatisch | — | 0,40 m (`L₂ = 2·l₂`) |
| Gelenkdämpfung `c` | Slider | 0,00–0,50 Nms/rad | 0,00 Nms/rad |
| Gravitation `g` | Slider | 1,0–35,0 m/s² | 9,81 m/s² |
| Adaptation | Toggle | AN/AUS | AN |
| Schienenbegrenzung | Toggle | AN/AUS | AN |

*Hinweis: Die LQR-Verstärkungen werden bei Parameteränderungen automatisch
neu berechnet, sofern die Adaptation aktiviert ist.*

---

## 9. Szenarien

| Szenario | Anfangsbedingung | Erwartung |
|----------|-----------------|-----------|
| **A** – Stabilisierung | `x=0, θ₁=0,05, θ₂=−0,03` | LQR stabilisiert aufrechte Lage |
| **B** – Freier Fall | `x=0, θ₁=1,57, θ₂=0` | Pendel fallt, Swing-Up fangt |
| **C** – Swing-Up | `x=0, θ₁=π, θ₂=π` | PFL-Swing-Up schwingt auf |

---

## 10. Interaktion

| Aktion | Beschreibung |
|--------|-------------|
| **Power ON/OFF** | Regler ein-/ausschalten |
| **Reset State** | Zustand auf Szenario-Anfangswert zurücksetzen |
| **Shock Test** | Zufälligen Impuls auf Wagen und Pendel |
| **Download CSV** | Daten als CSV exportieren |
| **Adaptation ON/OFF** | Adaptive LQR-Verstärkung ein-/ausschalten |
| **Limits ON/OFF** | Schienenbegrenzung ein-/ausschalten |
| **Mausrad** | Animation zoomen (10 % – 300 %) |
| **Klick auf Schiene** | Sollposition `x_soll` setzen |
| **Zielwinkel-Buttons** | `(0,0)`, `(0,π)`, `(π,0)`, `(π,π)` setzen |

---

## 11. Datenausgabe (CSV)

Ein Klick auf „Download CSV" exportiert bis zu 10 000 Zeilen mit:

```
time, x, dx, theta1, dtheta1, theta2, dtheta2, force_N, Epot_J, Ekin_J
```

Die Zeit ist der Frame-Index. Alle Winkel sind auf
`(−π, π]` normiert. Die Energien sind relativ zur aufrechten Ruhelage.

---

## 12. Implementierungsdetails

- **Simulationsschrittweite:** `dt = 62,5 µs` (16 kHz)
- **LQR-Verstärkung:** Adaptiv aus Massenmatrix berechnet
- **Regelgesetz (LQR):** Zustandsrückführung mit 6 Verstärkungen
- **PFL:** Entkoppelt Wagenbeschleunigung von Pendeldynamik
- **Zustandsautomat:** 3 Zustande (LQR, Stabilisierung, Swing-Up)
- **Löser:** RK4 für die nichtlineare Simulation
- **Anzeige:** 6 Einzelplots mit 150 Punkten Gleitfenster
- **Sprache:** Reines JavaScript, keine Bibliotheken