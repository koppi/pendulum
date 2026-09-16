# Doppel-Invertpendulum — Regelung im Browser

Diese Simulation bringt ein doppeltes invertiertes Pendel auf einem fahrbaren
Wagen aus der hängenden Ruhelage in die aufrechte **Up-Up-Lage** und hält es
dort. Der Regler besteht aus drei Teilen:

* **LQR-Zustandsrückführung** um die obere Ruhelage, als echte Lösung der
  algebraischen Riccati-Gleichung (Matrix-Signumfunktion),
* **Aufschwingen über Trajektorienoptimierung** (iLQR auf dem vollständigen
  nichtlinearen Modell), geführt von einem **zeitvarianten LQR**,
* **Beruhigung** durch Energieentzug, die das System vor jedem Versuch in den
  Startzustand der geplanten Trajektorie bringt.

Der Übergabepunkt zwischen Aufschwingen und Balancieren ist keine Heuristik,
sondern eine Niveaumenge der Riccati-Matrix.

Die gesamte Berechnung — Riccati-Löser und Trajektorienoptimierung
eingeschlossen — läuft in JavaScript im Browser, ohne Bibliotheken, in
Echtzeit.

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

Potenzielle Energie (relativ zur aufrechten Ruhelage). Pendel 1 trägt die
Masse von Pendel 2 an seiner Spitze, weshalb `m₂·L₁` mit eingeht:

```
Epot = (m₁·l₁ + m₂·L₁)·g·(cos(θ₁)−1) + m₂·l₂·g·(cos(θ₂)−1)
```

Kinetische Energie der Glieder — mit dem Kopplungsterm `M₂₃`:

```
Ekin = ½·(M₂₂·θ̇₁² + 2·M₂₃·θ̇₁·θ̇₂ + M₃₃·θ̇₂²)
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
M₁₁·ẍ + M₁₂·θ̈₁ + M₁₃·θ̈₂ + C₁ + G₁ = u
```

Die Pendelbeschleunigungen `θ̈₁, θ̈₂` hängen affin von `v = ẍ` ab:

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

Dabei sind `C₁, C₂, C₃` die Komponenten von `C(q,q̇)·q̇` und damit
**quadratisch** in den Gelenkgeschwindigkeiten:

```
C₁ = −(m₁·l₁+m₂·L₁)·sin(θ₁)·θ̇₁²  −  m₂·l₂·sin(θ₂)·θ̇₂²
C₂ =  m₂·L₁·l₂·sin(θ₁−θ₂)·θ̇₂²
C₃ = −m₂·L₁·l₂·sin(θ₁−θ₂)·θ̇₁²
```

### 2.2 PFL-Stellgesetz

```
u = α·v + φ

α = M₁₁ + M₁₂·θ̈₁ᵛ + M₁₃·θ̈₂ᵛ
φ = M₁₂·θ̈₁ᵖ + M₁₃·θ̈₂ᵖ + C₁ + G₁
```

Damit wird jeder gewünschte `v`-Wert exakt in eine Kraft `u` umgesetzt.

### 2.3 Energiebilanz

Mit `v` als Eingang gilt für die Energie des Pendelpaares

```
E   = ½·θ̇ᵀ·M_p·θ̇ + a₁·g·cos(θ₁) + a₂·g·cos(θ₂)
a₁  = m₁·l₁ + m₂·L₁,   a₂ = m₂·l₂
M_p = [[M₂₂, M₂₃], [M₂₃, M₃₃]]
```

exakt die Beziehung

```
Ė = −v·W − c·(θ̇₁² + θ̇₂²),   W = a₁·cos(θ₁)·θ̇₁ + a₂·cos(θ₂)·θ̇₂
```

(die Coriolis-Terme heben sich gegen `½·θ̇ᵀ·Ṁ_p·θ̇` auf). `W` ist der
einzige Kanal, über den der Wagen Energie zu- oder abführen kann —
diese Identität trägt sowohl die Beruhigung (Abschnitt 6) als auch die
Energiepumpen-Rückfallebene.

---

## 3. Zustandsautomat

Der Regler arbeitet in drei Zuständen:

| `recoveryState` | Modus | Beschreibung |
|------------------|-------|--------------|
| 0 | `LQR Balance` | Zustandsrückführung um die obere Ruhelage |
| 1 | `Settling` | Energieentzug bis zur hängenden Ruhelage, Wagen zentrieren |
| 2 | `Swing-up` | Abfahren der optimierten Aufschwing-Trajektorie |

### Zustandsübergänge

```
Power ON / Reset → 0, falls der Zustand im Einfangbereich liegt, sonst → 1
Zustand 1: hängende Ruhelage erreicht ∧ Trajektorie vorhanden → 2
Zustand 2: Einfangbereich erreicht → 0
Zustand 2: Trajektorie abgelaufen ohne Einfang → 1 (neuer Versuch)
Zustand 0: |θ₁|>0,8 oder |θ₂|>0,8 für >20 Zyklen → 1
Schienen-Anschlag (±2,5 m): → 1
```

Der Einfangbereich ist **nicht** heuristisch, sondern wird aus der
Riccati-Matrix des LQR abgeleitet (Abschnitt 4.3).

---

## 4. LQR-Zustandsrückführung (RecoveryState = 0)

### 4.1 Linearisierung

Bei `θ₁ = θ₂ = 0, q̇ = 0` ist `M(q)` konstant, `C` verschwindet und `G`
ist linear. Mit `N = M₀⁻¹` ergibt sich

```
q̈ = N·([u,0,0]ᵀ − ∂G/∂q·q − D·q̇)
```

und daraus das Zustandsmodell `ẋ = A·x + B·u` mit
`x = [x, ẋ, θ₁, θ̇₁, θ₂, θ̇₂]ᵀ`. Die offene Strecke hat für die
Standardparameter die Eigenwerte

```
0, 0, ±4,96, ±11,44   (1/s)
```

— zwei instabile Pole, wie für ein doppeltes Invertpendel erwartet.

### 4.2 Lösung der Riccati-Gleichung

Die Verstärkung wird als echte LQR-Lösung berechnet, nicht geschätzt.
Gelöst wird

```
AᵀP + P·A − P·B·R⁻¹·Bᵀ·P + Q = 0
```

über die **Matrix-Signumfunktion** des Hamilton-Operators

```
H = [[A, −B·R⁻¹·Bᵀ], [−Q, −Aᵀ]]
```

Die Newton-Iteration `Z ← ½·(c·Z + (c·Z)⁻¹)` konvergiert quadratisch
gegen `sign(H)`; der stabile invariante Unterraum wird von den Spalten
von `I − sign(H)` aufgespannt und hat die Form `[u; P·u]`, woraus `P`
per Ausgleichsrechnung folgt. Einige Newton/Kleinman-Schritte
(Lyapunov-Gleichung als 36×36-System) polieren `P` auf Maschinen-
genauigkeit. Das Residuum liegt bei `≈10⁻¹⁰`, die Rechenzeit bei
wenigen Millisekunden — die Synthese läuft daher bei jeder
Parameteränderung neu.

Gewichte:

```
Q = diag(5,0 ; 0,1 ; 10 ; 50 ; 10 ; 50)
R = 0,5
```

Die geringe Gewichtung der Wagenposition und die hohe Gewichtung der
Winkelgeschwindigkeiten maximieren den Einzugsbereich, der innerhalb
der Stellgrößenbegrenzung von ±50 N erreichbar ist — und genau diesen
Bereich muss der Aufschwingvorgang treffen.

Für die Standardparameter ergibt sich

```
K ≈ [3,16 ; 6,49 ; −303,4 ; −8,91 ; 343,9 ; 41,4]
u  = −K·(x − x_soll-Vektor),   u = clamp(u, −50, 50)
```

mit den Regelkreispolen `−24,60`, `−8,10 ± 4,31i`, `−1,70`, `−0,71 ± 0,81i`.

Gemessener Einzugsbereich (nichtlineare Simulation, ±50 N):

| Auslenkung | Grenze |
|------------|--------|
| `θ₁ = θ₂`  | ±0,59 rad |
| nur `θ₁`   | ±0,36 rad |
| nur `θ₂`   | ±0,41 rad |
| `θ̇₁ = θ̇₂`  | ±2,1 rad/s |

### 4.3 Einfangbereich

Übergeben wird an den LQR, sobald

```
eᵀ·P·e < catchLevel   ∧   |θ₁|,|θ₂| < 0,6   ∧   |θ̇₁|,|θ̇₂| < 4,0   ∧   |x−x_soll| < 2,0
```

Die expliziten Schranken sind nur eine Plausibilitätssicherung (kein Einfang
auf der falschen Seite eines Winkelumlaufs oder direkt am Anschlag) —
entscheidend ist die Niveaumenge.

Die Wagenposition geht **absichtlich nicht** in die quadratische Form ein: Im
Moment der Übergabe wird die LQR-Sollposition auf die aktuelle Wagenposition
gesetzt und danach mit einer Zeitkonstante von 1,5 s auf `x_soll` zurückgeführt.
So steht das gesamte Kraftbudget für das Einfangen der Glieder zur Verfügung
statt für einen Positionsfehler, der ohnehin gleich verschwindet.

Die Schwelle ist die größte Niveaumenge von `eᵀPe`, auf der die
unbegrenzte LQR-Kraft die Stellgrenze noch einhält:

```
max |K·e| über eᵀPe ≤ c  =  √(c · K·P⁻¹·Kᵀ)
⟹ catchLevel = U_max² / (K·P⁻¹·Kᵀ)     (≈ 28,5 für die Standardparameter)
```

Monte-Carlo-Stichproben (je ~100 akzeptierte Zustände, anschließend 12 s
nachsimuliert) bestätigen die Menge für die meisten Strecken vollständig. Zwei
Ausreißer blieben: bei `l₂ = 0,5 m` ging 1 von 106 akzeptierten Zuständen
verloren, bei `g = 1 m/s²` 4 von 135 — dort ist die Schranke also um wenige
Prozent optimistisch, weil sie aus dem *linearisierten* Modell stammt.

Zwei Dinge fangen das ab:

* Die **geplante** Aufschwingbahn hängt nicht von dieser Schranke ab. Sie wird
  vor der Übernahme komplett durchsimuliert — Aufschwingen *und* anschließendes
  Balancieren — und nur akzeptiert, wenn der LQR das Pendel danach tatsächlich
  aufrecht zur Ruhe bringt (Abschnitt 5.3).
* Greift der Test nach einer **Störung** doch einmal daneben, erkennt die
  Divergenzüberwachung das binnen 20 Rechenzyklen, der Regler beruhigt und
  schwingt erneut auf.

---

## 5. Aufschwingen durch Trajektorienoptimierung (RecoveryState = 2)

### 5.1 Warum kein reines Energiepumpen

Energiepumpen (`v ∝ (E−E_d)·W`) treibt die **Gesamtenergie** zuverlässig
auf ihren Wert in der oberen Ruhelage, sagt aber nichts darüber aus,
**wie** sich diese Energie auf die beiden Glieder verteilt. Die Menge
`E = E_d` ist fünfdimensional, der Einfangbereich darin verschwindend
klein. Messung an der Standardstrecke: in 120 s Energiepumpen liegen
beide Winkel zusammen nur **0,06 s** lang unterhalb von 0,5 rad, und
der beste erreichte Wert von `eᵀPe` betrug 178 gegenüber einer Schwelle
von 28,5. Das Aufschwingen gelang so in weniger als einem Drittel der
Versuche, mit Medianzeiten um 24 s.

### 5.2 iLQR-Trajektorienoptimierung

Stattdessen wird das Aufschwingen als Optimalsteuerungsproblem gelöst —
**iLQR** auf dem vollständigen nichtlinearen Modell, ausgehend von der
hängenden Ruhelage:

```
min  ½·w_f·e_Nᵀ·P·e_N  +  Σ_k [ ½·r·w_k² + ½·q_x·x_k²
                               + ρ(k)·(a₁·(1−cosθ₁) + a₂·(1−cosθ₂))
                               + ½·ρ_v(k)·(θ̇₁² + θ̇₂²)  + Schienenbarriere ]
```

| Element | Wert / Zweck |
|---------|--------------|
| Stützstellen `N` | 125 |
| Horizont `T` | `2,5 s · 4,22 / ω_hang`, begrenzt auf 1,5–6,0 s |
| Kraft | `u = 50·tanh(w)` — die Stellgrenze gilt konstruktiv, das Problem bleibt unbeschränkt |
| Endkosten | `eᵀPe` = LQR-Restkosten, also genau das Maß des Einfangbereichs |
| `ρ(k), ρ_v(k)` | „Aufrichten"-Formterme mit Gewicht `∝ (k/N)³` |
| Schienenbarriere | weich ab `|x| > 1,2 m` |

Der Horizont skaliert mit der langsamen Eigenfrequenz `ω_hang` des
hängenden Doppelpendels (aus `det(K_g − ω²M_p) = 0`), damit lange oder
schwach beschleunigte Pendel proportional mehr Zeit bekommen. Da `N`
fest bleibt, bleiben auch die Kosten eines Optimierungsschrittes gleich.

Der Endkostenterm allein hat einen zu kleinen Einzugsbereich für den
Optimierer; die glatten „Aufrichten"-Terme führen ihn dorthin. Gemessen
über je fünf Zufallsstarts mit identischem Startrauschen: ohne die
Formterme erreichte ein Start Endkosten unter 1, mit ihnen drei. Über
14 Starts mit den Formtermen erreichten neun eine brauchbare Lösung.

**Zeitscheiben:** Die Optimierung läuft in Scheiben von ≤ 8 ms je
Animationsbild, während der Regler noch beruhigt — die Oberfläche
bleibt flüssig. Ein Neustart, der nach 30 Iterationen noch weit entfernt
ist (mehr als das 22-fache des Einfangbereichs), wird verworfen statt
austherapiert; die letzten Versuche einer Runde laufen dagegen immer
voll durch. Die Neustarts wechseln zwischen 125 und bis zu 220
Stützstellen ab: das grobe Gitter konvergiert bei den Standardparametern
häufiger, sehr lange Glieder liefern brauchbare Bahnen dagegen nur auf
dem feinen.

### 5.3 Abnahme durch geschlossene Simulation

Eine Aufschwing-Trajektorie ist **offen instabil**. Die Endkosten der
Optimierung allein sagen daher wenig aus: dieselben Kräfte auf einem
anderen Zeitgitter nachintegriert laufen weg — gemessen wurde eine
Kandidatin mit `eᵀPe = 0,1` auf ihrem eigenen Gitter, die offen
nachgerechnet bei `eᵀPe > 12 000` endete.

Maßgeblich ist das **geregelte** Verhalten. Jede Kandidatin wird deshalb
vor der Übernahme einmal vollständig durchsimuliert — Schrittweite 0,5 ms,
Stellgrenze und Schienenabstand aktiv:

1. Aufschwingen entlang der Bahn mit dem zeitvarianten LQR, bis der
   Einfangtest anspricht;
2. anschließend acht Sekunden mit dem balancierenden LQR, inklusive
   Rückführung der Wagen-Sollposition.

Übernommen wird nur, was beides besteht — das Pendel muss am Ende aufrecht
und in Ruhe sein. Damit hängt die ausgelieferte Trajektorie nicht mehr davon
ab, ob die Einfangschranke im Grenzfall exakt stimmt. Der angezeigte Wert
`eᵀPe` ist der beim Einfang erreichte. Eine abgelehnte Kandidatin wird
verworfen und weitergesucht.

Schlägt ein Aufschwingversuch trotzdem fehl (Trajektorie abgelaufen oder
Schienenanschlag), startet im Hintergrund eine weitere Suchrunde, während
der Regler erneut beruhigt — bis zu sechs Runden.

### 5.4 Zeitvariantes LQR zur Bahnführung

Die geplante Trajektorie wird nicht gesteuert, sondern **geregelt**
abgefahren. Rückwärts entlang der Bahn wird die diskrete
Riccati-Rekursion gelöst

```
S_k = Q_t + AₖᵀS_{k+1}Aₖ − AₖᵀS_{k+1}Bₖ·(R_t + BₖᵀS_{k+1}Bₖ)⁻¹·BₖᵀS_{k+1}Aₖ
K_k = (R_t + BₖᵀS_{k+1}Bₖ)⁻¹·BₖᵀS_{k+1}Aₖ
```

mit `S_N = P` (der LQR-Matrix), `Q_t = diag(1,1,20,2,20,2)`, `R_t = 0,05`.
`Aₖ, Bₖ` folgen aus finiten Differenzen eines RK4-Schrittes. Gestellt wird

```
u = u_nom[k] − K_k·(x − x_nom[k]),   u = clamp(u, −50, 50)
```

Erst diese Rückführung macht das Verfahren robust: Der Plan wird mit
`T/125 ≈ 20 ms` gerechnet, die Simulation läuft mit 62,5 µs. Gemessene
Toleranz gegenüber Abweichungen im Startzustand: ±0,4 m, ±0,4 m/s,
±0,4 rad und ±0,4 rad/s werden noch eingefangen (15/15 Versuche).

---

## 6. Beruhigung (RecoveryState = 1)

Die Trajektorie wird aus der hängenden Ruhelage geplant. Dorthin bringt
der Regler das System mit der Energiebilanz aus Abschnitt 2.3: `v = +k_w·W`
liefert `Ė = −k_w·W² ≤ 0`, entzieht also bedingungslos Energie.

Das allein bleibt jedoch stecken. Das hängende Doppelpendel hat zwei
Eigenmoden — für die Standardparameter bei 4,22 rad/s mit der Form
`[1; 1,448]` (gleichphasig) und 10,94 rad/s mit `[1; −2,073]`
(gegenphasig). Deren Eingangskopplung `M_pxᵀ·u` beträgt −0,4448 bzw.
−0,0927: die gegenphasige Mode ist rund **fünfmal schwächer** an den
Wagen gekoppelt. Sie klingt daher kaum ab, während `W ≈ 0` bleibt.

Ein zweiter Term dämpft genau diese Mode direkt, indem die Kopplung
`∂(θ̈₁−θ̈₂)/∂v` invertiert wird:

```
v = k_w·W
    − clamp( (k_rel·(θ̇₁−θ̇₂)) / (θ̈₁ᵛ − θ̈₂ᵛ), ±8 )
    − k_x·(x − x_soll) − k_d·ẋ
v = clamp(v, −12, 12)

k_w = 6,0   k_rel = 3,0   k_x = 1,5   k_d = 1,5
```

Damit erreicht das System die Ruhelage aus allen geprüften Anfangs-
zuständen in 1–4 s, ohne die Schienenanschläge zu berühren.

---

## 7. Schienenbegrenzung

Die Schiene hat harte Anschläge bei `x = ±2,5 m`. Bei Erreichen
wird der Wagen auf die Grenze gesetzt (`ẋ = 0`) und in den
Beruhigungsmodus (Zustand 1) geschaltet. Der Test läuft in jedem
Integrationsschritt, nicht einmal pro Bild, sodass der Wagen die Grenze nicht
überfahren kann.

Die geplante Trajektorie hält von sich aus Abstand: die weiche Schienenbarriere
der Optimierung greift schon ab `|x| > 1,2 m`, und in allen geprüften Fällen
blieb die Bahn unter 1,25 m.

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

*Hinweis: Bei „Adaptation: AN" werden die LQR-Verstärkungen bei jeder
Parameteränderung neu berechnet. Bei „AUS" bleibt die Verstärkung auf dem
Stand des Abschaltzeitpunkts eingefroren — man kann so beobachten, wie der
Regler auf eine falsch modellierte Strecke reagiert. Die Riccati-Matrix `P`
(und damit Einfangtest und Endkosten der Planung) folgt in beiden Fällen der
tatsächlichen Strecke.*

*Jede Parameteränderung verwirft außerdem die Aufschwing-Trajektorie und
startet die Planung neu — sie läuft im Hintergrund weiter, während der Regler
balanciert oder beruhigt.*

### Geprüfter Arbeitsbereich

Jeder Lauf startet aus der hängenden Ruhelage und gilt als bestanden, wenn das
Pendel danach in der Up-Up-Lage gehalten wird (`|θ| < 0,03 rad`, Wagen auf
`x_soll`):

| Parameter | geprüft | Ergebnis |
|-----------|---------|----------|
| `M` | 0,5 – 5,0 kg | vollständig |
| `m₁` | 0,05 – 2,0 kg | vollständig |
| `m₂` | 0,05 – 3,5 kg | vollständig |
| `l₁` | 0,10 – 1,00 m | vollständig |
| `l₂` | 0,10 – 0,50 m | vollständig |
| `g` | 1,0 – 20,0 m/s² | vollständig |
| `c` | 0,00 – 0,20 Nms/rad | vollständig |

Ebenfalls geprüft: alle drei Szenarien, die vier Zielwinkel-Knöpfe, Störimpulse
(bis `Δẋ = −2 m/s`, `Δθ̇₁ = Δθ̇₂ = −3 rad/s` aus der balancierten Lage) und das
Verschieben von `x_soll` während des Balancierens.

**Bekannte Grenzen.** An drei Enden der Schieberegler findet die Optimierung
innerhalb ihres Versuchsbudgets keine Bahn, die die Abnahme besteht:

* `l₂ = 1,0 m` — ein 2 m langes oberes Glied,
* `g = 35 m/s²` — größte einstellbare Gravitation,
* `c = 0,50 Nms/rad` — größte Gelenkdämpfung.

Das Verhalten bleibt dort gutartig: der Regler beruhigt das Pendel, plant im
Hintergrund weiter und versucht es erneut; die Anzeige steht auf
„No plan — energy pumping". Das **Balancieren** selbst — Szenario A und die
Störungsversuche — funktioniert auch bei diesen Einstellungen.

---

## 9. Szenarien

| Szenario | Anfangsbedingung | Erwartung |
|----------|-----------------|-----------|
| **A** – Stabilisierung | `x=0, θ₁=0,05, θ₂=−0,03` | LQR hält direkt; aufrecht nach ≈3 s |
| **B** – Freier Fall | `x=0, θ₁=1,57, θ₂=0` | Pendel fällt, wird beruhigt, schwingt auf; aufrecht nach ≈8 s |
| **C** – Swing-Up | `x=0, θ₁=π, θ₂=π` | Aufschwingen aus der Ruhelage; aufrecht nach ≈3 s (Plan bereit) |

In allen drei Fällen wird die Up-Up-Lage auf besser als `10⁻³ rad` gehalten,
der Wagen kehrt auf `x_soll` zurück.

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

- **Simulationsschrittweite:** `dt = 62,5 µs` (16 kHz), RK4
- **Stellgrenze:** `|u| ≤ 50 N`, Schiene `|x| ≤ 2,5 m`
- **LQR:** exakte CARE-Lösung über die Matrix-Signumfunktion des
  Hamilton-Operators, nachpoliert mit Newton/Kleinman; Residuum `≈10⁻¹⁰`,
  Rechenzeit wenige Millisekunden
- **Einfangtest:** `eᵀPe < U_max²/(K·P⁻¹·Kᵀ)` — die größte Niveaumenge ohne
  Stellgrößensättigung
- **Aufschwingen:** iLQR, 125 Stützstellen, Horizont an die langsame
  Eigenfrequenz des hängenden Pendels gekoppelt, Kraft über `tanh` beschränkt
- **Bahnführung:** zeitvariantes LQR aus der diskreten Riccati-Rekursion
  entlang der Bahn, Endgewicht `S_N = P`
- **Planung im Hintergrund:** ≤ 8 ms je Animationsbild, mindestens ein
  Optimierungsschritt pro Bild; Neuplanung bei jeder Parameteränderung
- **PFL:** setzt eine gewünschte Wagenbeschleunigung exakt in eine Kraft um;
  Coriolis-Terme quadratisch in den Gelenkraten
- **Zustandsautomat:** 3 Zustände (Balance, Beruhigung, Aufschwingen)
- **Anzeige:** 6 Einzelplots mit 150 Punkten Gleitfenster
- **Sprache:** Reines JavaScript, keine Bibliotheken

### Numerik im Detail

| Baustein | Verfahren | Größe |
|----------|-----------|-------|
| Matrixinversion | Gauß-Jordan mit Spaltenpivotisierung | bis 12×12 |
| `sign(H)` | Newton-Iteration mit Normskalierung | 12×12 |
| Lyapunov-Gleichung | direktes lineares System | 36×36 |
| iLQR-Jacobi-Matrizen | zentrale Differenzen eines RK4-Schrittes | 6×6, 6×1 |
| Liniensuche | Rückwärtssuche über 9 Schrittweiten | — |