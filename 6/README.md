# Experiment 6: Alternating Current (AC) Circuits

This directory contains the experimental data, phasor analysis routines, error propagation scripts, and LaTeX reports for the **Alternating Current (AC) Circuits** experiment in Physics Laboratory II at Sharif University of Technology.

## Overview

The primary objective of this experiment is to investigate the amplitude and phase relationships between voltage and current in series **RL**, **RC**, and **RLC** circuits driven by a sinusoidal AC voltage source at the line frequency ($f = 50.0\ \text{Hz}$).

Key physical phenomena investigated include:
- Ohm's law and complex impedance representations in the frequency domain.
- Vector phasor diagram construction using Kirchhoff's voltage law: $\tilde{V}_Z = \tilde{V}_R + \tilde{V}_L + \tilde{V}_C$.
- Extraction of non-ideal inductor parameters: decomposing total coil voltage $V_L$ into an in-phase ohmic resistive drop ($V_{r_L} = I r_L$) and a quadrature inductive drop ($V_{X_L} = I X_L$).
- Verification of capacitor ideality and determination of capacitance $C$.
- Series resonance prediction and verification that at $f = 50\ \text{Hz} > f_0 \approx 32.7\ \text{Hz}$, the series RLC circuit operates in the inductive regime ($X_L > X_C$).

---

## Theoretical Background

### 1. Sinusoidal Quantities & RMS Values

A sinusoidal alternating voltage is represented in the time domain as:
```
v(t) = V_m sin(omega * t + theta)
```
where $V_m$ is the peak amplitude, $\omega = 2\pi f$ is the angular frequency, and $\theta$ is the initial phase. Standard AC meters measure the root-mean-square (RMS) value:
```
V_rms = sqrt((1/T) * int_0^T v^2(t) dt) = V_m / sqrt(2)
I_rms = I_m / sqrt(2)
```

### 2. Complex Impedance & Phasors

Using Euler's formula $e^{i\omega t}$, steady-state currents and voltages are expressed as complex phasors:
```
v(t) = Im{ V_tilde * e^(i*omega*t) },    i(t) = Im{ I_tilde * e^(i*omega*t) }
```
The complex impedance of each element is defined by $\tilde{V} = \tilde{I} Z$:
- **Resistor ($R$):**
  ```
  Z_R = R,    v_R(t) is in-phase with i(t)
  ```
- **Ideal Inductor ($L$):**
  ```
  v_L(t) = L * (di/dt)  -->  Z_L = i * omega * L = i * X_L    (voltage leads current by +90 deg)
  ```
- **Real Inductor ($L, r_L$):**
  ```
  Z_coil = r_L + i * omega * L = r_L + i * X_L
  ```
- **Ideal Capacitor ($C$):**
  ```
  i_C(t) = C * (dv_C/dt)  -->  Z_C = 1 / (i * omega * C) = -i * X_C    (voltage lags current by -90 deg)
  where X_C = 1 / (omega * C)
  ```

### 3. Phasor Decomposition for Non-Ideal Inductors

Because a physical coil possesses winding resistance $r_L$, the voltage across the inductor $\vec{V}_L$ is not strictly perpendicular to $\vec{V}_R$. By the law of cosines in the voltage phasor triangle $\vec{V}_Z = \vec{V}_R + \vec{V}_L$:
```
V_Z^2 = V_R^2 + V_L^2 + 2 * V_R * V_L * cos(delta)
cos(delta) = (V_Z^2 - V_R^2 - V_L^2) / (2 * V_R * V_L)
```
The orthogonal components of coil voltage are:
```
V_rL = V_L * cos(delta) = (V_Z^2 - V_R^2 - V_L^2) / (2 * V_R)
V_XL = V_L * sin(delta) = sqrt(V_L^2 - V_rL^2)
```
From which the coil internal resistance $r_L$, inductive reactance $X_L$, and self-inductance $L$ are obtained:
```
r_L = V_rL / I
X_L = V_XL / I
L   = X_L / (2 * pi * f)
```

### 4. Series RLC Circuit & Resonance

For a series connection of $R$, $L$ (with $r_L$), and $C$:
```
Z_eq = (R + r_L) + i * (omega * L - 1 / (omega * C)) = R_tot + i * (X_L - X_C)
|Z| = sqrt(R_tot^2 + (X_L - X_C)^2)
phi = arctan((X_L - X_C) / R_tot)
```
The undamped series resonant frequency occurs when the imaginary part vanishes ($X_L = X_C$):
```
omega_0 * L = 1 / (omega_0 * C)  -->  f_0 = 1 / (2 * pi * sqrt(L * C))
```
- For $f < f_0$: $X_C > X_L \implies$ Circuit is capacitive ($\varphi < 0$).
- For $f = f_0$: $X_L = X_C \implies |Z| = R_\text{tot}$, current is maximal and in-phase with source.
- For $f > f_0$: $X_L > X_C \implies$ Circuit is inductive ($\varphi > 0$).

---

## Experimental Setup

The circuit is powered by a variable AC voltage supply adjusted to provide a stable input voltage ($\approx 5.34\ \text{V}$ to $5.37\ \text{V}$ RMS) at mains line frequency ($f = 50.0 \pm 0.1\ \text{Hz}$).

```
       +---( ~ AC Source )---+
       |                     |
       +---[  Resistor R ]---+
       |                     |
       +---[ Inductor  L ]---+  (in RL / RLC)
       |   [ + coil r_L  ]   |
       +---[ Capacitor C ]---+  (in RC / RLC)
```

**Instruments Used:**
- AC Variable Power Supply ($50\ \text{Hz}$, $0\text{--}20\ \text{V}$)
- Calibrated Digital Multimeter (AC Voltmeter mode, resolution $\pm 0.01\ \text{V}$)
- Calibrated Digital Multimeter (AC Ammeter mode, resolution $\pm 0.01\ \text{mA}$)
- Decade / fixed power resistor ($R \approx 104\ \Omega$)
- Heavy iron-core induction coil ($L \approx 1.18\ \text{H}$, $r_L \approx 40\text{--}50\ \Omega$)
- AC bipolar capacitor ($C \approx 20\ \mu\text{F}$)

---

## Key Analyses

1. **Cross-Circuit Resistor Linearity (Ohm's Law Verification):**
   - Regressing $V_R$ against current $I$ across the three configurations (RL: $13.54\ \text{mA}$, RLC: $20.04\ \text{mA}$, RC: $28.29\ \text{mA}$) using Ordinary Least Squares (OLS).
   - Confirms high linear correlation ($R^2 = 0.99999$) and yields $R = 102.35 \pm 0.36\ \Omega$.
2. **Series RL Phasor Analysis:**
   - Determination of the non-orthogonality angle $\delta = 83.89^\circ \pm 0.62^\circ$, proving the presence of internal winding resistance $r_L = 39.52 \pm 3.96\ \Omega$.
   - Extraction of $L = 1.1758 \pm 0.0037\ \text{H}$ and phase lead $\varphi = 68.65^\circ \pm 0.54^\circ$.
3. **Series RC Phasor Analysis:**
   - Calculation of $\delta_{RC} = 89.91^\circ \pm 0.33^\circ \approx 90^\circ$, demonstrating ideal capacitive behavior and negligible Equivalent Series Resistance (ESR).
   - Extraction of capacitance $C = 20.19 \pm 0.06\ \mu\text{F}$ and phase lag $\varphi = -56.70^\circ \pm 0.11^\circ$.
4. **Series RLC Synthesis & Cross-Validation:**
   - Complete vector loop evaluation yielding independent estimates $L = 1.1945 \pm 0.0032\ \text{H}$ and $C = 19.81 \pm 0.07\ \mu\text{F}$.
   - Comparison between measured impedance $|Z|_\text{meas} = 266.47 \pm 0.52\ \Omega$ and cross-predicted theoretical impedance $|Z|_\text{cross} = 255.55 \pm 2.37\ \Omega$ (agreement within $4.1\%$).
   - Identification of natural resonance at $f_0 = 32.66 \pm 0.07\ \text{Hz}$, explaining why the circuit behaves inductively at $50\ \text{Hz}$.

---

## Results

### Summary of Component Parameters and Impedances

| Quantity | Symbol | Measured / Derived Value | Unit | Method |
| :--- | :---: | :---: | :---: | :--- |
| **Resistor (OLS slope)** | $R$ | $102.35 \pm 0.36$ | $\Omega$ | Linear fit $V_R$ vs $I$ ($R^2 = 0.99999$) |
| **Resistor (RL circuit)** | $R_{RL}$ | $104.87 \pm 0.74$ | $\Omega$ | $V_R / I$ |
| **Resistor (RC circuit)** | $R_{RC}$ | $103.57 \pm 0.36$ | $\Omega$ | $V_R / I$ |
| **Resistor (RLC circuit)** | $R_{RLC}$ | $104.29 \pm 0.50$ | $\Omega$ | $V_R / I$ |
| **Coil Internal Resistance (RL)** | $r_{L,RL}$ | $39.52 \pm 3.96$ | $\Omega$ | Phasor decomposition |
| **Coil Internal Resistance (RLC)** | $r_{L,RLC}$ | $51.13 \pm 2.76$ | $\Omega$ | Phasor decomposition from $V_{RL}$ |
| **Inductive Reactance (RL)** | $X_{L,RL}$ | $369.38 \pm 0.90$ | $\Omega$ | $V_{X_L} / I$ |
| **Self-Inductance (RL)** | $L_{RL}$ | $1.1758 \pm 0.0037$ | $\text{H}$ | $X_L / (2\pi f)$ |
| **Self-Inductance (RLC)** | $L_{RLC}$ | $1.1945 \pm 0.0032$ | $\text{H}$ | $X_L / (2\pi f)$ |
| **Capacitive Reactance (RC)** | $X_{C,RC}$ | $157.65 \pm 0.36$ | $\Omega$ | $V_C / I$ |
| **Capacitance (RC)** | $C_{RC}$ | $20.19 \pm 0.06$ | $\mu\text{F}$ | $1 / (2\pi f X_C)$ |
| **Capacitance (RLC)** | $C_{RLC}$ | $19.81 \pm 0.07$ | $\mu\text{F}$ | $1 / (2\pi f X_C)$ |
| **Resonant Frequency** | $f_0$ | $32.66 \pm 0.07$ | $\text{Hz}$ | $1 / (2\pi\sqrt{LC})$ |

### Circuit Impedances & Phase Angles

| Configuration | $|Z|_\text{meas}\ [\Omega]$ | $|Z|_\text{th}\ [\Omega]$ | Relative Diff. | Phase Angle $\varphi$ | Behavior |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Series RL** | $396.60 \pm 0.79$ | $396.60 \pm 1.69$ | $0.00\%$ | $+68.65^\circ \pm 0.54^\circ$ | Inductive ($V$ leads $I$) |
| **Series RC** | $188.76 \pm 0.36$ | $188.63 \pm 0.36$ | $0.07\%$ | $-56.70^\circ \pm 0.11^\circ$ | Capacitive ($I$ leads $V$) |
| **Series RLC** | $266.47 \pm 0.52$ | $264.97 \pm 1.78$ | $0.56\%$ | $+54.09^\circ \pm 0.50^\circ$ | Inductive ($f > f_0$) |

*(Note: When cross-predicting $|Z|_\text{RLC}$ using parameters strictly measured in the isolated RL and RC circuits, $|Z|_\text{cross} = 255.55 \pm 2.37\ \Omega$, yielding a discrepancy of only $4.1\%$.)*

---

## Files

- `6.tex` — Original concise Persian LaTeX lab report.
- `6-1.tex` — Comprehensive extended Persian LaTeX report with complete derivations and TikZ phasor plots.
- `6.pdf` / `6-1.pdf` — Compiled PDF reports.
- `analysis.py` — Python numerical analysis script containing OLS regressions, full partial-derivative error propagation, and vector plot generation.
- `plots/` — Directory containing generated vector PDF figures:
  - `VR_vs_I_resistor.pdf` — Linear regression of $V_R$ vs $I$ verifying Ohm's law.
  - `phasor_RL.pdf` — Detailed vector diagram for the RL circuit with coil decomposition.
  - `phasor_RC.pdf` — Phasor diagram for the RC circuit showing orthogonal capacitive lag.
  - `phasor_RLC.pdf` — Full phasor diagram for the series RLC network.
  - `phasor_summary.pdf` — Multi-panel side-by-side comparison of all three phasors.
  - `impedance_comparison.pdf` — Bar chart comparing measured and theoretical impedances.
  - `resonance_curve.pdf` — Theoretical frequency response curves ($I(f)$ and $|Z|(f)$) showing resonance at $f_0 = 32.7\ \text{Hz}$.
- `task1.md`, `task2.md`, `task3.md` — Experimental assignment prompts and guide notes.
- `README.md` — This documentation file.

---

## How to Compile

To compile the LaTeX reports, run `xelatex` (required for `xepersian` and Persian font rendering):

```bash
xelatex 6-1.tex
xelatex 6-1.tex
```

To re-run the Python analysis and regenerate all figures:

```bash
python analysis.py
```

---

## Notes

- **Why is $V_L$ not perpendicular to $V_R$?** A real inductor consists of many meters of copper wire wound around a ferromagnetic core. The finite conductivity of copper introduces a non-zero DC/AC winding resistance $r_L \approx 40\ \Omega$. Consequently, $V_L$ has both an in-phase ohmic component ($V_{r_L} = I r_L$) and a quadrature inductive component ($V_{X_L} = I \omega L$). The angle is $\delta \approx 83.9^\circ < 90^\circ$.
- **Why can capacitor ESR not be measured with this method?** In modern capacitors, the equivalent series resistance (ESR) is negligible ($\ll 1\ \Omega$), producing an in-phase voltage drop of less than a few millivolts. This drop is completely buried within the digital voltmeter quantization error ($\pm 10\ \text{mV}$). Measuring ESR requires dedicated AC bridges or high-frequency impedance spectroscopy.
- **Inductor Discrepancy ($r_L = 39.5\ \Omega$ vs $51.1\ \Omega$):** The apparent increase in coil resistance in the RLC circuit arises from Joule heating ($I^2 r_L$) during sustained operation, magnetic core hysteresis/eddy-current losses, and small loading effects of the measurement probes.
