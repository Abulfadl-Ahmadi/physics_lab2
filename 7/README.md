# Experiment 7: Capacitor Characteristics & Dielectric Constant Measurement (with Oscilloscope Analysis)

This directory contains the experimental data, analysis scripts, generated plots, and LaTeX reports for **Experiment 7** in the Physics Lab II curriculum at Sharif University of Technology.

---

## Overview

The purpose of this experiment is to systematically investigate the properties of a circular parallel-plate capacitor under alternating current (AC) excitation. The primary physical objectives include:
1. **Oscilloscope Waveform & Frequency Verification:** Measuring signal amplitude, period, and frequency of a function generator using an oscilloscope, and evaluating reading uncertainties.
2. **Capacitance & Dielectric Constant of Plexiglass ($K_{\mathrm{plexi}}$):** Measuring the AC current as a function of frequency ($1\text{ to }25\text{ kHz}$) with a plexiglass sheet between the plates to extract capacitance and the relative permittivity.
3. **Capacitance & Permittivity of Air ($\varepsilon_{\mathrm{air}}$):** Repeating frequency measurements with air as dielectric to deduce $\varepsilon_{\mathrm{air}}$, comparing with vacuum permittivity $\varepsilon_0$, and calculating the geometry-independent capacitance ratio $C_{\mathrm{plexi}}/C_{\mathrm{air}}$.
4. **Distance Dependence ($C \propto 1/d$):** Varying the plate separation $d$ ($3\text{ to }9\text{ mm}$) at fixed frequency ($14\text{ kHz}$) to verify the inverse relationship and assess parasitic (stray) capacitance.

---

## Theoretical Background

### 1. Alternating Current through an Ideal Capacitor

When a sinusoidal AC voltage $V(t) = V_m \sin(\omega t)$ is applied across an ideal capacitor with capacitance $C$, the stored electric charge is:

```
q(t) = C · V(t) = C · V_m sin(ωt)
```

Differentiating with respect to time yields the instantaneous current:

```
I(t) = dq/dt = ω · C · V_m cos(ωt) = ω · C · V_m sin(ωt + π/2)
```

The current leads the voltage by a phase angle of $\pi/2$ ($90^\circ$). The peak current amplitude is $I_m = \omega C V_m = 2\pi f C V_m$. In terms of effective root-mean-square (RMS) quantities:

```
V_rms = V_m / √2,      I_rms = I_m / √2
I_rms = 2π · f · C · V_rms
```

For constant $V_{\mathrm{rms}}$ and fixed geometry, $I_{\mathrm{rms}}$ is strictly proportional to frequency $f$. The slope of $I$ versus $f$ provides a direct measurement of capacitance:

```
m = dI/df = 2π · C · V_rms   ===>   C = m / (2π · V_rms)
```

### 2. Parallel-Plate Geometry and Relative Permittivity

The capacitance of a parallel-plate capacitor with plate area $A = \pi r^2$ and separation $d$ filled with a dielectric of relative permittivity $K$ is:

```
C = (K · ε_0 · A) / d
```

From the measured capacitance $C$, the dielectric constant is obtained as:

```
K = (C · d) / (ε_0 · A)
```

where $\varepsilon_0 = 8.8541878 \times 10^{-12}\text{ F/m}$.

### 3. Geometry-Independent Permittivity Ratio

Since both the plexiglass and air measurements are performed at identical plate spacing $d = 2.8\text{ mm}$ and identical driving amplitude $V_m = 4\text{ V}$, the ratio of capacitances is independent of plate area and spacing:

```
C_plexi / C_air = m_plexi / m_air = K_plexi / K_air
```

### 4. Distance Dependence & Stray Capacitance

At constant frequency $f$ and voltage $V$, varying plate distance $d$ yields:

```
I(d) = (2π · f · V_rms · ε · A) · (1/d) + I_stray
```

The linear fit of $I$ versus $1/d$ yields an intercept $b_3$ attributable to fixture and cable stray capacitance:

```
C_stray = b_3 / (2π · f · V_rms)
```

### 5. Oscilloscope Calibration and Signal Measurements

The oscilloscope visualizes time-varying voltages:
- **Period & Frequency:** The horizontal division count $N_x$ for one cycle gives $T = N_x \times (\text{Time/Div})$, with $f = 1/T$.
- **Peak-to-Peak & RMS Voltage:** Vertical division count $N_y$ gives $V_{pp} = N_y \times (\text{Volts/Div})$, $V_m = V_{pp}/2$, and $V_{\mathrm{rms}} = V_m / \sqrt{2}$.
- **Current via Shunt Resistor:** Because oscilloscopes measure voltage, current is monitored across a small sensing resistor $R_s$:
  $$I(t) = \frac{V_s(t)}{R_s}$$
- **Lissajous Method (X-Y Mode):** Superposing two sinusoidal signals $x(t) = X_m \sin(\omega t)$ and $y(t) = Y_m \sin(\omega t + \Delta\phi)$ produces an ellipse. The phase shift satisfies:
  $$\sin(\Delta\phi) = \frac{y_0}{y_{\max}}$$
  where $y_0$ is the $y$-intercept when $x=0$.

---

## Experimental Setup

- **Circular Parallel-Plate Capacitor:** Two circular aluminum plates, radius $r = 10.0 \pm 0.1\text{ cm}$ ($A = 314.16 \pm 6.28\text{ cm}^2$), with a micrometric screw for plate distance adjustment.
- **Dielectric Insert:** Sheet of plexiglass, thickness $d = 2.80 \pm 0.05\text{ mm}$.
- **Function Generator:** Sinusoidal output, variable frequency ($1\text{ to }25\text{ kHz}$), regulated amplitude ($V_m = 4.0\text{ V}$, $10.0\text{ V}$).
- **Oscilloscope:** Dual-channel analog/digital oscilloscope for waveform inspection, time base calibration ($50\,\mu\text{s/div}$), and amplitude validation.
- **AC Microammeter:** Precision multirange AC microammeter for RMS current determination.

---

## Key Analyses

- **OLS Regression of $I$ vs. $f$ (Plexiglass):** Linear fit over 8 data points ($1\text{--}25\text{ kHz}$) yielding $C_{\mathrm{plexi}} = 238.89 \pm 0.73\text{ (stat)} \pm 3.07\text{ (total) pF}$ and $K_{\mathrm{plexi}} = 2.405 \pm 0.007\text{ (stat)} \pm 0.072\text{ (total)}$.
- **OLS Regression of $I$ vs. $f$ (Air):** Linear fit over 7 data points ($1\text{--}25\text{ kHz}$) yielding $C_{\mathrm{air}} = 58.60 \pm 0.72\text{ pF}$, $\varepsilon_{\mathrm{air}} = (5.223 \pm 0.065) \times 10^{-12}\text{ F/m}$, and $K_{\mathrm{air}} = 0.590 \pm 0.007$.
- **Geometry-Independent Ratio:** $C_{\mathrm{plexi}}/C_{\mathrm{air}} = 4.077 \pm 0.052$, proving a 4-fold capacitance enhancement by the dielectric.
- **Distance Dependence Regression ($I$ vs. $1/d$):** Confirms $I \propto 1/d$ ($R^2 = 0.9878$). The nonzero intercept $4.88 \pm 0.29\,\mu\text{A}$ isolates the fixture stray capacitance $C_{\mathrm{stray}} = 19.63 \pm 1.17\text{ pF}$.
- **Oscilloscope Frequency Verification:** Calibration at $10\text{ kHz}$ with $\text{Time/Div} = 50\,\mu\text{s}$ confirmed $T = 100\,\mu\text{s}$ ($N = 2.0\text{ div}$), with a reading limit of $\approx 5\%$.

---

## Results

| Parameter | Value | Statistical Error | Total Uncertainty | Unit |
| :--- | :---: | :---: | :---: | :---: |
| Plate Radius $r$ | $10.0$ | — | $\pm 0.1$ | $\mathrm{cm}$ |
| Plate Area $A$ | $314.16$ | — | $\pm 6.28$ | $\mathrm{cm^2}$ |
| Fixed Spacing $d$ | $2.80$ | — | $\pm 0.05$ | $\mathrm{mm}$ |
| Peak Driving Voltage $V_m$ | $4.00$ | — | $\pm 0.05$ | $\mathrm{V}$ |
| RMS Driving Voltage $V_{\mathrm{rms}}$ | $2.8284$ | — | $\pm 0.0354$ | $\mathrm{V}$ |
| **Slope $m_1$ (Plexiglass)** | **$4.2455$** | $\pm 0.0130$ | $\pm 0.0130$ | $\mu\mathrm{A/kHz}$ |
| $I$-Intercept $c_1$ (Plexiglass) | $0.7033$ | $\pm 0.1913$ | $\pm 0.1913$ | $\mu\mathrm{A}$ |
| Determination Coeff. $R_1^2$ | $0.99994$ | — | — | — |
| **Capacitance $C_{\mathrm{plexi}}$** | **$238.89$** | $\pm 0.73$ | $\pm 3.07$ | $\mathrm{pF}$ |
| **Dielectric Const. $K_{\mathrm{plexi}}$** | **$2.405$** | $\pm 0.007$ | $\pm 0.072$ | — |
| **Slope $m_2$ (Air)** | **$1.0414$** | $\pm 0.0129$ | $\pm 0.0129$ | $\mu\mathrm{A/kHz}$ |
| $I$-Intercept $c_2$ (Air) | $0.0143$ | $\pm 0.1967$ | $\pm 0.1967$ | $\mu\mathrm{A}$ |
| Determination Coeff. $R_2^2$ | $0.99924$ | — | — | — |
| **Capacitance $C_{\mathrm{air}}$** | **$58.60$** | $\pm 0.72$ | $\pm 1.03$ | $\mathrm{pF}$ |
| **Permittivity $\varepsilon_{\mathrm{air}}$** | **$5.223 \times 10^{-12}$** | $\pm 0.065 \times 10^{-12}$ | $\pm 0.169 \times 10^{-12}$ | $\mathrm{F/m}$ |
| **Relative Permittivity $K_{\mathrm{air}}$**| **$0.590$** | $\pm 0.007$ | $\pm 0.019$ | — |
| **Capacitance Ratio $C_{\mathrm{plexi}}/C_{\mathrm{air}}$** | **$4.077$** | $\pm 0.052$ | $\pm 0.052$ | — |
| **Distance Slope $a_3$ ($I$ vs $1/d$)** | **$28.92$** | $\pm 1.44$ | $\pm 1.44$ | $\mu\mathrm{A\cdot mm}$ |
| Distance Intercept $b_3$ | $4.88$ | $\pm 0.29$ | $\pm 0.29$ | $\mu\mathrm{A}$ |
| Determination Coeff. $R_3^2$ | $0.98784$ | — | — | — |
| **Stray Capacitance $C_{\mathrm{stray}}$** | **$19.63$** | $\pm 1.17$ | $\pm 1.19$ | $\mathrm{pF}$ |

---

## Files

- `7.tex` — Primary Persian LaTeX report source code
- `report7_final.tex` — Polished and verified Persian LaTeX report
- `analysis.py` — Python analysis script (OLS regressions, error propagation, plot generation)
- `plots/` — Generated publication-quality figures:
  - `I_vs_f_plexi.pdf` — Current vs. Frequency for Plexiglass with linear fit & error bars
  - `I_vs_f_air.pdf` — Current vs. Frequency for Air with linear fit & error bars
  - `I_vs_inv_d.pdf` — Current vs. $1/d$ confirming distance scaling and stray offset
  - `oscilloscope_waveform.pdf` — Simulated oscilloscope graticule showing period and $V_{pp}$
  - `lissajous_figures.pdf` — Simulated Lissajous figures for phase shift measurement
- `README.md` — This documentation file

---

## How to Compile

To compile the LaTeX report:

```bash
# Requires XeLaTeX due to xepersian and font requirements
xelatex report7_final.tex
xelatex report7_final.tex
```

To run the data analysis and regenerate all plots:

```bash
python analysis.py
```

---

## Notes & Error Sources

1. **Edge Effects & Fringing Fields:** The finite diameter of the circular plates ($2r = 20\text{ cm}$) produces fringe electric fields at the perimeter, causing the field lines to bow outwards. For small gaps ($d/r \ll 1$), this adds a small positive edge capacitance.
2. **Effective Permittivity of Air ($K_{\mathrm{air}} \approx 0.59$):** The apparent reduction in air permittivity is primarily due to systematic micro-screw zero calibration error, series cable resistance, meter loading, and voltage division across internal source impedance.
3. **Stray & Fixture Capacitance:** Connecting coaxial leads and terminals possess parallel parasitic capacitance ($C_{\mathrm{stray}} \approx 19.6\text{ pF}$), beautifully isolated by the intercept of the $I$ vs. $1/d$ plot.
4. **Oscilloscope vs. Meter Current Reading:** The oscilloscope input has high impedance ($1\text{ M}\Omega \parallel 20\text{ pF}$), whereas the AC microammeter inserts a small shunt burden. To measure currents on an oscilloscope, a calibrated shunt resistor must be placed in series.
