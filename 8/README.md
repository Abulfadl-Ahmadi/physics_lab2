# Experiment 8: Phase Difference in AC Circuits via Lissajous Figures & Determination of L, C, and Resonance Frequency

This directory contains the experimental data, analysis scripts, generated plots, and LaTeX report for **Experiment 8: AC Circuits and Lissajous Figures** in Physics Laboratory II at Sharif University of Technology.

---

## Overview

The primary objective of this experiment is to explore the application of **Lissajous figures** on an oscilloscope in XY mode to analyze alternating current (AC) circuits. Specifically:
1. Determine an unknown signal frequency by comparing it against a reference frequency using tangency point ratios.
2. Measure the phase difference $\phi$ between voltage and current across reactive components.
3. Determine the self-inductance $L$ of an inductor in a series $RL$ circuit from the linear dependence of $\tan\phi$ on frequency $f$.
4. Determine the capacitance $C$ of a capacitor in a series $RC$ circuit from the linear dependence of $\tan\phi$ on inverse frequency $1/f$.
5. Investigate resonance in a series $RLC$ circuit, observe the phase zero-crossing ($\phi = 0$) where the Lissajous ellipse collapses to a straight line, and compare the resonance inductance $L_{\mathrm{res}}$ against the single-component measurements.

---

## Theoretical Background

### 1. Lissajous Ellipse and Phase Measurement
When two sinusoidal signals of identical angular frequency $\omega$ are fed into the oscilloscope in XY mode:
```math
x(t) = A \sin(\omega t), \qquad y(t) = B \sin(\omega t + \phi)
```
Eliminating time $t$ yields the equation of an ellipse:
```math
\left(\frac{x}{A}\right)^2 + \left(\frac{y}{B}\right)^2 - 2\left(\frac{x}{A}\right)\left(\frac{y}{B}\right)\cos\phi = \sin^2\phi
```
At $x = 0$, $\sin(\omega t) = 0$, giving $y(t) = \pm B \sin\phi$. Hence, by measuring the vertical intercept $y_0$ (at $x=0$) and the maximum vertical excursion $B$:
```math
|\sin\phi| = \frac{|y_0|}{B} \implies \tan\phi = \frac{\sin\phi}{\sqrt{1 - \sin^2\phi}}
```

### 2. Frequency Ratio from Tangent Intercepts
When two signals with different rational frequencies $f_x$ and $f_y$ are displayed in XY mode, a stationary closed curve forms. Counting the tangency points or axis crossings gives:
```math
\frac{f_y}{f_x} = \frac{N_x}{N_y} \implies f_y = f_x \left(\frac{N_x}{N_y}\right)
```

### 3. Series RL Circuit
The complex impedance is $Z = R + j\omega L$. The phase angle of the total impedance satisfies:
```math
\tan\phi = \frac{\omega L}{R} = \left(\frac{2\pi L}{R}\right) f
```
A plot of $\tan\phi$ versus $f$ is a straight line through the origin with slope $a = \frac{2\pi L}{R}$, yielding:
```math
L = \frac{a R}{2\pi}, \qquad \sigma_L = L \sqrt{\left(\frac{\sigma_a}{a}\right)^2 + \left(\frac{\sigma_R}{R}\right)^2}
```

### 4. Series RC Circuit
The complex impedance is $Z = R - \frac{j}{\omega C}$. The phase angle is negative (current leads voltage):
```math
\tan\phi = -\frac{1}{\omega R C} = -\left(\frac{1}{2\pi R C}\right)\frac{1}{f}
```
Plotting $\tan\phi$ versus $1/f$ yields a straight line through the origin with negative slope $m = -\frac{1}{2\pi R C}$, giving:
```math
C = -\frac{1}{2\pi R m}, \qquad \sigma_C = C \sqrt{\left(\frac{\sigma_m}{|m|}\right)^2 + \left(\frac{\sigma_R}{R}\right)^2}
```

### 5. Series RLC Circuit and Resonance
The total impedance is:
```math
Z = R + j\left(\omega L - \frac{1}{\omega C}\right) \implies \tan\phi = \frac{\omega L - \frac{1}{\omega C}}{R}
```
At the resonance frequency $f_{\mathrm{res}}$, the reactive components cancel ($\omega_{\mathrm{res}} L = \frac{1}{\omega_{\mathrm{res}} C}$), yielding $\Im(Z) = 0$ and $\phi = 0$:
```math
f_{\mathrm{res}} = \frac{1}{2\pi\sqrt{LC}} \implies L_{\mathrm{res}} = \frac{1}{(2\pi f_{\mathrm{res}})^2 C}
```
On the oscilloscope screen, resonance is observed when the Lissajous ellipse collapses into a single straight diagonal line with zero phase difference.

---

## Experimental Setup

- **Dual-Channel Oscilloscope**: Operating in XY display mode.
- **Function Generators**: Two calibrated AC signal generators providing variable frequency sinusoidal waveforms ($20\ \text{Hz} - 150\ \text{Hz}$).
- **Resistor Decade Box**: Set to a precision reference resistance $R = 300\ \Omega \pm 1\%$.
- **Capacitor**: Electrolytic/film capacitor with nominal rating $C_{\mathrm{nom}} = 10\ \mu\text{F}$.
- **Inductor**: Iron-core coil / laboratory inductor.
- **BNC Cables and Banana Leads**: Low-noise shielded patch cables.

---

## Key Analyses

1. **Unknown Frequency Determination**: Reference signal $f_x = 100\ \text{Hz}$, tangency ratio $N_x/N_y = 2/4$, yielding $f_y = 50.0 \pm 0.3\ \text{Hz}$.
2. **RL Impedance Regression**: Constrained least-squares regression through the origin of $\tan\phi$ vs. $f$ over $f \in [30, 120]\ \text{Hz}$ to extract self-inductance $L$.
3. **RC Impedance Regression**: Constrained least-squares regression through the origin of $\tan\phi$ vs. $1/f$ to determine capacitance $C$.
4. **Resonance Identification**: Swept excitation frequency around $50\ \text{Hz}$ until $\phi = 0$ is detected at $f_{\mathrm{res}} = 48.8 \pm 0.5\ \text{Hz}$.
5. **Cross-Validation & Error Propagation**:
   - Theoretical resonance predicted from independent $L_{\mathrm{RL}}$ and $C_{\mathrm{RC}}$ measurements: $f_0 = 56.8 \pm 2.4\ \text{Hz}$.
   - Discrepancy analysis ($8.0\ \text{Hz}$ difference) accounting for the inductor's internal winding resistance $R_L$ and capacitor equivalent series resistance (ESR).

---

## Results

| Parameter | Experimental Value | Reference / Theory | Unit | Notes |
|:---|:---|:---|:---:|:---|
| Reference Frequency $f_x$ | $100.0 \pm 0.5$ | — | $\text{Hz}$ | Calibrated generator CH1 |
| Tangency Ratio $N_x / N_y$ | $2 / 4$ | $1 / 2$ | — | 2-lobe figure 8 pattern |
| Unknown Frequency $f_y$ | $50.0 \pm 0.3$ | $50.0$ | $\text{Hz}$ | Measured via Lissajous |
| Reference Resistor $R$ | $300.0 \pm 3.0$ | $300$ | $\Omega$ | Decade resistance box |
| RL Origin Slope $a$ | $(1.438 \pm 0.117) \times 10^{-2}$ | — | $\text{s}$ | $R^2 = 0.8270$ |
| Inductance $L_{\mathrm{RL}}$ | $0.686 \pm 0.056$ | $(0.69 \pm 0.06)$ | $\text{H}$ | Phase slope method |
| RC Origin Slope $m$ | $-46.31 \pm 0.80$ | — | $\text{s}^{-1}$ | $R^2 = 0.9958$ |
| Capacitance $C_{\mathrm{RC}}$ | $11.46 \pm 0.23$ | $10.0$ (nominal) | $\mu\text{F}$ | $+14.6\%$ electrolytic tolerance |
| Resonance Frequency $f_{\mathrm{res}}$ | $48.8 \pm 0.5$ | — | $\text{Hz}$ | In-phase line ($\phi = 0$) |
| $L_{\mathrm{res}}$ (using $C_{\mathrm{nom}} = 10\ \mu\text{F}$) | $1.064 \pm 0.022$ | $(1.06 \pm 0.02)$ | $\text{H}$ | Standard report formula |
| $L_{\mathrm{res}}$ (using $C_{\mathrm{RC}} = 11.46\ \mu\text{F}$) | $0.929 \pm 0.027$ | — | $\text{H}$ | Corrected with measured $C$ |
| Predicted Resonance $f_0(L_{\mathrm{RL}}, C_{\mathrm{RC}})$ | $56.8 \pm 2.4$ | — | $\text{Hz}$ | Discrepancy due to $R_L$ & ESR |

---

## Files

- `8.tex` — Full Persian LaTeX report with complete theoretical derivations, data tables, circuit diagrams (`circuitikz`), and TikZ plots.
- `8.pdf` — Compiled LaTeX report.
- `analysis.py` — Python analysis pipeline performing OLS regressions, uncertainty propagation, table summaries, and vector PDF plotting.
- `plots/RL_tan_vs_f.pdf` — Vector plot of $\tan\phi$ vs. $f$ with error bars and fitted regression lines.
- `plots/RC_tan_vs_inv_f.pdf` — Vector plot of $\tan\phi$ vs. $1/f$ with error bars and fitted regression lines.
- `plots/RLC_tan_vs_f.pdf` — Vector plot showing the continuous phase transition across resonance in the $RLC$ circuit.
- `plots/Lissajous_simulation.pdf` — Synthetic Lissajous figures demonstrating the phase evolution from $0^\circ$ to $180^\circ$ and the $2:1$ frequency ratio.
- `README.md` — This experiment documentation file.

---

## How to Compile

To compile the Persian LaTeX report:
```bash
xelatex 8.tex
# Or compile twice to ensure accurate cross-references and page labels:
xelatex 8.tex
xelatex 8.tex
```

To re-run the Python data analysis and regenerate all publication figures:
```bash
python analysis.py
```

---

## Notes & Systematic Error Analysis

1. **High Phase Angle Sensitivity**:
   The derivative of $\tan\phi$ with respect to $s = \sin\phi$ is:
   ```math
   \frac{d(\tan\phi)}{ds} = \frac{1}{(1 - s^2)^{3/2}} = \frac{1}{\cos^3\phi}
   ```
   As $\phi \to 90^\circ$ ($s \to 1$), the denominator approaches zero, leading to drastic error amplification in $\tan\phi$ from millimeter-level scale uncertainties on the oscilloscope CRT.
2. **Inductor Internal Resistance ($R_L$)**:
   In reality, a wire-wound inductor possesses non-negligible ohmic resistance $R_L$. The total series resistance is $R_{\mathrm{tot}} = R + R_L$, making $\tan\phi = \frac{\omega L}{R + R_L}$. Ignoring $R_L$ underestimates the true value of $L$ derived from the RL phase slope.
3. **Capacitor Non-Ideality & ESR**:
   Electrolytic capacitors exhibit non-negligible dielectric dissipation and Equivalent Series Resistance (ESR) at power frequencies ($30 - 120\ \text{Hz}$), shifting the effective phase angle and contributing to the observed $14.6\%$ deviation above nominal capacitance.
