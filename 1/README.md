# Experiment 1: Ohm's Law and Resistance Parameters

This directory contains the experimental data, analysis scripts, publication-grade vector plots, and the LaTeX report for the **Ohm's Law and Resistance Parameters** experiment at Sharif University of Technology.

## Overview

The primary objectives of this experiment are:
1. **Verifying Ohm's Law ($V = R \cdot I$)**: Investigating the current-voltage relationship of a cylindrical metallic conductor to confirm linearity and Ohmic behavior.
2. **Studying Length Dependence ($R = f(l)$)**: Demonstrating the direct proportionality between resistance and wire length at constant cross-section.
3. **Studying Cross-Sectional Area Dependence ($R = f(1/S)$)**: Verifying the inverse proportionality between resistance and wire cross-sectional area.
4. **Determining Material Resistivity ($\rho$)**: Measuring the electrical resistivity of three distinct engineering materials (Nickel-Chromium alloy, Galvanized steel, and Pure Chromium).

---

## Theoretical Background

### Ohm's Law and Linear Conduction

In a microscopic framework, when an electric field $\mathbf{E}$ is applied across a metallic conductor, free conduction electrons experience a drift velocity $\mathbf{v}_d = -\mu_e \mathbf{E}$. The current density is given by Ohm's law in microscopic form:

```
J = σ · E = (1 / ρ) · E
```

For a macroscopic cylindrical conductor of length $l$, cross-sectional area $S$, and uniform electric field $E = V / l$, the total electric current is:

```
I = J · S = (V / (ρ · l)) · S = V / R
```

Rearranging gives the macroscopic form of Ohm's law:

```
V = R · I
```

where $R$ is the electrical resistance of the wire. A conductor is defined as **Ohmic** if its resistance $R$ remains constant across varying applied currents and voltages, producing a straight line passing through the origin in the $V-I$ plane ($y$-intercept $V_0 = 0$).

### Geometric and Material Factors Influencing Resistance

The resistance of a uniform wire depends on its geometry and the intrinsic electronic properties of the material:

```
R = ρ · (l / S) = ρ · (l / (π · (d / 2)²)) = (4 · ρ · l) / (π · d²)
```

where:
- $\rho$ is the electrical resistivity ($\Omega \cdot \mathrm{m}$),
- $l$ is the length of the wire ($\mathrm{m}$),
- $S = \pi (d/2)^2$ is the circular cross-sectional area ($\mathrm{m}^2$),
- $d$ is the wire diameter ($\mathrm{m}$).

From this relation, three distinct linear scaling laws are experimentally tested:
1. **Length dependence**: $R = m_l \cdot l$, where the slope is $m_l = \frac{\rho}{S}$.
2. **Cross-sectional area dependence**: $R = m_S \cdot \left(\frac{1}{S}\right)$, where the slope is $m_S = \rho \cdot l$.
3. **Resistivity calculation**:
   ```
   ρ = R · S / l = (V · π · d²) / (4 · I · l)
   ```

---

## Experimental Setup

The setup consists of an experimental wire board equipped with five 1.00-meter calibrated wires, mounted alongside metric measurement scales:

- **Wire 1**: Nickel-Chromium (NiCr, Nichrome), diameter $d_1 = 0.25\ \mathrm{mm}$
- **Wire 2**: Nickel-Chromium (NiCr, Nichrome), diameter $d_2 = 0.40\ \mathrm{mm}$
- **Wire 3**: Nickel-Chromium (NiCr, Nichrome), diameter $d_3 = 0.30\ \mathrm{mm}$
- **Wire 4**: Galvanized iron/steel, diameter $d_4 = 0.30\ \mathrm{mm}$
- **Wire 5**: Pure Chromium (Cr), diameter $d_5 = 0.40\ \mathrm{mm}$
- **DC Regulated Power Supply**: Continuously variable $0 - 30\ \mathrm{V}$ output with fine current limiting.
- **Digital Multimeter (DMM)**:
  - Voltage mode: DC voltmeter, resolution $\delta V = 0.01\ \mathrm{V}$
  - Current mode: DC ammeter, resolution $\delta I = 1.0\ \mathrm{mA}$ ($0.001\ \mathrm{A}$)
- **Measurement Scale & Caliper**: Length resolution $\delta l = 1.0\ \mathrm{mm}$ ($0.001\ \mathrm{m}$), diameter tolerance $\delta d = 0.01\ \mathrm{mm}$ ($1 \times 10^{-5}\ \mathrm{m}$).

---

## Key Analyses

### 1. $V$ vs. $I$ Linear Regression (Part 1)
Current through Wire 1 is stepped from $100\ \mathrm{mA}$ to $500\ \mathrm{mA}$ while recording the potential drop across the wire.
- Ordinary Least Squares (OLS) regression: $V = R \cdot I + V_0$
- Extraction of $R = 26.620 \pm 0.059\ \Omega$ with coefficient of determination $R^2 = 0.999985$.
- The intercept $V_0 = -0.024 \pm 0.020\ \mathrm{V}$ is consistent with zero within $1.2\sigma$, confirming strictly Ohmic behavior without thermal runaway.

### 2. Resistance vs. Wire Length $R(l)$ (Part 2)
With Wire 1 energized at a fixed current $I = 250\ \mathrm{mA}$, voltage taps at $l \in \{10, 27, 50, 80, 100\}\ \mathrm{cm}$ yield:
- OLS regression: $R = m_l \cdot l + R_0$
- Measured slope: $m_l = 26.433 \pm 0.094\ \Omega/\mathrm{m}$ ($R^2 = 0.999962$).
- Zero offset: $R_0 = 0.045 \pm 0.059\ \Omega$ confirming zero lead/contact resistance.
- Resistivity from length slope:
  $$\rho_1 = m_l \cdot S_1 = (1.298 \pm 0.104) \times 10^{-6}\ \Omega\cdot\mathrm{m}$$

### 3. Resistance vs. Inverse Cross-Sectional Area $R(1/S)$ (Part 3)
Three identical-length ($L = 1.00\ \mathrm{m}$) Nichrome wires of different diameters ($0.25, 0.30, 0.40\ \mathrm{mm}$) are connected in series at $I = 250\ \mathrm{mA}$:
- OLS regression of $R$ vs. $1/S$:
  $$R = m_S \cdot (1/S) + c_S$$
  $$m_S = (1.447 \pm 0.212) \times 10^{-6}\ \Omega\cdot\mathrm{m}^2 \quad (R^2 = 0.9791)$$
- Resistivity from area slope:
  $$\rho_2 = \frac{m_S}{L} = (1.447 \pm 0.212) \times 10^{-6}\ \Omega\cdot\mathrm{m}$$

### 4. Material Resistivity Comparison (Part 4)
Wires of equal length ($L = 1.00\ \mathrm{m}$) made of Nichrome, Galvanized steel, and Pure Chromium are evaluated at $I = 250\ \mathrm{mA}$:
- Direct calculation with full partial-derivative error propagation:
  $$\left(\frac{\delta \rho}{\rho}\right) = \sqrt{\left(\frac{\delta V}{V}\right)^2 + \left(\frac{\delta I}{I}\right)^2 + \left(2\frac{\delta d}{d}\right)^2 + \left(\frac{\delta L}{L}\right)^2}$$

---

## Results

| Quantity | Experimental Value | Fit Quality ($R^2$) | Literature / Benchmark | Unit |
| :--- | :--- | :--- | :--- | :--- |
| **Wire 1 Resistance ($R_{\text{OLS}}$)** | $26.620 \pm 0.059$ | $0.999985$ | $26.48$ (pointwise mean) | $\Omega$ |
| **Voltage Intercept ($V_0$)** | $-0.024 \pm 0.020$ | — | $0.000$ (ideal Ohmic) | $\mathrm{V}$ |
| **Length Gradient ($dR/dl$)** | $26.433 \pm 0.094$ | $0.999962$ | — | $\Omega/\mathrm{m}$ |
| **Length Intercept ($R_0$)** | $0.045 \pm 0.059$ | — | $0.000$ (no contact res.) | $\Omega$ |
| **Area Gradient ($m_S = \rho L$)** | $(1.447 \pm 0.212) \times 10^{-6}$ | $0.979075$ | — | $\Omega\cdot\mathrm{m}^2$ |
| **$\rho$ (NiCr, Length Fit)** | $(1.298 \pm 0.104) \times 10^{-6}$ | — | $1.05 - 1.50 \times 10^{-6}$ | $\Omega\cdot\mathrm{m}$ |
| **$\rho$ (NiCr, Area Fit)** | $(1.447 \pm 0.212) \times 10^{-6}$ | — | $1.05 - 1.50 \times 10^{-6}$ | $\Omega\cdot\mathrm{m}$ |
| **$\rho$ (NiCr, Direct Wire 3)** | $(1.077 \pm 0.072) \times 10^{-6}$ | — | $1.05 - 1.50 \times 10^{-6}$ | $\Omega\cdot\mathrm{m}$ |
| **$\rho$ (NiCr, Weighted Mean)** | $\mathbf{(1.170 \pm 0.057) \times 10^{-6}}$ | — | **$1.10 - 1.50 \times 10^{-6}$** | $\Omega\cdot\mathrm{m}$ |
| **$\rho$ (Galvanized Steel)** | $(0.243 \pm 0.016) \times 10^{-6}$ | — | $0.15 - 0.25 \times 10^{-6}$ | $\Omega\cdot\mathrm{m}$ |
| **$\rho$ (Pure Chromium)** | $(1.508 \pm 0.076) \times 10^{-6}$ | — | $1.25 - 1.60 \times 10^{-6}$ | $\Omega\cdot\mathrm{m}$ |

---

## Files

- `1.tex` — LaTeX source of the full Persian laboratory report
- `1.pdf` — Compiled PDF document
- `analysis.py` — Complete Python script performing OLS regressions, uncertainty propagation, and generating publication plots
- `plots/` — Generated publication-quality vector PDF figures:
  - `V_vs_I_ohms_law.pdf`: Voltage vs. current with error bars, linear fit, and residual panel
  - `R_vs_length.pdf`: Resistance as a function of wire length $R(l)$
  - `R_vs_inv_area.pdf`: Resistance as a function of inverse cross-sectional area $R(1/S)$
  - `resistivity_comparison.pdf`: Bar chart comparing experimental resistivities across materials
  - `ohm_law_summary.pdf`: Four-panel comprehensive summary figure
- `README.md` — This technical documentation

---

## How to Compile

To run the complete data analysis and generate all figures:
```bash
python analysis.py
```

To compile the LaTeX report:
```bash
xelatex 1.tex
```
*(Note: Use `xelatex` due to the Persian typesetting package `xepersian` and Persian font requirements).*

---

## Notes

- **Linearity and Joule Heating**: In Part 1, the linear model holds up to $500\ \mathrm{mA}$ with maximum deviation below $1.6\%$. Joule heating ($P = I^2 R \approx 3.3\ \mathrm{W}$ at $500\ \mathrm{mA}$) causes slight temperature elevation, but the low temperature coefficient of resistivity for Nichrome ($\alpha \approx 4 \times 10^{-4}\ \mathrm{K}^{-1}$) ensures minimal thermal drift.
- **Zero Intercepts**: Both the $V-I$ intercept ($-0.024 \pm 0.020\ \mathrm{V}$) and the $R-l$ intercept ($0.045 \pm 0.059\ \Omega$) are statistically zero within $1-2$ standard errors, demonstrating high instrument calibration accuracy and negligible contact resistance.
- **Material Classification**: The measured resistivities clearly place all three tested specimens in the metallic conductor regime ($\sim 10^{-7}$ to $10^{-6}\ \Omega\cdot\mathrm{m}$), demonstrating the effect of alloying (NiCr) in drastically increasing resistivity compared to basic structural metals.
