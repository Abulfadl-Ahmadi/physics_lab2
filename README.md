# Physics Lab II — Electricity & Magnetism Laboratory

**Sharif University of Technology**  
**Department of Physics**  
**Author:** Abolfazl Ahmadi  
**Academic Year:** Fall 2025 (1404)  
**Curriculum:** Electricity, Magnetism & AC Circuit Fundamentals  

---

## Overview

This repository contains the complete experimental data, statistical analysis scripts, publication-grade vector plots, and typeset academic reports for the **General Physics Lab II (Electricity & Magnetism)** course at Sharif University of Technology.

Each experiment directory represents an end-to-end scientific study: raw bench measurements, Python-based Ordinary Least Squares (OLS) regression, Gaussian quadrature error propagation, vector graphics in PDF format, and modular Obsidian-compatible study guides tailored for oral exam mastery.

### Physics Covered

The experiments trace the progression of classical electrodynamics, DC network analysis, and AC steady-state circuit physics:

1. **Ohm's Law & Conductor Resistivity (`1`)** — Linear $V-I$ characteristics, microscopic Drude conduction, geometric scaling ($R \propto l/S$), and resistivity determination for Nichrome, Galvanized Steel, and Chromium.
2. **Kirchhoff's Circuit Laws & Wheatstone Bridge (`2`)** — Verification of KCL nodal current conservation, KVL loop potential sums, nodal admittance matrix solving, and null-method bridge resistance metrology.
3. **Earth's Magnetic Field & Solenoid Distribution (`3`)** — Axial Helmholtz coil field profile, vacuum permeability $\mu_0$, and absolute geomagnetic field components ($B_E^h, B_E^v, |B_E|$) via tangent galvanometer.
4. **Transient RC & RL Circuits (`4`)** — Exponential capacitor charging and discharging, time constant $\tau = RC$, internal resistance of digital voltmeters, and series/parallel equivalent capacitances.
5. **Magnetic Force & Current Balance (`5`)** — Lorentz force on current-carrying conductors ($F = I L B$), current balance calibration, and field homogeneity across pole pieces.
6. **Alternating Current Circuits (`6`)** — Steady-state $50\text{ Hz}$ AC analysis in RL, RC, and RLC series networks, phasor diagrams, inductive/capacitive phase shifts, coil internal resistance ($r_L$), and resonance frequency.
7. **Oscilloscope & Dielectric Permittivity (`7`)** — Oscilloscope waveform metrology, parallel plate capacitor impedance ($I$ vs. $f$ and $I$ vs. $1/d$), relative dielectric constant of Plexiglass ($K \approx 2.41$), and stray fixture capacitance.
8. **AC Circuits, Resonance & Lissajous Figures (`8`)** — Oscilloscope X-Y phase shift determination via Lissajous ellipses, $2:1$ frequency figure-8 patterns, non-linear phase response $\tan\phi$, and series RLC resonance.

---

## Repository Structure

```
Lab-Phy-II/
├── README.md                          # Root repository documentation
├── .gitignore                         # LaTeX & Python build artifact exclusion
├── StudyGuides/                       # Obsidian-compatible Persian oral exam guides
│   ├── README.md                      # Guide index and oral exam overview
│   ├── 1_Ohm_Law.md
│   ├── 2_Kirchhoff.md
│   ├── 3_Magnetic_Field.md
│   ├── 4_RC_RL_Circuits.md
│   ├── 5_Magnetic_Force.md
│   ├── 5_RLC_Resonance.md
│   ├── 6_AC_Circuits.md
│   ├── 7_Oscilloscope.md
│   └── 8_Experiment8.md
└── Lab-Phy-2/                         # Experiment subdirectories
    ├── 1/ … 8/
    │   ├── analysis.py                # Standalone OLS regression & error propagation
    │   ├── README.md                  # Detailed per-experiment documentation
    │   ├── plots/                     # High-resolution vector figures (PDF)
    │   └── *.tex                      # Academic LaTeX report source
```

---

## Key Experimental Results Summary

| Experiment | Key Derived Parameter | Experimental Result | Reference / Nominal | Method / Notes |
|:---|:---|:---:|:---:|:---|
| **01. Ohm's Law** | NiCr Resistivity $\rho$<br>Galvanized Resistivity $\rho$ | $(1.170 \pm 0.057) \times 10^{-6}\,\Omega\cdot\text{m}$<br>$(0.243 \pm 0.016) \times 10^{-6}\,\Omega\cdot\text{m}$ | $1.10 - 1.50 \times 10^{-6}$<br>$0.15 - 0.25 \times 10^{-6}$ | Weighted mean across length & area fits<br>Ohmic linearity $R^2 > 0.99998$ |
| **02. Kirchhoff's Laws** | KCL Current Imbalance $\Delta I$<br>Bridge Resistance $R_x$ | $0.26 \pm 0.017\text{ mA}$ ($<0.5\%$ error)<br>$83.27 \pm 1.82\ \Omega$ | $0.00\text{ mA}$<br>$81.06 \pm 0.12\ \Omega$ (Direct) | Verified across 2 nodes & 4 loops<br>Wheatstone null method ($z=1.21\sigma$) |
| **03. Magnetic Field** | Vacuum Permeability $\mu_0$<br>Horizontal Earth Field $B_E^h$ | $(1.545 \pm 0.022) \times 10^{-6}\text{ H/m}$<br>$24.77 \pm 0.94\ \mu\text{T}$ | $1.257 \times 10^{-6}\text{ H/m}$<br>$\approx 25.0\ \mu\text{T}$ (Tehran) | Finite coil thickness correction<br>Tangent galvanometer linear fit |
| **04. RC & RL Circuits** | Voltmeter Resistance $\bar{R}_v$<br>Series Capacitance $C_{\mathrm{eq}}$ | $10.52 \pm 0.38\text{ M}\Omega$<br>$3.675 \pm 0.134\ \mu\text{F}$ | $10.0\text{ M}\Omega$ (Nominal)<br>$3.333\ \mu\text{F}$ | Semilogarithmic slope ($R^2 > 0.9995$)<br>Agrees within $10\%$ tolerance |
| **05. Magnetic Force** | Magnetic Induction $B(2\text{ A})$<br>Force-Length Slope | $11.75 \pm 0.31\text{ mT}$<br>$0.0475 \pm 0.0011\text{ N/m}$ | —<br>— | Current balance weighted mean<br>Length variation fit ($R^2 = 0.9989$) |
| **06. AC Circuits** | Coil Inductance $L_{\mathrm{RL}}$<br>Capacitance $C_{\mathrm{RC}}$ | $1.176 \pm 0.004\text{ H}$<br>$20.19 \pm 0.06\ \mu\text{F}$ | —<br>$20.0\ \mu\text{F}$ (Nominal) | $50\text{ Hz}$ phasor decomposition<br>Coil loss angle $\delta = 83.89^\circ$ |
| **07. Oscilloscope & Dielectric** | Plexiglass Constant $K_{\mathrm{plexi}}$<br>Stray Capacitance $C_{\mathrm{stray}}$ | $2.405 \pm 0.072$<br>$19.63 \pm 1.17\text{ pF}$ | $2.40 - 3.40$<br>— | Linear current vs frequency fit<br>$I(1/d)$ intercept isolation |
| **08. AC Resonance & Lissajous** | Resonance Frequency $f_{\mathrm{res}}$<br>Inductance $L_{\mathrm{res}}$ | $48.8 \pm 0.5\text{ Hz}$<br>$1.064 \pm 0.022\text{ H}$ | —<br>$1.06\text{ H}$ (Nominal) | Zero-phase Lissajous line collapse<br>Resonant LC balance condition |

---

## How to Use This Repository

### Prerequisites
- Python 3.8+ with `numpy`, `scipy`, and `matplotlib`.
- A modern LaTeX distribution (TeX Live / MacTeX / MiKTeX) supporting XeLaTeX and Persian typography (`xepersian`, `fouriernc`, `booktabs`, `siunitx`).

### Running the Analysis
Navigate to any experiment directory inside `Lab-Phy-2/` and execute the Python analysis script:

```bash
cd Lab-Phy-2/1
python analysis.py
```

This will:
1. Parse the circuit measurements and instrument uncertainty logs.
2. Perform OLS regressions with statistical error bounds on slopes and intercepts.
3. Compute derived impedances, reactances, and material properties via error quadrature.
4. Export publication-quality vector figures to the local `plots/` directory.

### Compiling Reports
Reports can be compiled via XeLaTeX:

```bash
cd Lab-Phy-2/1
xelatex 1.tex
```

---

## Methodology & Rigorous Error Analysis

Uncertainty calculations throughout this repository adhere strictly to the BIPM *Guide to the Expression of Uncertainty in Measurement* (GUM):

1. **Ordinary Least Squares (OLS) Linear Regression:** Standard error of fit parameters:
   $$s_m = \frac{s_{y/x}}{\sqrt{\sum (x_i - \bar{x})^2}}, \qquad s_c = s_{y/x} \sqrt{\frac{1}{N} + \frac{\bar{x}^2}{\sum (x_i - \bar{x})^2}}$$
2. **First-Order Quadrature Propagation:** For any derived quantity $q(x_1, \dots, x_k)$:
   $$\delta q = \sqrt{\sum_{i=1}^k \left( \frac{\partial q}{\partial x_i} \delta x_i \right)^2}$$
3. **Instrumental Uncertainty Limits:** Digital multimeter precision ($\pm 0.5\% + 2\text{ digits}$), analog scale tolerances ($\pm \frac{1}{2}\text{ div}$), and phase uncertainty via Lissajous screen reading.

---

## Academic Integrity & Authorship

All reports, derivations, and code implementations in this repository reflect rigorous academic scholarship at Sharif University of Technology. Real-world laboratory challenges (internal coil losses $r_L$, instrument burden voltage, stray capacitance, finite coil geometry) are modeled and resolved transparently.
