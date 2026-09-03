# Experiment 4: RC Circuits — Charging, Discharging, Time Constant, and Capacitance

This directory contains experimental data, analysis scripts, generated plots, and LaTeX reports for the **RC Circuits (Charging and Discharging of Capacitors)** experiment in Physics Laboratory II at Sharif University of Technology.

## Overview

The primary objective of this experiment is to investigate the transient response of resistor-capacitor ($RC$) circuits during charging and discharging cycles. By measuring the time evolution of voltage across a high-impedance DC voltmeter acting as the resistive load ($R_v \sim 10\text{ M}\Omega$), we:
1. Verify the exponential decay law $V(t) = V_0 e^{-t/\tau}$.
2. Extract the circuit time constant $\tau = R C$ using Ordinary Least Squares (OLS) regression on semilogarithmic coordinates.
3. Determine the internal resistance of the laboratory voltmeter $R_v$ through two independent measurements with calibrated capacitors ($C_1 = 20\,\mu\text{F}$ and $C_2 = 4\,\mu\text{F}$).
4. Experimentally measure the equivalent capacitance of series and parallel combinations and compare them against theoretical formulas:
   $$\frac{1}{C_{\text{series}}} = \frac{1}{C_1} + \frac{1}{C_2}, \qquad C_{\text{parallel}} = C_1 + C_2$$

---

## Theoretical Background

### 1. Charging of a Capacitor through a Voltmeter

Consider a DC voltage source $V_0$ connected in series with an uncharged capacitor $C$ and a voltmeter of internal resistance $R_v$. Applying Kirchhoff's Voltage Law (KVL):

```
V_0 - V_C(t) - V_R(t) = 0
```

Since the charging current is $I(t) = \frac{dq}{dt} = C \frac{dV_C}{dt}$, the differential equation governing the capacitor voltage $V_C(t)$ is:

```
R_v C \frac{dV_C}{dt} + V_C(t) = V_0
```

With the initial condition $V_C(0) = 0$, the solution is:

```
V_C(t) = V_0 \left(1 - e^{-t / \tau}\right), \qquad \tau = R_v C
```

The voltmeter measures the potential drop across its own internal resistance $V(t) = V_R(t) = V_0 - V_C(t)$:

```
V(t) = V_0 e^{-t / \tau}
```

Thus, while the capacitor charges asymptotically to $V_0$, the voltmeter reading decays exponentially from $V_0$ to zero with time constant $\tau$.

---

### 2. Discharging of a Capacitor through a Voltmeter

When a capacitor pre-charged to $V_0$ is disconnected from the power supply and placed across the voltmeter terminals, KVL gives:

```
V_C(t) + I(t) R_v = 0 \implies R_v C \frac{dV}{dt} + V(t) = 0
```

Integrating with $V(0) = V_0$ yields:

```
V(t) = V_0 e^{-t / \tau}, \qquad \tau = R_v C
```

At $t = \tau$, the voltage drops to $V(\tau) = V_0 / e \approx 0.3679\, V_0$.

---

### 3. Linearization and Ordinary Least Squares (OLS)

Taking the natural logarithm of the normalized voltage ratio:

```
\ln\left(\frac{V(t)}{V_0}\right) = -\frac{1}{\tau}\, t
```

In terms of decimal (base-10) logarithms:

```
\log_{10}\left(\frac{V(t)}{V_0}\right) = -\frac{1}{\tau \ln 10}\, t = b\, t + a
```

Defining $x_i = t_i$ and $y_i = \log_{10}(V_i / V_0)$, the slope $b$ and intercept $a$ from OLS regression are:

```
b = \frac{N \sum x_i y_i - \sum x_i \sum y_i}{N \sum x_i^2 - (\sum x_i)^2}, \qquad a = \frac{1}{N}\left(\sum y_i - b \sum x_i\right)
```

The time constant and internal resistance are extracted via:

```
\tau = -\frac{1}{b \ln 10} = -\frac{1}{m_{\ln}}, \qquad R_v = \frac{\tau}{C}
```

Using first-order error propagation:

```
\delta \tau = \tau \frac{\mathrm{SE}_m}{|m|}, \qquad \delta R_v = R_v \sqrt{\left(\frac{\delta \tau}{\tau}\right)^2 + \left(\frac{\delta C}{C}\right)^2}
```

---

### 4. Capacitor Combinations

- **Series Combination:**
  $$\frac{1}{C_{\text{eq}}^{(\mathrm{th})}} = \frac{1}{C_1} + \frac{1}{C_2} \implies C_{\text{eq}}^{(\mathrm{th})} = \frac{C_1 C_2}{C_1 + C_2} = \frac{20 \times 4}{20 + 4} = 3.333\,\mu\text{F}$$
  $$C_{\text{eq}}^{(\mathrm{exp})} = \frac{\tau_s}{\bar{R}_v}$$

- **Parallel Combination:**
  $$C_{\text{eq}}^{(\mathrm{th})} = C_1 + C_2 = 20 + 4 = 24.000\,\mu\text{F}$$
  $$C_{\text{eq}}^{(\mathrm{exp})} = \frac{\tau_p}{\bar{R}_v}$$

- **Relative Error:**
  $$\delta C = \frac{|C_{\text{eq}}^{(\mathrm{exp})} - C_{\text{eq}}^{(\mathrm{th})}|}{C_{\text{eq}}^{(\mathrm{th})}} \times 100\%$$

---

### 5. Theoretical Analogy with RL Circuits

In dual fashion, an inductor-resistor ($RL$) circuit governed by $\mathcal{E} - L \frac{dI}{dt} - IR = 0$ exhibits transient current growth and decay characterized by the inductive time constant:

```
\tau_L = \frac{L}{R}
```

During current build-up: $I(t) = \frac{\mathcal{E}}{R} (1 - e^{-t / \tau_L})$.  
During current decay: $I(t) = I_0 e^{-t / \tau_L}$.

---

## Experimental Setup

- **DC Power Supply:** Regulated DC supply set to $V_0 = 10.0\text{ V}$.
- **Electrolytic / Film Capacitors:**
  - $C_1 = 20\,\mu\text{F}$ (nominal)
  - $C_2 = 4\,\mu\text{F}$ (nominal)
- **High-Impedance Voltmeter / Multimeter:** Digital voltmeter with input impedance $R_v \approx 10\text{ M}\Omega$.
- **Digital Stopwatch:** Precision $\pm 0.5\text{ s}$.
- **Switching Key and Low-Resistance Leads.**

---

## Key Analyses

1. **Part 1 (Charging $C_1 = 20\,\mu\text{F}$):**
   - Normalized voltage $V/V_0$ recorded over $t = 0$ to $135\text{ s}$ ($\Delta t = 15\text{ s}$).
   - OLS fit on $\ln(V/V_0)$ vs. $t$: slope $m_1 = (-5.539 \pm 0.004) \times 10^{-3}\text{ s}^{-1}$, $R^2 = 0.999995$.
   - Time constant $\tau_1 = 180.53 \pm 0.14\text{ s}$.
   - Voltmeter resistance $R_{v1} = 9.026 \pm 0.007\text{ (fit)} \pm 0.451\text{ (tot) M}\Omega$.

2. **Part 2 (Discharging $C_2 = 4\,\mu\text{F}$):**
   - Capacitor pre-charged to $10.0\text{ V}$, discharged through voltmeter over $135\text{ s}$.
   - OLS fit: slope $m_2 = (-2.081 \pm 0.003) \times 10^{-2}\text{ s}^{-1}$, $R^2 = 0.999986$.
   - Time constant $\tau_2 = 48.06 \pm 0.06\text{ s}$.
   - Voltmeter resistance $R_{v2} = 12.014 \pm 0.016\text{ (fit)} \pm 0.601\text{ (tot) M}\Omega$.
   - Mean voltmeter resistance: $\bar{R}_v = 10.52 \pm 0.38\text{ M}\Omega$.

3. **Part 3 (Series Discharge $C_1$ & $C_2$):**
   - Discharge of series pair over $135\text{ s}$.
   - OLS fit: $\tau_s = 38.66 \pm 0.30\text{ s}$ ($R^2 = 0.999520$).
   - (Linear interpolation at $V_0/e = 3.68\text{ V}$ yields $\tau \approx 39.45\text{ s}$).
   - Equivalent capacitance: $C_{\text{eq}}^{(\mathrm{exp})} = 3.68 \pm 0.13\,\mu\text{F}$ vs. $C_{\text{th}} = 3.33\,\mu\text{F}$ (relative error: $10.25\%$).

4. **Part 4 (Parallel Discharge $C_1$ & $C_2$):**
   - Discharge of parallel pair over $300\text{ s}$ ($\Delta t = 30\text{ s}$).
   - OLS fit: $\tau_p = 228.60 \pm 0.27\text{ s}$ ($R^2 = 0.999987$).
   - (Linear interpolation at $V_0/e = 3.62\text{ V}$ yields $\tau \approx 230.03\text{ s}$).
   - Equivalent capacitance: $C_{\text{eq}}^{(\mathrm{exp})} = 21.73 \pm 0.78\,\mu\text{F}$ vs. $C_{\text{th}} = 24.00\,\mu\text{F}$ (relative error: $9.46\%$).

---

## Results

| Parameter | Theoretical | Experimental (OLS) | Uncertainty | Unit | Rel. Error |
|:---|:---:|:---:|:---:|:---:|:---:|
| $\tau_1$ ($C_1$ charging) | — | 180.53 | $\pm 0.14$ | s | — |
| $R_{v1}$ (from $C_1$) | — | 9.026 | $\pm 0.451$ | $\text{M}\Omega$ | — |
| $\tau_2$ ($C_2$ discharging) | — | 48.06 | $\pm 0.06$ | s | — |
| $R_{v2}$ (from $C_2$) | — | 12.014 | $\pm 0.601$ | $\text{M}\Omega$ | — |
| $\bar{R}_v$ (mean internal res.) | — | 10.520 | $\pm 0.376$ | $\text{M}\Omega$ | — |
| $\tau_{\text{series}}$ ($C_1 \| C_2$) | 35.07 | 38.66 | $\pm 0.30$ | s | $10.25\%$ |
| $C_{\text{eq, series}}$ | 3.333 | 3.675 | $\pm 0.134$ | $\mu\text{F}$ | $10.25\%$ |
| $\tau_{\text{parallel}}$ ($C_1 + C_2$) | 252.48 | 228.60 | $\pm 0.27$ | s | $9.46\%$ |
| $C_{\text{eq, parallel}}$ | 24.000 | 21.730 | $\pm 0.777$ | $\mu\text{F}$ | $9.46\%$ |

---

## Files

- `4.tex` — Full Persian LaTeX report with embedded TikZ/pgfplots figures.
- `4-1.tex` — Alternative LaTeX report draft.
- `4.pdf` — Compiled PDF report.
- `analysis.py` — Python analysis script containing OLS regressions, error propagation, and plot generation.
- `plots/` — Generated publication-quality PDF figures:
  - `part1_C1_charging.pdf`: Semilogarithmic decay of voltmeter reading and dual-curve charging trajectory.
  - `part2_C2_discharging.pdf`: Exponential and semilogarithmic discharge of $C_2$.
  - `part3_series_discharging.pdf`: Series configuration discharge curve with $V_0/e$ marker.
  - `part4_parallel_discharging.pdf`: Parallel configuration discharge curve with $V_0/e$ marker.
  - `all_circuits_log_comparison.pdf`: Unified $\log_{10}(V/V_0)$ vs. $t$ comparison for all 4 configurations.
  - `capacitance_comparison.pdf`: Bar chart comparing experimental and theoretical equivalent capacitances.
- `README.md` — This documentation file.

---

## How to Compile

Because the report is written in Persian and uses the `xepersian` package, compile using `xelatex`:

```bash
xelatex 4.tex
xelatex 4.tex
```

To run the Python analysis and regenerate all figures:

```bash
python analysis.py
```

---

## Notes & Discussion

1. **Why does the voltmeter reading decay during charging?**
   In the charging circuit, the voltmeter is placed in series with the capacitor. By KVL, the voltmeter reads the voltage across its own resistance: $V_R(t) = V_0 - V_C(t)$. As charge accumulates on the capacitor plates, $V_C(t) \to V_0$, reducing the series current $I(t) \to 0$ and causing the voltmeter reading to decay exponentially to zero.

2. **Difference between $R_{v1}$ ($9.0\text{ M}\Omega$) and $R_{v2}$ ($12.0\text{ M}\Omega$):**
   Digital multimeters have an internal resistive divider network. Small voltage-dependent input impedances, capacitor leakage currents (dielectric absorption in electrolytic capacitors), and the nominal tolerance of the capacitor rating ($\pm 5\%$ to $\pm 10\%$) account for the apparent difference between the two calculated voltmeter resistances.

3. **Comparison between OLS Regression and Interpolation at $V_0/e$:**
   The report uses a two-point linear interpolation to find the moment when $V(t) = V_0/e$ ($\tau \approx 39\text{ s}$ for series and $\tau \approx 230\text{ s}$ for parallel). In contrast, OLS regression utilizes all 10–11 data points, averaging out random timing jitter and voltmeter read errors. The resulting values ($\tau_s = 38.66\text{ s}$, $\tau_p = 228.60\text{ s}$) match the non-linear Levenberg-Marquardt fits ($38.66\text{ s}$ and $228.96\text{ s}$) with extreme precision ($R^2 > 0.9995$).
