# Experiment 2: Kirchhoff's Laws and Wheatstone Bridge

This directory contains experimental data, numerical analysis scripts, vector plots, and LaTeX reports for the **Kirchhoff's Laws and Wheatstone Bridge** experiment in Physics Laboratory II at Sharif University of Technology.

---

## Overview

The purpose of this experiment is to experimentally investigate and verify fundamental circuit laws in direct current (DC) networks:
1. **Ohm's Law:** Verifying the linear constitutive relation $V = I R$ on individual resistors and across all network branches simultaneously.
2. **Kirchhoff's Current Law (KCL):** Verifying algebraic charge conservation at multi-branch circuit nodes:
   $$\sum I_{\text{in}} = \sum I_{\text{out}} \iff \sum_{k} I_k = 0$$
3. **Kirchhoff's Voltage Law (KVL):** Verifying energy conservation and path independence of electric potential around multiple closed loops:
   $$\sum \Delta V = 0 \iff \sum V_{\text{source}} = \sum V_{\text{drop}}$$
4. **Network Theory & Nodal Simulation:** Solving the full two-source, five-resistor bridge network analytically using nodal admittance formulation and comparing simulated currents and voltages with experimental measurements.
5. **Unknown Resistance Determination:** Measuring an unknown resistance $R_x$ (rheostat) via two independent methods:
   - Direct method: Simultaneous voltage and current measurement using Ohm's Law ($R_{x,1} = V_{R_x} / I_{R_x}$).
   - Null method: Balanced Wheatstone bridge condition ($R_{x,2} = \frac{R_1 R_3}{R_2}$).
   - Evaluating accuracy, sensitivity, and instrument loading effects between the two methods.

---

## Theoretical Background

### 1. Kirchhoff's Current Law (KCL)

Kirchhoff's First Law is a direct manifestation of the **continuity equation** and **conservation of electric charge**:

```
∇ · J + ∂ρ/∂t = 0
```

In steady-state DC conditions ($\partial\rho/\partial t = 0$), the divergence of current density vanishes ($\nabla \cdot \mathbf{J} = 0$). Integrating over a Gaussian surface enclosing an electrical node (junction) gives:

```
∮_S J · dA = 0  →  ∑_k I_k = 0  →  ∑ I_in = ∑ I_out
```

No net charge accumulates at any circuit junction.

---

### 2. Kirchhoff's Voltage Law (KVL)

Kirchhoff's Second Law is a consequence of the **conservative nature of static electric fields** (Faraday's law under magnetostatic conditions, $\partial \mathbf{B}/\partial t = 0$):

```
∇ × E = 0  →  ∮_C E · dl = 0
```

Since the electric field is the negative gradient of electrostatic potential ($\mathbf{E} = -\nabla V$), the line integral of $\mathbf{E}$ along any closed loop $C$ is identically zero:

```
∮_C (-∇V) · dl = 0  →  ∑_k ΔV_k = 0
```

Separating active sources (EMFs) and passive resistive drops:

```
∑ V_sources = ∑ I_k R_k
```

---

### 3. Balanced Wheatstone Bridge Principle

A Wheatstone bridge consists of four resistive arms arranged in a diamond loop ($R_1, R_2, R_3, R_x$) energized by an external voltage source $V_1$. A null detector (galvanometer) connects the two central nodes $B$ and $D$.

```
         A
        / \
       /   \
     R_1   R_3
     /       \
    /    G    \
   B --------- D
    \         /
     \       /
     R_2   R_x
       \   /
        \ /
         C
```

When the bridge is balanced, no current flows through the detector ($I_g = 0$), implying that nodes $B$ and $D$ are equipotential:

```
V_B = V_D
```

Because $I_g = 0$, the branches form simple voltage dividers:
- Current in branch $A \to B \to C$ is $I_L$:
  ```
  V_A - V_B = I_L R_1,    V_B - V_C = I_L R_2
  ```
- Current in branch $A \to D \to C$ is $I_R$:
  ```
  V_A - V_D = I_R R_3,    V_D - V_C = I_R R_x
  ```

Equating potential drops:

```
V_A - V_B = V_A - V_D  →  I_L R_1 = I_R R_3
V_B - V_C = V_D - V_C  →  I_L R_2 = I_R R_x
```

Dividing the two equations eliminates branch currents:

```
(I_L R_1) / (I_L R_2) = (I_R R_3) / (I_R R_x)  →  R_1 / R_2 = R_3 / R_x
```

Solving for the unknown resistance:

```
R_x = (R_1 · R_3) / R_2
```

**Key Advantage of Null Measurement:** Because the galvanometer draws zero current at balance, the measurement is completely immune to the detector's internal resistance, lead resistances of the galvanometer, and voltmeter loading errors.

---

### 4. Error Propagation Formulas

For all derived experimental quantities $f(x_1, x_2, \dots, x_n)$, standard first-order quadrature error propagation is applied:

```
δf = √[ ∑_i (∂f/∂x_i · δx_i)² ]
```

Specifically:
- **Theoretical Current:** $I = V / R$
  ```
  δI = I · √[ (δV / V)² + (δR / R)² ]
  ```
- **Theoretical Voltage:** $V = I · R$
  ```
  δV = V · √[ (δI / I)² + (δR / R)² ]
  ```
- **Direct Resistance:** $R_{x,1} = V_{R_x} / I_{R_x}$
  ```
  δR_{x,1} = R_{x,1} · √[ (δV_{R_x} / V_{R_x})² + (δI_{R_x} / I_{R_x})² ]
  ```
- **Wheatstone Resistance:** $R_{x,2} = (R_1 R_3) / R_2$
  ```
  δR_{x,2} = R_{x,2} · √[ (δR_1 / R_1)² + (δR_2 / R_2)² + (δR_3 / R_3)² ]
  ```
- **Discrepancy and Compatibility Index (z-score):**
  ```
  ΔR = |R_{x,1} - R_{x,2}|,    δ(ΔR) = √(δR_{x,1}² + δR_{x,2}²)
  z = ΔR / δ(ΔR)
  ```

---

## Experimental Setup

- **Fixed Resistors (5 elements):**
  - $R_1 = 388 \pm 1\,\Omega$
  - $R_2 = 219 \pm 1\,\Omega$
  - $R_3 = 47 \pm 1\,\Omega$
  - $R_4 = 46 \pm 1\,\Omega$
  - $R_5 = 99 \pm 1\,\Omega$
- **Variable Resistor / Rheostat:** Used as the unknown element $R_x$.
- **DC Power Supplies:**
  - Channel 1: $V_1 = 5.00 \pm 0.01\,\text{V}$
  - Channel 2: $V_2 = 8.00 \pm 0.01\,\text{V}$
- **Digital Multimeter (DMM):** Used as ohmmeter, DC milliammeter ($\pm 0.01\,\text{mA}$ resolution), and DC voltmeter ($\pm 0.01\,\text{V}$ resolution).
- **Galvanometer:** Zero-center microammeter for bridge null balance detection.
- **Circuit Breadboard & Connecting Leads.**

---

## Key Analyses

1. **Component-Level Ohm's Law Verification:**
   - Evaluated theoretical current $I_{\text{theo}} = V_{\text{exp}} / R$ and theoretical voltage $V_{\text{theo}} = I_{\text{exp}} R$ for all 5 resistors.
   - Performed global Ordinary Least Squares (OLS) regression of measured potential drop $V_{\text{exp}}$ against the calculated product $I_{\text{exp}} R$:
     $$\text{Slope} = 0.9860 \pm 0.0027, \quad \text{Intercept} = 0.0460 \pm 0.0115\,\text{V}, \quad R^2 = 0.999977$$
   - Confirmed Ohm's law linearity across all branches with maximum error under $0.89\%$.

2. **KCL Node Verification:**
   - **Node N:** $\sum I_{\text{in}} = I_1 + I_5 = 49.85 \pm 0.014\,\text{mA}$, $\sum I_{\text{out}} = I_4 = 49.57 \pm 0.010\,\text{mA}$.
     Residual $\Delta I = 0.28 \pm 0.017\,\text{mA}$ (relative error: $0.56\%$).
   - **Node M:** $\sum I_{\text{in}} = I_2 + I_4 = 59.73 \pm 0.014\,\text{mA}$, $\sum I_{\text{out}} = I_3 = 59.47 \pm 0.010\,\text{mA}$.
     Residual $\Delta I = 0.26 \pm 0.017\,\text{mA}$ (relative error: $0.44\%$).

3. **KVL Loop Verification:**
   - **Loop 1 ($V_2$ Right):** $+V_2 - V_{R4} + V_{R2} - V_{R1} = +0.000 \pm 0.020\,\text{V}$ (error: $0.00\%$).
   - **Loop 2 ($V_2$ Left):** $+V_2 - V_{R4} - V_{R3} - V_{R5} = +0.010 \pm 0.020\,\text{V}$ (error: $0.125\%$).
   - **Loop 3 ($V_1$ Top):** $+V_1 - V_{R2} - V_{R3} = -0.060 \pm 0.017\,\text{V}$ (error: $1.20\%$).
   - **Loop 4 ($V_1$ Bottom):** $+V_1 - V_{R1} + V_{R5} = -0.070 \pm 0.017\,\text{V}$ (error: $1.40\%$).

4. **Nodal Circuit Simulation:**
   - Solved the linear nodal admittance matrix for node potentials: $V_N = -2.9053\,\text{V}$, $V_M = +2.8075\,\text{V}$.
   - Simulated branch currents match measurements within $\pm 0.32\,\text{mA}$ ($< 1.6\%$).

5. **Unknown Resistor Determination ($R_x$):**
   - Direct method: $R_{x,1} = 81.06 \pm 0.12\,\Omega$.
   - Wheatstone bridge null method: $R_{x,2} = 83.27 \pm 1.82\,\Omega$.
   - Discrepancy: $\Delta R = 2.21 \pm 1.83\,\Omega$ ($2.65\%$), compatibility index $z = 1.21$ (consistent within $1.2\sigma$).

---

## Results

### Component Measurements and Ohm's Law Verification

| Resistor | $R\;[\Omega]$ | $V_{\text{exp}}\;[\text{V}]$ | $I_{\text{exp}}\;[\text{mA}]$ | $I_{\text{theo}}\;[\text{mA}]$ | Rel. Err ($I$) | $V_{\text{theo}}\;[\text{V}]$ | Rel. Err ($V$) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $R_1$ | $388 \pm 1$ | $7.96 \pm 0.01$ | $20.69 \pm 0.01$ | $20.515 \pm 0.059$ | $0.851\%$ | $8.028 \pm 0.021$ | $0.844\%$ |
| $R_2$ | $219 \pm 1$ | $2.24 \pm 0.01$ | $10.16 \pm 0.01$ | $10.228 \pm 0.065$ | $0.668\%$ | $2.225 \pm 0.010$ | $0.672\%$ |
| $R_3$ | $47 \pm 1$  | $2.82 \pm 0.01$ | $59.47 \pm 0.01$ | $60.000 \pm 1.294$ | $0.883\%$ | $2.795 \pm 0.059$ | $0.891\%$ |
| $R_4$ | $46 \pm 1$  | $2.28 \pm 0.01$ | $49.57 \pm 0.01$ | $49.565 \pm 1.099$ | $0.010\%$ | $2.280 \pm 0.050$ | $0.010\%$ |
| $R_5$ | $99 \pm 1$  | $2.89 \pm 0.01$ | $29.16 \pm 0.01$ | $29.192 \pm 0.312$ | $0.109\%$ | $2.887 \pm 0.029$ | $0.109\%$ |

### Kirchhoff's Current Law (KCL) Summary

| Node | Inflowing Branches | $\sum I_{\text{in}}\;[\text{mA}]$ | Outflowing Branches | $\sum I_{\text{out}}\;[\text{mA}]$ | $\Delta I\;[\text{mA}]$ | Rel. Discrepancy |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Node N** | $I_1 + I_5$ | $49.85 \pm 0.014$ | $I_4$ | $49.57 \pm 0.010$ | $0.28 \pm 0.017$ | $0.56\%$ |
| **Node M** | $I_2 + I_4$ | $59.73 \pm 0.014$ | $I_3$ | $59.47 \pm 0.010$ | $0.26 \pm 0.017$ | $0.44\%$ |

*Note: In the student LaTeX draft, an arithmetic transcription error ($20.69 + 29.16 = 49.16\,\text{mA}$) was recorded for Node N. Python re-computation yields the exact sum $49.85\,\text{mA}$ and residual $0.28\,\text{mA}$.*

### Kirchhoff's Voltage Law (KVL) Summary

| Loop | Loop Equation | Algebraic Sum $\sum \Delta V\;[\text{V}]$ | Loop EMF Reference | Rel. Error |
|:---|:---|:---:|:---:|:---:|
| **Loop 1** | $+V_2 - V_{R4} + V_{R2} - V_{R1}$ | $+0.000 \pm 0.020$ | $V_2 = 8.00\,\text{V}$ | $0.00\%$ |
| **Loop 2** | $+V_2 - V_{R4} - V_{R3} - V_{R5}$ | $+0.010 \pm 0.020$ | $V_2 = 8.00\,\text{V}$ | $0.125\%$ |
| **Loop 3** | $+V_1 - V_{R2} - V_{R3}$ | $-0.060 \pm 0.017$ | $V_1 = 5.00\,\text{V}$ | $1.20\%$ |
| **Loop 4** | $+V_1 - V_{R1} + V_{R5}$ | $-0.070 \pm 0.017$ | $V_1 = 5.00\,\text{V}$ | $1.40\%$ |

### Unknown Resistance Comparison

| Method | Governing Formula | Value $[\Omega]$ | Uncertainty $[\Omega]$ | Rel. Difference | Compatibility ($z$) |
|:---|:---|:---:|:---:|:---:|:---:|
| **Direct Method** | $R_{x,1} = V_{R_x} / I_{R_x}$ | $81.06$ | $\pm 0.12$ | — | — |
| **Wheatstone Bridge** | $R_{x,2} = (R_1 R_3) / R_2$ | $83.27$ | $\pm 1.82$ | $2.65\%$ | $1.21\sigma$ |

---

## Files

- `2.tex` — Complete Persian LaTeX source code of the lab report.
- `2.pdf` — Compiled report PDF.
- `analysis.py` — Python analysis script implementing OLS regressions, error propagation, network simulation, and plot generation.
- `plots/` — Directory containing generated publication-quality vector plots:
  - `ohms_law_verification.pdf`: Dual-panel plot of $V$ vs $I \cdot R$ with linear fit, ideal line, and mV residuals.
  - `kcl_node_balance.pdf`: Grouped bar chart comparing entering vs leaving currents at circuit nodes with error bars.
  - `kvl_loop_balance.pdf`: Residual voltage sums across independent loops with $1\sigma$ uncertainty intervals.
  - `wheatstone_bridge_comparison.pdf`: Bar chart comparing direct Ohm's law vs balanced Wheatstone bridge resistance with compatibility zone.
  - `circuit_simulation_comparison.pdf`: Comparison of theoretical nodal analysis vs experimental currents and voltages.
- `fig1.png`, `fig2.png`, `fig3.png` — Circuit schematics (general bridge, two-source network, and Wheatstone bridge configuration).
- `README.md` — This documentation file.

---

## How to Compile

To compile the LaTeX report, use `xelatex` (required for `xepersian` and Persian font typesetting):

```bash
xelatex 2.tex
xelatex 2.tex
```

To run the full numerical analysis, propagate errors, and regenerate all figures:

```bash
python analysis.py
```

---

## Notes & Error Discussion

1. **Ammeter Insertion Loading (Burden Voltage):**
   When measuring branch currents, the ammeter introduces a small non-zero internal resistance in series with the branch. For low-resistance branches (such as $R_3 = 47\,\Omega$ and $R_4 = 46\,\Omega$), an ammeter resistance of $0.5$ to $1\,\Omega$ contributes $1$–$2\%$ systematic attenuation to the measured current.

2. **Voltmeter Loading:**
   Although modern digital multimeters have an input impedance of $10\text{ M}\Omega$, which is very large compared to the hundred-ohm network resistances here, small contact resistances at breadboard clips can introduce millivolt-level potential drops.

3. **Joule Heating & Resistance Drift:**
   Branch $R_3$ carries approximately $59.5\,\text{mA}$, dissipating:
   $$P = I^2 R = (0.0595)^2 \times 47 \approx 0.166\,\text{W}$$
   For standard $1/4\,\text{W}$ carbon-film resistors, operating at $> 60\%$ rated power causes measurable temperature rise and resistance drift due to the positive/negative temperature coefficient of resistance (TCR).

4. **Why Wheatstone Bridge is Superior to the Direct Method:**
   The direct measurement of $R_x$ using separate voltmeter and ammeter probes is vulnerable to voltmeter burden current or ammeter series voltage drop. In contrast, the Wheatstone bridge relies on a **null condition** ($I_g = 0$), so the measurement is inherently immune to detector resistance, lead resistances, and source fluctuations.
