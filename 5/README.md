# Experiment 5: Magnetic Force on a Current-Carrying Conductor (Current Balance)

This directory contains the experimental data, Python analysis scripts, plots, and the final LaTeX report for the **Magnetic Force on a Current-Carrying Conductor** experiment in Physics Laboratory II at Sharif University of Technology.

> [!NOTE]
> In some earlier syllabus listings or tracking sheets, Experiment 5 was indexed as *RLC Resonance*. In this laboratory manual and codebase, Experiment 5 is the **Magnetic Force / Current Balance** experiment, corresponding directly to `5.tex`.

---

## Overview

The purpose of this experiment is to study the magnetic force (Lorentz force) exerted on a straight current-carrying conductor placed inside an external magnetic field:
$$\vec{F} = i\,(\vec{L} \times \vec{B})$$

Using a precision current balance setup with an electromagnet, the experiment systematically investigates:
1. **Angular dependence:** Verifying $F \propto \sin\theta$ and finding maximal force when $\vec{L} \perp \vec{B}$ ($\theta = 90^\circ$).
2. **Length dependence:** Verifying $F \propto L$ at fixed currents ($i = 4.0\text{ A}$, $I_m = 2.0\text{ A}$) to extract the magnetic field $B$.
3. **Loop current dependence:** Verifying $F \propto i$ at fixed length ($L = 10\text{ cm}$) and coil current ($I_m = 2.0\text{ A}$) to extract $B$.
4. **Electromagnet current dependence:** Verifying $B \propto I_m$ at fixed loop current ($i = 4.0\text{ A}$) and length ($L = 10\text{ cm}$) to determine the coil calibration constant $\alpha = B/I_m$.

---

## Theoretical Background

### Lorentz Force on a Current-Carrying Wire

A charge carrier $q$ moving with drift velocity $\vec{v}_d$ in a magnetic field $\vec{B}$ experiences the Lorentz force:
```
F_L = q(v_d × B)
```
Integrating over all conduction electrons in a straight wire segment of length $L$ and cross section $A$ carrying steady current $i = n q v_d A$:
```
F = i (L × B)
```
If the conductor is oriented at an angle $\theta$ relative to $\vec{B}$, the scalar magnitude of the force is:
```
F = i L B sin(θ)
```
For perpendicular orientation ($\theta = 90^\circ$, $\sin\theta = 1$):
```
F = i L B
```

### Electromagnet Field Relation

For an iron-core electromagnet operating in the linear, unsaturated region of its magnetization curve, the magnetic induction $B$ between the pole pieces is directly proportional to the coil current $I_m$:
```
B = α · I_m
```
where $\alpha$ is the coil geometry and permeability factor (in $\text{T/A}$). Consequently, the force equation becomes:
```
F = α · i · L · I_m
```

### Linear Regressions and Sensitivity

1. **Varying length $L$ (fixed $i, I_m$):**
   $$\Delta F = a_L L + b_L, \quad a_L = i B \implies B = \frac{a_L}{i}$$
2. **Varying loop current $i$ (fixed $L, I_m$):**
   $$\Delta F = a_i i + b_i, \quad a_i = L B \implies B = \frac{a_i}{L}$$
3. **Varying magnet current $I_m$ (fixed $L, i$):**
   $$\Delta F = a_m I_m + b_m, \quad a_m = i L \alpha \implies \alpha = \frac{a_m}{i L}, \quad B(I_m) = \alpha I_m$$

---

## Experimental Setup

- **U-shaped Electromagnet:** Iron core with adjustable pole shoe gap (set to $\approx 1\text{ cm}$ or $4\text{ cm}$), energized by a regulated DC power supply ($I_m = 0\text{--}2.0\text{ A}$).
- **Conductor Loops:** Rigid U-shaped printed circuit boards with horizontal active segment lengths $L = 1.25, 2.50, 5.00, 10.00\text{ cm}$.
- **Current Balance Support:** Counterbalanced beam mounted on an analytical balance or force sensor measuring force deviations $\Delta F$ with a resolution of $\pm 0.05\text{ mN}$.
- **DC Current Source:** Regulated high-current supply feeding loop current $i = 0\text{--}4.0\text{ A}$.
- **Flexible Lead Wires:** Extremely flexible bare copper leads connected with a gentle catenary slack to prevent mechanical torque or spring forces on the balance beam.

---

## Key Analyses

- **OLS Fits:** Ordinary least-squares linear fits implemented in `analysis.py` with standard errors for slope, intercept, and determination coefficient $R^2$.
- **Error Propagation:** Complete partial-derivative quadrature for derived magnetic induction:
  $$\delta B_A = B_A \sqrt{\left(\frac{\delta a_L}{a_L}\right)^2 + \left(\frac{\delta i}{i}\right)^2}$$
  $$\delta B_B = B_B \sqrt{\left(\frac{\delta a_i}{a_i}\right)^2 + \left(\frac{\delta L}{L}\right)^2}$$
  $$\delta \alpha = \alpha \sqrt{\left(\frac{\delta a_m}{a_m}\right)^2 + \left(\frac{\delta i}{i}\right)^2 + \left(\frac{\delta L}{L}\right)^2}$$
- **Weighted Average:** Inverse-variance weighting across the three independent determinations of $B(I_m = 2.0\text{ A})$:
  $$w_k = \frac{1}{\delta B_k^2}, \quad \bar{B} = \frac{\sum w_k B_k}{\sum w_k}, \quad \delta \bar{B} = \frac{1}{\sqrt{\sum w_k}}$$

---

## Results

| Method / Independent Variable | Condition | Slope | Intercept | $R^2$ | Extracted $B$ or $\alpha$ |
|:---|:---|:---|:---|:---|:---|
| **Part A: $\Delta F$ vs. $L$** | $i = 4.0\text{ A}, I_m = 2.0\text{ A}$ | $0.0475 \pm 0.0011\text{ N/m}$ | $0.407 \pm 0.066\text{ mN}$ | $0.9989$ | $B = 11.87 \pm 0.32\text{ mT}$ |
| **Part B: $\Delta F$ vs. $i$** | $L = 10\text{ cm}, I_m = 2.0\text{ A}$ | $1.024 \pm 0.190\text{ mN/A}$ | $-1.120 \pm 0.520\text{ mN}$ | $0.9356$ | $B = 10.24 \pm 1.90\text{ mT}$ |
| **Part C: $\Delta F$ vs. $I_m$** | $L = 10\text{ cm}, i = 4.0\text{ A}$ | $1.972 \pm 0.308\text{ mN/A}$ | $-0.930 \pm 0.422\text{ mN}$ | $0.9534$ | $\alpha = 4.93 \pm 0.77\text{ mT/A}$<br>$B(2\text{ A}) = 9.86 \pm 1.55\text{ mT}$ |
| **Weighted Average** | $I_m = 2.0\text{ A}$ | — | — | — | $\mathbf{B = 11.75 \pm 0.31\text{ mT}}$ |

---

## Files

- `5.tex` — Full Persian LaTeX lab report source (XePersian)
- `5.pdf` — Compiled PDF document
- `analysis.py` — Python script for OLS regressions, uncertainty propagation, and figure generation
- `plots/F_vs_L.pdf` — Plot of $\Delta F$ vs. $L$ with error bars and linear regression line
- `plots/F_vs_i.pdf` — Plot of $\Delta F$ vs. $i$ with error bars and linear regression line
- `plots/F_vs_Im.pdf` — Plot of $\Delta F$ vs. $I_m$ with error bars and linear regression line
- `README.md` — This documentation file

---

## How to Compile

Because the report is written in Persian using the `xepersian` package and `Dibaj` font, compile with **XeLaTeX**:

```bash
xelatex 5.tex
```

To run the data analysis and regenerate all plots:

```bash
python analysis.py
```

---

## Notes

- **Direction of Magnetic Force:** The current direction in the loop is intentionally chosen such that $\vec{F} = i(\vec{L} \times \vec{B})$ points downward. A downward force adds directly to the apparent weight of the balance beam, ensuring stable contact with the knife-edge support and higher balance sensitivity. An upward force risks unloading the beam, introducing mechanical instability or chatter.
- **Role of Flexible Wires:** Bare, thin, highly flexible copper leads are used without insulation to minimize mechanical stiffness and thermal expansion distortion. The leads must hang with a slight catenary curve ("sag") so that vertical beam deflections are unconstrained.
- **Zero Intercept Discrepancies:** Small non-zero intercepts ($b \approx 0.41\text{ mN}$ in Part A, $-1.12\text{ mN}$ in Part B, $-0.93\text{ mN}$ in Part C) are attributed to remnant magnetic field in the core, balance zero-drift, and edge fringing fields at the loop corners. Taking the slope $\Delta F / \Delta x$ rather than single-point ratios effectively eliminates these constant zero offsets.
