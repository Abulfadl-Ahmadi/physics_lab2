"""
Ohm's Law and Wire Resistance Parameters Experiment Analysis
Sharif University of Technology - Department of Physics
Lab-Phy-II: Experiment 1

This script performs:
1. Verification of Ohm's Law (V vs I) using OLS linear regression and error analysis.
2. Study of resistance dependence on wire length: R = f(l).
3. Study of resistance dependence on cross-sectional area: R = f(1/S).
4. Determination of electrical resistivity (rho) across different materials (Nichrome, Galvanized steel, Pure Chromium).
5. Comprehensive error propagation for all derived quantities using partial derivatives.
6. Publication-grade vector PDF plots saved to the 'plots/' directory.
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 0. Setup and Helpers
# ----------------------------------------------------------------------
plots_dir = Path('plots')
plots_dir.mkdir(exist_ok=True)

# Matplotlib publication-grade styling
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.family': 'serif',
    'mathtext.fontset': 'dejavuserif',
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'grid.alpha': 0.4,
    'grid.linestyle': '--'
})


def ols_fit(x, y):
    """
    Ordinary Least Squares (OLS) linear regression.
    Model: y = m * x + c
    
    Returns:
        slope (m), se_slope (se_m), intercept (c), se_intercept (se_c), r2
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    N = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    s_yx = np.sqrt(ss_res / (N - 2)) if N > 2 else 0.0
    ss_x = np.sum((x - np.mean(x))**2)
    se_m = s_yx / np.sqrt(ss_x) if ss_x > 0 else 0.0
    se_c = s_yx * np.sqrt(1.0 / N + np.mean(x)**2 / ss_x) if ss_x > 0 else 0.0
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    return m, se_m, c, se_c, r2


# ----------------------------------------------------------------------
# 1. Experimental Uncertainties (Instrument Resolution)
# ----------------------------------------------------------------------
delta_V = 0.01        # Voltmeter uncertainty: 0.01 V
delta_I = 1.0e-3      # Ammeter uncertainty: 1 mA = 0.001 A
delta_l = 1.0e-3      # Ruler length uncertainty: 1 mm = 0.001 m
delta_d = 0.01e-3     # Micrometer diameter uncertainty: 0.01 mm = 1e-5 m

print("=" * 70)
print("EXPERIMENT 1: OHM'S LAW AND RESISTANCE PARAMETERS")
print("Department of Physics, Sharif University of Technology")
print("=" * 70)

# ----------------------------------------------------------------------
# Part 1: Verification of Ohm's Law (V vs I)
# Wire 1: Nichrome (NiCr), Length L = 1.00 m, Diameter d = 0.25 mm
# ----------------------------------------------------------------------
I1_mA = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
I1 = I1_mA * 1e-3  # Convert to Amperes
V1 = np.array([2.62, 5.31, 7.98, 10.63, 13.27])  # Volts

dI1 = np.full_like(I1, delta_I)
dV1 = np.full_like(V1, delta_V)

# Pointwise resistances and uncertainty propagation:
# R = V / I  ==>  delta_R = R * sqrt((delta_V/V)^2 + (delta_I/I)^2)
R_pointwise = V1 / I1
dR_pointwise = R_pointwise * np.sqrt((dV1 / V1)**2 + (dI1 / I1)**2)

# OLS Fit: V = R * I + V0
m1, se_m1, c1, se_c1, r2_1 = ols_fit(I1, V1)

print("\n--- Part 1: V vs I (Ohm's Law Verification) ---")
print(f"OLS Slope (Resistance R) : {m1:.4f} +/- {se_m1:.4f} Ohm")
print(f"OLS Intercept (V0)       : {c1:.4f} +/- {se_c1:.4f} V")
print(f"Determination Coeff (R2) : {r2_1:.6f}")
print("Pointwise R values:")
for i, (i_val, v_val, r_val, dr_val) in enumerate(zip(I1_mA, V1, R_pointwise, dR_pointwise)):
    rel_err = (abs(r_val - m1) / m1) * 100
    print(f"  Point {i+1}: I = {i_val:5.1f} mA, V = {v_val:5.2f} V -> R = {r_val:5.2f} +/- {dr_val:4.2f} Ohm (Dev: {rel_err:4.2f}%)")

# Plot Part 1: V vs I
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

# Fit line
I_fit = np.linspace(0.08, 0.52, 200)
V_fit = m1 * I_fit + c1

ax1.errorbar(I1 * 1e3, V1, xerr=dI1 * 1e3, yerr=dV1, fmt='o', color='navy',
             ecolor='cornflowerblue', elinewidth=1.5, capsize=4, capthick=1.5,
             label=r'Measured Data ($Wire\ 1:\ \mathrm{NiCr},\ d=0.25\,\mathrm{mm}$)')
ax1.plot(I_fit * 1e3, V_fit, 'r-', linewidth=2.0,
         label=f'OLS Fit: $V = ({m1:.2f} \\pm {se_m1:.2f})I + ({c1:.3f} \\pm {se_c1:.3f})$\n$R^2 = {r2_1:.6f}$')
ax1.set_ylabel(r'Voltage $V$ [$\mathrm{V}$]')
ax1.set_title(r"Ohm's Law Verification: Voltage vs. Current ($V - I$)")
ax1.grid(True)
ax1.legend(loc='upper left', frameon=True, framealpha=0.9)

# Residuals subplot
res1 = V1 - (m1 * I1 + c1)
ax2.axhline(0, color='gray', linestyle='--', linewidth=1.2)
ax2.errorbar(I1 * 1e3, res1 * 1e3, yerr=dV1 * 1e3, fmt='s', color='crimson',
             ecolor='darkred', elinewidth=1.5, capsize=4, capthick=1.5)
ax2.set_xlabel(r'Current $I$ [$\mathrm{mA}$]')
ax2.set_ylabel(r'Residuals [$\mathrm{mV}$]')
ax2.grid(True)
plt.tight_layout()
plt.savefig(plots_dir / 'V_vs_I_ohms_law.pdf', bbox_inches='tight')
plt.close()


# ----------------------------------------------------------------------
# Part 2: Resistance vs. Wire Length R = f(l)
# Wire 1: Nichrome, d = 0.25 mm, I_const = 250 mA
# ----------------------------------------------------------------------
l2_cm = np.array([10.0, 27.0, 50.0, 80.0, 100.0])
l2 = l2_cm * 1e-2  # Convert to meters
V2 = np.array([0.68, 1.79, 3.30, 5.32, 6.61])  # Volts
I_const2 = 0.250  # Amperes
dI_const2 = delta_I

# Resistance: R = V / I_const
R2 = V2 / I_const2
dR2 = R2 * np.sqrt((delta_V / V2)**2 + (dI_const2 / I_const2)**2)
dl2 = np.full_like(l2, delta_l)

# OLS Fit: R = m_l * l + c_l
m_l, se_m_l, c_l, se_c_l, r2_l = ols_fit(l2, R2)

# Resistivity from length slope:
# R = rho * l / S  ==>  m_l = rho / S  ==>  rho_1 = m_l * S
# S1 = pi * (d1/2)^2
d1 = 0.25e-3  # m
S1 = np.pi * (d1 / 2.0)**2
dS1 = 2.0 * S1 * (delta_d / d1)

rho_1 = m_l * S1
# Uncertainty propagation:
# delta_rho / rho = sqrt((se_m_l / m_l)^2 + (dS1 / S1)^2) = sqrt((se_m_l / m_l)^2 + (2 * delta_d / d1)^2)
drho_1 = rho_1 * np.sqrt((se_m_l / m_l)**2 + (2.0 * delta_d / d1)**2)

print("\n--- Part 2: Resistance vs. Length R = f(l) ---")
print(f"OLS Slope (dR/dl)        : {m_l:.4f} +/- {se_m_l:.4f} Ohm/m")
print(f"OLS Intercept (R0)       : {c_l:.4f} +/- {se_c_l:.4f} Ohm")
print(f"Determination Coeff (R2) : {r2_l:.6f}")
print(f"Calculated Wire Area S1  : {S1:.4e} +/- {dS1:.4e} m^2")
print(f"Resistivity (rho_1, NiCr): ({rho_1*1e6:.3f} +/- {drho_1*1e6:.3f}) x 10^-6 Ohm*m")

# Plot Part 2: R vs Length
plt.figure(figsize=(7, 5))
l_fit = np.linspace(0.05, 1.05, 200)
R_fit_l = m_l * l_fit + c_l

plt.errorbar(l2_cm, R2, xerr=dl2 * 1e2, yerr=dR2, fmt='o', color='darkgreen',
             ecolor='limegreen', elinewidth=1.5, capsize=4, capthick=1.5,
             label=r'Measured Data ($I = 250\,\mathrm{mA}$)')
plt.plot(l_fit * 1e2, R_fit_l, 'r-', linewidth=2.0,
         label=f'OLS Fit: $R = ({m_l:.2f} \\pm {se_m_l:.2f})l + ({c_l:.3f} \\pm {se_c_l:.3f})$\n$R^2 = {r2_l:.6f}$')
plt.xlabel(r'Wire Length $l$ [$\mathrm{cm}$]')
plt.ylabel(r'Resistance $R$ [$\Omega$]')
plt.title(r'Resistance vs. Wire Length ($R - l$)')
plt.grid(True)
plt.legend(loc='upper left', frameon=True, framealpha=0.9)
plt.tight_layout()
plt.savefig(plots_dir / 'R_vs_length.pdf', bbox_inches='tight')
plt.close()


# ----------------------------------------------------------------------
# Part 3: Resistance vs. Cross-Sectional Area R = f(1/S)
# Wires 1, 2, 3: All Nichrome, Length L = 1.00 m, I_const = 250 mA
# Wire 1: d = 0.25 mm; Wire 2: d = 0.40 mm; Wire 3: d = 0.30 mm
# ----------------------------------------------------------------------
wire_diameters_mm = np.array([0.25, 0.40, 0.30])
wire_diameters = wire_diameters_mm * 1e-3  # meters
V3 = np.array([6.62, 2.13, 3.80])          # Volts
I_const3 = 0.250
L_const = 1.00                             # meters
dL_const = delta_l

R3 = V3 / I_const3
dR3 = R3 * np.sqrt((delta_V / V3)**2 + (delta_I / I_const3)**2)

# Cross-sectional areas: S = pi * (d / 2)^2
S3 = np.pi * (wire_diameters / 2.0)**2
dS3 = 2.0 * S3 * (delta_d / wire_diameters)

# Inverse areas: xi = 1 / S
inv_S = 1.0 / S3
d_inv_S = 2.0 * inv_S * (delta_d / wire_diameters)

# OLS Fit: R = m_S * (1/S) + c_S
m_S, se_m_S, c_S, se_c_S, r2_S = ols_fit(inv_S, R3)

# Resistivity from area slope:
# R = rho * L * (1/S)  ==>  m_S = rho * L  ==>  rho_2 = m_S / L
rho_2 = m_S / L_const
drho_2 = rho_2 * np.sqrt((se_m_S / m_S)**2 + (dL_const / L_const)**2)

print("\n--- Part 3: Resistance vs. Inverse Cross-Sectional Area R = f(1/S) ---")
print(f"OLS Slope (m_S = rho*L)  : {m_S:.4e} +/- {se_m_S:.4e} Ohm*m^2")
print(f"OLS Intercept (c_S)      : {c_S:.4f} +/- {se_c_S:.4f} Ohm")
print(f"Determination Coeff (R2) : {r2_S:.6f}")
print(f"Resistivity (rho_2, NiCr): ({rho_2*1e6:.3f} +/- {drho_2*1e6:.3f}) x 10^-6 Ohm*m")

for i, (d_mm, s_val, inv_s_val, r_val, dr_val) in enumerate(zip(wire_diameters_mm, S3, inv_S, R3, dR3)):
    rho_direct = r_val * s_val / L_const
    drho_direct = rho_direct * np.sqrt((dr_val / r_val)**2 + (2 * delta_d / (d_mm * 1e-3))**2 + (dL_const / L_const)**2)
    print(f"  Wire {i+1} (d = {d_mm:.2f} mm): S = {s_val*1e6:.4f} mm^2, 1/S = {inv_s_val*1e-6:5.2f} mm^-2 -> R = {r_val:5.2f} +/- {dr_val:4.2f} Ohm, rho = ({rho_direct*1e6:.3f} +/- {drho_direct*1e6:.3f}) x 10^-6 Ohm*m")

# Plot Part 3: R vs 1/S
plt.figure(figsize=(7, 5))
inv_S_fit = np.linspace(5e6, 23e6, 200)
R_fit_S = m_S * inv_S_fit + c_S

plt.errorbar(inv_S * 1e-6, R3, xerr=d_inv_S * 1e-6, yerr=dR3, fmt='s', color='purple',
             ecolor='mediumorchid', elinewidth=1.5, capsize=4, capthick=1.5,
             label=r'Measured Wires ($\mathrm{NiCr},\ L = 1.00\,\mathrm{m}$)')
plt.plot(inv_S_fit * 1e-6, R_fit_S, 'r-', linewidth=2.0,
         label=f'OLS Fit: $R = ({m_S*1e6:.2f} \\pm {se_m_S*1e6:.2f}) \\cdot (1/S) + ({c_S:.2f} \\pm {se_c_S:.2f})$\n$R^2 = {r2_S:.4f}$')
plt.xlabel(r'Inverse Cross-Sectional Area $1/S$ [$\mathrm{mm}^{-2}$]')
plt.ylabel(r'Resistance $R$ [$\Omega$]')
plt.title(r'Resistance vs. Inverse Area ($R - 1/S$)')
plt.grid(True)
plt.legend(loc='upper left', frameon=True, framealpha=0.9)
plt.tight_layout()
plt.savefig(plots_dir / 'R_vs_inv_area.pdf', bbox_inches='tight')
plt.close()


# ----------------------------------------------------------------------
# Part 4: Resistivity across Materials R = f(rho)
# Wires 3, 4, 5: Length L = 1.00 m, I_const = 250 mA
# Wire 3: Nichrome (NiCr), d = 0.30 mm
# Wire 4: Galvanized Iron (گالوانیزه), d = 0.30 mm
# Wire 5: Pure Chromium (کروم خالص), d = 0.40 mm
# ----------------------------------------------------------------------
material_names = ['Nichrome (NiCr)', 'Galvanized steel', 'Pure Chromium (Cr)']
d_mat_mm = np.array([0.30, 0.30, 0.40])
d_mat = d_mat_mm * 1e-3
V_mat = np.array([3.81, 0.86, 3.00])
I_const4 = 0.250

R_mat = V_mat / I_const4
dR_mat = R_mat * np.sqrt((delta_V / V_mat)**2 + (delta_I / I_const4)**2)

S_mat = np.pi * (d_mat / 2.0)**2
dS_mat = 2.0 * S_mat * (delta_d / d_mat)

rho_mat = R_mat * S_mat / L_const
drho_mat = rho_mat * np.sqrt(
    (delta_V / V_mat)**2 +
    (delta_I / I_const4)**2 +
    (2.0 * delta_d / d_mat)**2 +
    (dL_const / L_const)**2
)

print("\n--- Part 4: Material Resistivity Determination ---")
for name, d_val, r_val, dr_val, rho_val, drho_val in zip(material_names, d_mat_mm, R_mat, dR_mat, rho_mat, drho_mat):
    print(f"{name:20s} (d = {d_val:.2f} mm): R = {r_val:5.2f} +/- {dr_val:4.2f} Ohm -> rho = ({rho_val*1e6:.4f} +/- {drho_val*1e6:.4f}) x 10^-6 Ohm*m")

# Combined / Weighted Nichrome Resistivity
# We have rho_1 (from length slope), rho_2 (from area slope), and rho_direct (wire 3)
rho_nicr_estimates = np.array([rho_1, rho_2, rho_mat[0]])
drho_nicr_estimates = np.array([drho_1, drho_2, drho_mat[0]])
weights = 1.0 / drho_nicr_estimates**2
rho_nicr_weighted = np.sum(weights * rho_nicr_estimates) / np.sum(weights)
drho_nicr_weighted = 1.0 / np.sqrt(np.sum(weights))

print(f"\n--- Combined Nichrome Resistivity ---")
print(f"From Length Regression : ({rho_1*1e6:.3f} +/- {drho_1*1e6:.3f}) x 10^-6 Ohm*m")
print(f"From Area Regression   : ({rho_2*1e6:.3f} +/- {drho_2*1e6:.3f}) x 10^-6 Ohm*m")
print(f"From Direct Wire 3     : ({rho_mat[0]*1e6:.3f} +/- {drho_mat[0]*1e6:.3f}) x 10^-6 Ohm*m")
print(f"Weighted Mean (NiCr)   : ({rho_nicr_weighted*1e6:.3f} +/- {drho_nicr_weighted*1e6:.3f}) x 10^-6 Ohm*m")

# Plot Part 4: Resistivity Comparison
plt.figure(figsize=(7, 5))
x_pos = np.arange(len(material_names))
bars = plt.bar(x_pos, rho_mat * 1e6, yerr=drho_mat * 1e6, capsize=5,
               color=['#1f77b4', '#2ca02c', '#d62728'], alpha=0.85, edgecolor='black', linewidth=1.2)

# Add value labels above bars
for bar, rho_v, drho_v in zip(bars, rho_mat * 1e6, drho_mat * 1e6):
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2.0, yval + drho_v + 0.05,
             f"${rho_v:.3f} \\pm {drho_v:.3f}$", ha='center', va='bottom', fontsize=9.5)

plt.xticks(x_pos, material_names)
plt.ylabel(r'Resistivity $\rho$ [$10^{-6}\ \Omega\cdot\mathrm{m}$]')
plt.title(r'Experimental Electrical Resistivity by Material')
plt.ylim(0, 2.0)
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(plots_dir / 'resistivity_comparison.pdf', bbox_inches='tight')
plt.close()


# ----------------------------------------------------------------------
# Comprehensive 4-Panel Summary Figure
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel (a): V vs I
axes[0, 0].errorbar(I1 * 1e3, V1, xerr=dI1 * 1e3, yerr=dV1, fmt='o', color='navy',
                    ecolor='cornflowerblue', capsize=3, label='Data')
axes[0, 0].plot(I_fit * 1e3, V_fit, 'r-',
                label=f'$V = {m1:.2f}I {c1:+.3f}$\n$R^2 = {r2_1:.5f}$')
axes[0, 0].set_xlabel(r'$I$ [$\mathrm{mA}$]')
axes[0, 0].set_ylabel(r'$V$ [$\mathrm{V}$]')
axes[0, 0].set_title(r'(a) Ohm\'s Law: $V$ vs. $I$')
axes[0, 0].grid(True)
axes[0, 0].legend(loc='upper left')

# Panel (b): R vs l
axes[0, 1].errorbar(l2_cm, R2, xerr=dl2 * 1e2, yerr=dR2, fmt='o', color='darkgreen',
                    ecolor='limegreen', capsize=3, label='Data')
axes[0, 1].plot(l_fit * 1e2, R_fit_l, 'r-',
                label=f'$R = {m_l:.2f}l {c_l:+.3f}$\n$R^2 = {r2_l:.5f}$')
axes[0, 1].set_xlabel(r'Length $l$ [$\mathrm{cm}$]')
axes[0, 1].set_ylabel(r'$R$ [$\Omega$]')
axes[0, 1].set_title(r'(b) Length Dependence: $R$ vs. $l$')
axes[0, 1].grid(True)
axes[0, 1].legend(loc='upper left')

# Panel (c): R vs 1/S
axes[1, 0].errorbar(inv_S * 1e-6, R3, xerr=d_inv_S * 1e-6, yerr=dR3, fmt='s', color='purple',
                    ecolor='mediumorchid', capsize=3, label='Data')
axes[1, 0].plot(inv_S_fit * 1e-6, R_fit_S, 'r-',
                label=f'$R = {m_S*1e6:.2f}(1/S) {c_S:+.2f}$\n$R^2 = {r2_S:.4f}$')
axes[1, 0].set_xlabel(r'$1/S$ [$\mathrm{mm}^{-2}$]')
axes[1, 0].set_ylabel(r'$R$ [$\Omega$]')
axes[1, 0].set_title(r'(c) Area Dependence: $R$ vs. $1/S$')
axes[1, 0].grid(True)
axes[1, 0].legend(loc='upper left')

# Panel (d): Material Resistivity
bars = axes[1, 1].bar(x_pos, rho_mat * 1e6, yerr=drho_mat * 1e6, capsize=4,
                      color=['#1f77b4', '#2ca02c', '#d62728'], alpha=0.85, edgecolor='black')
axes[1, 1].set_xticks(x_pos)
axes[1, 1].set_xticklabels(['NiCr', 'Galv.', 'Pure Cr'])
axes[1, 1].set_ylabel(r'$\rho$ [$10^{-6}\ \Omega\cdot\mathrm{m}$]')
axes[1, 1].set_title(r'(d) Material Comparison: $\rho$')
axes[1, 1].grid(axis='y')
for bar, rho_v, drho_v in zip(bars, rho_mat * 1e6, drho_mat * 1e6):
    yval = bar.get_height()
    axes[1, 1].text(bar.get_x() + bar.get_width()/2.0, yval + drho_v + 0.04,
                    f'{rho_v:.2f}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(plots_dir / 'ohm_law_summary.pdf', bbox_inches='tight')
plt.close()

print("\nAll plots successfully generated and saved to 'plots/' directory:")
print("  1. plots/V_vs_I_ohms_law.pdf")
print("  2. plots/R_vs_length.pdf")
print("  3. plots/R_vs_inv_area.pdf")
print("  4. plots/resistivity_comparison.pdf")
print("  5. plots/ohm_law_summary.pdf")
print("=" * 70)
