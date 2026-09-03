"""
Kirchhoff's Laws and Wheatstone Bridge Experiment Analysis
Sharif University of Technology - Department of Physics
Lab-Phy-II: Experiment 2

This script performs:
1. Component-level Ohm's law verification with uncertainty propagation
2. Kirchhoff's Current Law (KCL) verification at circuit nodes
3. Kirchhoff's Voltage Law (KVL) verification around independent loops
4. Full theoretical nodal analysis of the 2-source 5-resistor network
5. Unknown resistance determination: Direct (V/I) vs Balanced Wheatstone Bridge (R1*R3/R2)
6. Publication-grade vector plots saved to plots/
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Ensure UTF-8 output on Windows console
if sys.stdout.encoding != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

# ----------------------------------------------------------------------
# 0. Setup and Helpers
# ----------------------------------------------------------------------
plots_dir = Path('plots')
plots_dir.mkdir(exist_ok=True)

# Set matplotlib publication style
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
# 1. Raw Data & Experimental Uncertainties
# ----------------------------------------------------------------------
# Multimeter measured resistances (Ohms)
resistor_names = ['R1', 'R2', 'R3', 'R4', 'R5']
R = np.array([388.0, 219.0, 47.0, 46.0, 99.0])
delta_R = np.full_like(R, 1.0)  # Digital multimeter 1 Ohm resolution/uncertainty

# Power supply voltages (Volts)
V1 = 5.00
delta_V1 = 0.01
V2 = 8.00
delta_V2 = 0.01

# Measured branch currents (mA and A)
I_exp_mA = np.array([20.69, 10.16, 59.47, 49.57, 29.16])
delta_I_mA = np.full_like(I_exp_mA, 0.01)  # 0.01 mA ammeter resolution
I_exp = I_exp_mA * 1e-3  # A
delta_I = delta_I_mA * 1e-3  # A

# Measured branch voltages across resistors (Volts)
V_exp = np.array([7.96, 2.24, 2.82, 2.28, 2.89])
delta_V = np.full_like(V_exp, 0.01)  # 0.01 V voltmeter resolution

# Unknown resistor measurements (Rheostat Rx in balanced bridge)
V_Rx = 0.886
delta_V_Rx = 0.001  # precision voltmeter
I_Rx_mA = 10.93
delta_I_Rx_mA = 0.01
I_Rx = I_Rx_mA * 1e-3
delta_I_Rx = delta_I_Rx_mA * 1e-3


# ----------------------------------------------------------------------
# 2. Ohm's Law Verification & Branch Discrepancies
# ----------------------------------------------------------------------
# Theoretical current: I_theo = V_exp / R
I_theo = V_exp / R  # A
I_theo_mA = I_theo * 1e3
# Error propagation: delta_I_theo = I_theo * sqrt((delta_V / V)^2 + (delta_R / R)^2)
delta_I_theo = I_theo * np.sqrt((delta_V / V_exp)**2 + (delta_R / R)**2)
delta_I_theo_mA = delta_I_theo * 1e3

# Theoretical voltage: V_theo = I_exp * R
V_theo = I_exp * R  # V
# Error propagation: delta_V_theo = V_theo * sqrt((delta_I / I)^2 + (delta_R / R)^2)
delta_V_theo = V_theo * np.sqrt((delta_I / I_exp)**2 + (delta_R / R)**2)

# Relative percentage discrepancies
pct_err_I = np.abs(I_exp_mA - I_theo_mA) / I_theo_mA * 100.0
pct_err_V = np.abs(V_exp - V_theo) / V_theo * 100.0

print("=" * 70)
print("1. OHM'S LAW COMPONENT VERIFICATION")
print("=" * 70)
header_str = f"{'Resistor':<8} | {'R (Ohm)':<10} | {'V_exp (V)':<12} | {'I_exp (mA)':<14} | {'I_theo (mA)':<14} | {'Err_I (%)':<10} | {'V_theo (V)':<12} | {'Err_V (%)':<10}"
print(header_str)
print("-" * len(header_str))
for i in range(len(R)):
    print(f"{resistor_names[i]:<8} | {R[i]:.0f} ± {delta_R[i]:.0f}    | "
          f"{V_exp[i]:.2f} ± {delta_V[i]:.2f}   | "
          f"{I_exp_mA[i]:.2f} ± {delta_I_mA[i]:.2f}   | "
          f"{I_theo_mA[i]:.3f} ± {delta_I_theo_mA[i]:.3f} | "
          f"{pct_err_I[i]:.3f}%     | "
          f"{V_theo[i]:.3f} ± {delta_V_theo[i]:.3f}  | "
          f"{pct_err_V[i]:.3f}%")
print()

# Global OLS fit of V_exp vs (I_exp * R)
# According to Ohm's Law: V = 1.0 * (I * R) + 0.0
IR_prod = I_exp * R
slope_ohm, se_slope_ohm, int_ohm, se_int_ohm, r2_ohm = ols_fit(IR_prod, V_exp)
print(f"Global Ohm's Law OLS Fit: V_meas = ({slope_ohm:.4f} ± {se_slope_ohm:.4f}) * (I*R) + ({int_ohm:.4f} ± {se_int_ohm:.4f})")
print(f"Coefficient of Determination R^2 = {r2_ohm:.6f}")
print()


# ----------------------------------------------------------------------
# 3. Kirchhoff's Current Law (KCL) Verification
# ----------------------------------------------------------------------
print("=" * 70)
print("2. KIRCHHOFF'S CURRENT LAW (KCL) VERIFICATION")
print("=" * 70)

# Node N:
# Entering: I_1 + I_5
# Leaving:  I_4
sum_I_in_N = I_exp_mA[0] + I_exp_mA[4]  # I1 + I5
delta_sum_I_in_N = np.sqrt(delta_I_mA[0]**2 + delta_I_mA[4]**2)
sum_I_out_N = I_exp_mA[3]  # I4
delta_sum_I_out_N = delta_I_mA[3]
delta_I_N = np.abs(sum_I_out_N - sum_I_in_N)
unc_diff_N = np.sqrt(delta_sum_I_in_N**2 + delta_sum_I_out_N**2)
pct_kcl_N = (delta_I_N / sum_I_out_N) * 100.0

print(f"Node N:")
print(f"  Sum(I_in)  = I1 + I5 = {I_exp_mA[0]:.2f} + {I_exp_mA[4]:.2f} = {sum_I_in_N:.2f} ± {delta_sum_I_in_N:.3f} mA")
print(f"  Sum(I_out) = I4      = {sum_I_out_N:.2f} ± {delta_sum_I_out_N:.3f} mA")
print(f"  Residual Delta_I     = {delta_I_N:.2f} ± {unc_diff_N:.3f} mA  (Relative Error: {pct_kcl_N:.2f}%)")

# Node M:
# Entering: I_2 + I_4
# Leaving:  I_3
sum_I_in_M = I_exp_mA[1] + I_exp_mA[3]  # I2 + I4
delta_sum_I_in_M = np.sqrt(delta_I_mA[1]**2 + delta_I_mA[3]**2)
sum_I_out_M = I_exp_mA[2]  # I3
delta_sum_I_out_M = delta_I_mA[2]
delta_I_M = np.abs(sum_I_out_M - sum_I_in_M)
unc_diff_M = np.sqrt(delta_sum_I_in_M**2 + delta_sum_I_out_M**2)
pct_kcl_M = (delta_I_M / sum_I_out_M) * 100.0

print(f"Node M:")
print(f"  Sum(I_in)  = I2 + I4 = {I_exp_mA[1]:.2f} + {I_exp_mA[3]:.2f} = {sum_I_in_M:.2f} ± {delta_sum_I_in_M:.3f} mA")
print(f"  Sum(I_out) = I3      = {sum_I_out_M:.2f} ± {delta_sum_I_out_M:.3f} mA")
print(f"  Residual Delta_I     = {delta_I_M:.2f} ± {unc_diff_M:.3f} mA  (Relative Error: {pct_kcl_M:.2f}%)")
print()


# ----------------------------------------------------------------------
# 4. Kirchhoff's Voltage Law (KVL) Verification
# ----------------------------------------------------------------------
print("=" * 70)
print("3. KIRCHHOFF'S VOLTAGE LAW (KVL) VERIFICATION")
print("=" * 70)

# Loop 1: Right loop (N -> M -> B -> N)
# Path: +V2 - VR4 + VR2 - VR1 = 0
kvl_1_sum = V2 - V_exp[3] + V_exp[1] - V_exp[0]
kvl_1_unc = np.sqrt(delta_V2**2 + delta_V[3]**2 + delta_V[1]**2 + delta_V[0]**2)
kvl_1_pct = (np.abs(kvl_1_sum) / V2) * 100.0

# Loop 2: Left loop (N -> M -> A -> N)
# Path: +V2 - VR4 - VR3 - VR5 = 0
kvl_2_sum = V2 - V_exp[3] - V_exp[2] - V_exp[4]
kvl_2_unc = np.sqrt(delta_V2**2 + delta_V[3]**2 + delta_V[2]**2 + delta_V[4]**2)
kvl_2_pct = (np.abs(kvl_2_sum) / V2) * 100.0

# Loop 3: Top loop with V1 (A -> B -> M -> A)
# Path: +V1 - VR2 - VR3 = 0
kvl_3_sum = V1 - V_exp[1] - V_exp[2]
kvl_3_unc = np.sqrt(delta_V1**2 + delta_V[1]**2 + delta_V[2]**2)
kvl_3_pct = (np.abs(kvl_3_sum) / V1) * 100.0

# Loop 4: Bottom loop with V1 (A -> B -> N -> A)
# Path: +V1 - VR1 + VR5 = 0
kvl_4_sum = V1 - V_exp[0] + V_exp[4]
kvl_4_unc = np.sqrt(delta_V1**2 + delta_V[0]**2 + delta_V[4]**2)
kvl_4_pct = (np.abs(kvl_4_sum) / V1) * 100.0

print(f"Loop 1 (Right inner loop with V2):")
print(f"  V2 - V_R4 + V_R2 - V_R1 = {V2:.2f} - {V_exp[3]:.2f} + {V_exp[1]:.2f} - {V_exp[0]:.2f} = {kvl_1_sum:+.3f} ± {kvl_1_unc:.3f} V (Err: {kvl_1_pct:.2f}%)")

print(f"Loop 2 (Left inner loop with V2):")
print(f"  V2 - V_R4 - V_R3 - V_R5 = {V2:.2f} - {V_exp[3]:.2f} - {V_exp[2]:.2f} - {V_exp[4]:.2f} = {kvl_2_sum:+.3f} ± {kvl_2_unc:.3f} V (Err: {kvl_2_pct:.3f}%)")

print(f"Loop 3 (Top loop with V1):")
print(f"  V1 - V_R2 - V_R3 = {V1:.2f} - {V_exp[1]:.2f} - {V_exp[2]:.2f} = {kvl_3_sum:+.3f} ± {kvl_3_unc:.3f} V (Err: {kvl_3_pct:.2f}%)")

print(f"Loop 4 (Bottom loop with V1):")
print(f"  V1 - V_R1 + V_R5 = {V1:.2f} - {V_exp[0]:.2f} + {V_exp[4]:.2f} = {kvl_4_sum:+.3f} ± {kvl_4_unc:.3f} V (Err: {kvl_4_pct:.2f}%)")
print()


# ----------------------------------------------------------------------
# 5. Full Theoretical Nodal Circuit Analysis
# ----------------------------------------------------------------------
print("=" * 70)
print("4. FULL THEORETICAL NODAL NETWORK ANALYSIS")
print("=" * 70)

# Reference: V_A = 0 V, V_B = V1 = 5.0 V
# Unknown node potentials: V_N, V_M
# Node N equation: (VB - VN)/R1 + (VA - VN)/R5 = (VN + V2 - VM)/R4
# Node M equation: (VB - VM)/R2 + (VN + V2 - VM)/R4 = (VM - VA)/R3
G = np.array([
    [-1.0/R[0] - 1.0/R[4] - 1.0/R[3], 1.0/R[3]],
    [1.0/R[3], -1.0/R[1] - 1.0/R[3] - 1.0/R[2]]
])
rhs = np.array([
    V2 / R[3] - V1 / R[0],
    -V2 / R[3] - V1 / R[1]
])
VN_theo, VM_theo = np.linalg.solve(G, rhs)

# Branch currents
I1_sim = (V1 - VN_theo) / R[0]
I2_sim = (V1 - VM_theo) / R[1]
I3_sim = (VM_theo - 0.0) / R[2]
I4_sim = (VN_theo + V2 - VM_theo) / R[3]
I5_sim = (0.0 - VN_theo) / R[4]
I_sim_mA = np.array([I1_sim, I2_sim, I3_sim, I4_sim, I5_sim]) * 1e3

# Branch voltages
V_sim = np.array([
    np.abs(V1 - VN_theo),
    np.abs(V1 - VM_theo),
    np.abs(VM_theo),
    np.abs(VN_theo + V2 - VM_theo),
    np.abs(VN_theo)
])

print(f"Theoretical Node Potentials: V_N = {VN_theo:.4f} V, V_M = {VM_theo:.4f} V")
print(f"{'Branch':<8} | {'I_sim (mA)':<12} | {'I_exp (mA)':<12} | {'Diff_I (mA)':<12} | {'V_sim (V)':<12} | {'V_exp (V)':<12} | {'Diff_V (V)':<12}")
print("-" * 86)
for i in range(5):
    print(f"{resistor_names[i]:<8} | {I_sim_mA[i]:<12.2f} | {I_exp_mA[i]:<12.2f} | "
          f"{I_exp_mA[i] - I_sim_mA[i]:<+12.2f} | {V_sim[i]:<12.3f} | {V_exp[i]:<12.2f} | "
          f"{V_exp[i] - V_sim[i]:<+12.3f}")
print()


# ----------------------------------------------------------------------
# 6. Unknown Resistance Determination: Ohm's Law vs Wheatstone Bridge
# ----------------------------------------------------------------------
print("=" * 70)
print("5. UNKNOWN RESISTANCE (Rx) DETERMINATION")
print("=" * 70)

# Method 1: Direct Ohm's Law Rx1 = V_Rx / I_Rx
Rx1 = V_Rx / I_Rx
delta_Rx1 = Rx1 * np.sqrt((delta_V_Rx / V_Rx)**2 + (delta_I_Rx / I_Rx)**2)

# Method 2: Balanced Wheatstone Bridge Rx2 = (R1 * R3) / R2
Rx2 = (R[0] * R[2]) / R[1]
delta_Rx2 = Rx2 * np.sqrt((delta_R[0] / R[0])**2 + (delta_R[2] / R[2])**2 + (delta_R[1] / R[1])**2)

# Comparison
delta_Rx = np.abs(Rx1 - Rx2)
unc_delta_Rx = np.sqrt(delta_Rx1**2 + delta_Rx2**2)
rel_diff_Rx = (delta_Rx / Rx2) * 100.0
z_score = delta_Rx / unc_delta_Rx

print(f"Method 1 (Direct Ohm's Law):")
print(f"  Rx1 = V_Rx / I_Rx = {V_Rx:.3f} V / {I_Rx_mA:.2f} mA = {Rx1:.2f} ± {delta_Rx1:.2f} Ω")
print(f"Method 2 (Wheatstone Bridge Balance):")
print(f"  Rx2 = (R1 * R3) / R2 = ({R[0]:.0f} * {R[2]:.0f}) / {R[1]:.0f} = {Rx2:.2f} ± {delta_Rx2:.2f} Ω")
print(f"Comparison:")
print(f"  Difference |Rx1 - Rx2|  = {delta_Rx:.2f} ± {unc_delta_Rx:.2f} Ω")
print(f"  Relative Percentage Diff = {rel_diff_Rx:.2f}%")
print(f"  Compatibility Index (z)  = {z_score:.2f} (values match within ~{z_score:.1f} sigma)")
print()


# ----------------------------------------------------------------------
# 7. Visualization & Plot Generation
# ----------------------------------------------------------------------
print("=" * 70)
print("6. GENERATING PUBLICATION PLOTS")
print("=" * 70)

# Plot 1: Ohm's Law Verification across Network Resistors
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6.5), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

# Main plot
x_fit = np.linspace(0, 8.5, 100)
y_fit = slope_ohm * x_fit + int_ohm

ax1.errorbar(IR_prod, V_exp, xerr=delta_V_theo, yerr=delta_V, fmt='o', color='#1f77b4',
             ecolor='#1f77b4', elinewidth=1.2, capsize=3, label='Experimental Branches (R1-R5)', zorder=4)

# Plot Rx1 as well
IR_prod_Rx = I_Rx * Rx1
ax1.errorbar([IR_prod_Rx], [V_Rx], xerr=[I_Rx*delta_Rx1], yerr=[delta_V_Rx], fmt='s', color='#d62728',
             ecolor='#d62728', elinewidth=1.2, capsize=3, label=f'Unknown Resistor $R_x$ ({Rx1:.1f} $\\Omega$)', zorder=4)

ax1.plot(x_fit, y_fit, 'r--', label=f'OLS Fit: $V = ({slope_ohm:.4f} \\pm {se_slope_ohm:.4f}) \\cdot (IR) + ({int_ohm:.3f} \\pm {se_int_ohm:.3f})$\n$R^2 = {r2_ohm:.6f}$')
ax1.plot(x_fit, x_fit, 'k:', alpha=0.6, label='Ideal Ohm\'s Law ($V = IR$, slope=1.0)')

for i in range(len(R)):
    ax1.annotate(resistor_names[i], (IR_prod[i], V_exp[i]), textcoords="offset points", xytext=(8, -5), fontsize=9, fontweight='bold')
ax1.annotate('$R_x$', (IR_prod_Rx, V_Rx), textcoords="offset points", xytext=(8, -5), fontsize=9, fontweight='bold', color='#d62728')

ax1.set_ylabel('Measured Potential Drop $V_{\\mathrm{exp}}$ [V]')
ax1.set_title("Ohm's Law Network Verification ($V$ vs $I \\cdot R$)")
ax1.legend(loc='upper left', framealpha=0.9)
ax1.grid(True)

# Residuals plot
residuals = V_exp - (slope_ohm * IR_prod + int_ohm)
ax2.axhline(0, color='r', linestyle='--', linewidth=1.2)
ax2.errorbar(IR_prod, residuals * 1e3, yerr=delta_V * 1e3, fmt='o', color='#1f77b4', capsize=3)
ax2.set_xlabel('Calculated Voltage Product $I_{\\mathrm{exp}} \\cdot R$ [V]')
ax2.set_ylabel('Residuals [mV]')
ax2.grid(True)

plt.tight_layout()
p1_path = plots_dir / 'ohms_law_verification.pdf'
plt.savefig(p1_path, dpi=300)
plt.close()
print(f"Saved: {p1_path}")


# Plot 2: KCL Node Current Conservation
fig, ax = plt.subplots(figsize=(6.5, 4.5))

nodes = ['Node N', 'Node M']
x_indices = np.arange(len(nodes))
width = 0.3

I_in_vals = [sum_I_in_N, sum_I_in_M]
I_in_errs = [delta_sum_I_in_N, delta_sum_I_in_M]
I_out_vals = [sum_I_out_N, sum_I_out_M]
I_out_errs = [delta_sum_I_out_N, delta_sum_I_out_M]

rects1 = ax.bar(x_indices - width/2, I_in_vals, width, yerr=I_in_errs, capsize=4,
               label=r'$\sum I_{\mathrm{in}}$ (Entering)', color='#2ca02c', alpha=0.85, edgecolor='black')
rects2 = ax.bar(x_indices + width/2, I_out_vals, width, yerr=I_out_errs, capsize=4,
               label=r'$\sum I_{\mathrm{out}}$ (Leaving)', color='#1f77b4', alpha=0.85, edgecolor='black')

ax.set_ylabel('Current [mA]')
ax.set_title("Kirchhoff's Current Law (KCL) Verification at Circuit Nodes")
ax.set_xticks(x_indices)
ax.set_xticklabels([f"Node N\n($I_1 + I_5$ vs $I_4$)", f"Node M\n($I_2 + I_4$ vs $I_3$)"])
ax.set_ylim(0, 75)
ax.legend(loc='upper left', framealpha=0.9)
ax.grid(True, axis='y')

# Annotate values and residuals
for i in range(len(nodes)):
    diff = np.abs(I_out_vals[i] - I_in_vals[i])
    pct = (diff / I_out_vals[i]) * 100
    ax.text(x_indices[i] - width/2, I_in_vals[i] + 3, f"{I_in_vals[i]:.2f} mA", ha='center', fontsize=9)
    ax.text(x_indices[i] + width/2, I_out_vals[i] + 3, f"{I_out_vals[i]:.2f} mA", ha='center', fontsize=9)
    ax.text(x_indices[i], max(I_in_vals[i], I_out_vals[i]) + 8, f"$\\Delta I = {diff:.2f}$ mA\n({pct:.2f}%)",
            ha='center', fontsize=9, bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.3))

plt.tight_layout()
p2_path = plots_dir / 'kcl_node_balance.pdf'
plt.savefig(p2_path, dpi=300)
plt.close()
print(f"Saved: {p2_path}")


# Plot 3: KVL Loop Voltage Balance
fig, ax = plt.subplots(figsize=(7.5, 4.5))

loops = [
    r'Loop 1 ($V_2$ Right)' + '\n' + r'$V_2 - V_{R4} + V_{R2} - V_{R1}$',
    r'Loop 2 ($V_2$ Left)' + '\n' + r'$V_2 - V_{R4} - V_{R3} - V_{R5}$',
    r'Loop 3 ($V_1$ Top)' + '\n' + r'$V_1 - V_{R2} - V_{R3}$',
    r'Loop 4 ($V_1$ Bottom)' + '\n' + r'$V_1 - V_{R1} + V_{R5}$'
]
x_l = np.arange(len(loops))

kvl_sums = [kvl_1_sum, kvl_2_sum, kvl_3_sum, kvl_4_sum]
kvl_uncs = [kvl_1_unc, kvl_2_unc, kvl_3_unc, kvl_4_unc]
kvl_pcts = [kvl_1_pct, kvl_2_pct, kvl_3_pct, kvl_4_pct]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
bars = ax.bar(x_l, kvl_sums, yerr=kvl_uncs, capsize=5, color=colors, alpha=0.75, edgecolor='black', width=0.45)
ax.axhline(0, color='black', linestyle='-', linewidth=1.0)

ax.set_ylabel(r'Loop Potential Algebraic Sum $\sum \Delta V$ [V]')
ax.set_title("Kirchhoff's Voltage Law (KVL) Residuals Across Independent Loops")
ax.set_xticks(x_l)
ax.set_xticklabels(loops, fontsize=8.5)
ax.set_ylim(-0.15, 0.15)
ax.grid(True)

for i, b in enumerate(bars):
    y_pos = kvl_sums[i] + (0.03 if kvl_sums[i] >= 0 else -0.04)
    ax.text(b.get_x() + b.get_width()/2, y_pos,
            f"{kvl_sums[i]:+.2f} ± {kvl_uncs[i]:.2f} V\n({kvl_pcts[i]:.2f}%)",
            ha='center', va='center', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='gray', alpha=0.85))

plt.tight_layout()
p3_path = plots_dir / 'kvl_loop_balance.pdf'
plt.savefig(p3_path, dpi=300)
plt.close()
print(f"Saved: {p3_path}")


# Plot 4: Unknown Resistance Determination Comparison
fig, ax = plt.subplots(figsize=(6, 4))

methods = ["Direct Method\n(Ohm's Law: $V_{R_x} / I_{R_x}$)", "Null Bridge Method\n(Wheatstone: $R_1 R_3 / R_2$)"]
vals = [Rx1, Rx2]
errs = [delta_Rx1, delta_Rx2]
cols = ['#3498db', '#e74c3c']

ax.bar(methods, vals, yerr=errs, capsize=6, color=cols, alpha=0.8, edgecolor='black', width=0.4)
ax.set_ylabel(r'Resistance $R_x$ [$\Omega$]')
ax.set_title(r'Unknown Resistance $R_x$: Direct vs Wheatstone Bridge')
ax.set_ylim(70, 92)
ax.grid(True, axis='y')

# Add shaded compatibility interval
ax.axhspan(Rx2 - delta_Rx2, Rx2 + delta_Rx2, color='red', alpha=0.15, label='Wheatstone 1$\\sigma$ Interval')

for i in range(len(vals)):
    ax.text(i, vals[i] + errs[i] + 1.2, f"{vals[i]:.2f} ± {errs[i]:.2f} $\\Omega$", ha='center', fontsize=9.5, fontweight='bold')

ax.text(0.5, 72.5, f"Discrepancy: $\\Delta R = {delta_Rx:.2f}$ $\\Omega$ ({rel_diff_Rx:.2f}%)\nCompatibility Index $z = {z_score:.2f}\\sigma$",
        ha='center', fontsize=9, bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffffcc', edgecolor='#cccc00'))
ax.legend(loc='upper right')

plt.tight_layout()
p4_path = plots_dir / 'wheatstone_bridge_comparison.pdf'
plt.savefig(p4_path, dpi=300)
plt.close()
print(f"Saved: {p4_path}")


# Plot 5: Circuit Network Simulation vs Experimental Measurements
fig, (ax_curr, ax_volt) = plt.subplots(1, 2, figsize=(10, 4.5))

x_pos = np.arange(5)
w = 0.35

# Currents
ax_curr.bar(x_pos - w/2, I_sim_mA, w, label='Nodal Theory', color='#34495e', alpha=0.8, edgecolor='black')
ax_curr.bar(x_pos + w/2, I_exp_mA, w, yerr=delta_I_mA, capsize=3, label='Experiment', color='#3498db', alpha=0.85, edgecolor='black')
ax_curr.set_xticks(x_pos)
ax_curr.set_xticklabels(resistor_names)
ax_curr.set_ylabel('Branch Current [mA]')
ax_curr.set_title('Branch Currents: Theory vs Experiment')
ax_curr.legend()
ax_curr.grid(True, axis='y')

# Voltages
ax_volt.bar(x_pos - w/2, V_sim, w, label='Nodal Theory', color='#2c3e50', alpha=0.8, edgecolor='black')
ax_volt.bar(x_pos + w/2, V_exp, w, yerr=delta_V, capsize=3, label='Experiment', color='#e67e22', alpha=0.85, edgecolor='black')
ax_volt.set_xticks(x_pos)
ax_volt.set_xticklabels(resistor_names)
ax_volt.set_ylabel('Resistor Voltage Drop [V]')
ax_volt.set_title('Resistor Voltages: Theory vs Experiment')
ax_volt.legend()
ax_volt.grid(True, axis='y')

plt.tight_layout()
p5_path = plots_dir / 'circuit_simulation_comparison.pdf'
plt.savefig(p5_path, dpi=300)
plt.close()
print(f"Saved: {p5_path}")

print("=" * 70)
print("ALL ANALYSES AND PLOTS SUCCESSFULLY COMPLETED.")
print("=" * 70)
