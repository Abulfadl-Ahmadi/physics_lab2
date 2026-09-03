"""
Experiment 6: Alternating Current (AC) Circuits — RL, RC, and RLC Analysis
Lab-Phy-II, Department of Physics, Sharif University of Technology

Analyses:
  - Part 1: Series RL Circuit at 50 Hz (Phasor decomposition, internal resistance r_L, inductance L, impedance Z, phase phi)
  - Part 2: Series RC Circuit at 50 Hz (Orthogonality check, capacitive reactance X_C, capacitance C, impedance Z, phase phi)
  - Part 3: Series RLC Circuit at 50 Hz (Full phasor loop, re-evaluation of r_L and L, comparison of measured vs theoretical Z)
  - Resistor Linearity: OLS regression of V_R vs I across all circuits to verify Ohm's law and determine R
  - Error propagation: Analytical partial-derivative quadrature for all derived physical quantities
  - Publication-quality vector PDF plots saved to plots/
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np
from pathlib import Path

# Create plots directory
Path('plots').mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# 1. Raw Experimental Data & Measurement Uncertainties
# -----------------------------------------------------------------------------
# Operating AC frequency (Sharif University lab mains supply)
f = 50.0  # Hz
df = 0.1  # Hz uncertainty
omega = 2.0 * np.pi * f  # rad/s
domega = 2.0 * np.pi * df

# Digital Multimeter (DMM) instrument uncertainties (AC mode)
dV = 0.01      # Volts (0.01 V resolution)
dI = 0.01e-3   # Amperes (0.01 mA resolution)

# Circuit 1: Series RL
VR_1 = 1.42    # V
VL_1 = 5.03    # V
VZ_1 = 5.37    # V
I_1  = 13.54e-3 # A (13.54 mA)

# Circuit 2: Series RC
VR_2 = 2.93    # V
VC_2 = 4.46    # V
VZ_2 = 5.34    # V
I_2  = 28.29e-3 # A (28.29 mA)

# Circuit 3: Series RLC
VR_3  = 2.09    # V
VL_3  = 7.59    # V
VC_3  = 3.22    # V
VRL_3 = 8.14    # V
VZ_3  = 5.34    # V
I_3   = 20.04e-3 # A (20.04 mA)


# -----------------------------------------------------------------------------
# 2. OLS Linear Regression Helper (Lab-Phys-IV standard)
# -----------------------------------------------------------------------------
def ols_fit(x, y):
    """
    Ordinary Least Squares (OLS) linear fit: y = m*x + c
    Returns:
      m: slope
      se_m: standard error of slope
      c: intercept
      se_c: standard error of intercept
      r2: coefficient of determination R^2
    """
    N = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    s_yx = np.sqrt(ss_res / (N - 2)) if N > 2 else 0.0
    ss_x = np.sum((x - np.mean(x))**2)
    se_m = s_yx / np.sqrt(ss_x) if ss_x > 0 else 0.0
    se_c = s_yx * np.sqrt(1.0 / N + (np.mean(x)**2) / ss_x) if ss_x > 0 else 0.0
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    return m, se_m, c, se_c, r2


def draw_vector(ax, x0, y0, dx, dy, color, label=None, linestyle='-', lw=1.8, mutation_scale=14):
    """Draws a crisp phasor vector with clean arrowhead and legend handle."""
    line, = ax.plot([x0, x0 + dx], [y0, y0 + dy], color=color, linestyle=linestyle, lw=lw, label=label)
    ax.annotate('', xy=(x0 + dx, y0 + dy), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, mutation_scale=mutation_scale))
    return line


# -----------------------------------------------------------------------------
# 3. Resistor Linearity (Cross-Circuit OLS Regression)
# -----------------------------------------------------------------------------
# The same nominal resistor R was used in all three circuits
I_arr_mA = np.array([13.54, 20.04, 28.29])
I_arr_A  = I_arr_mA * 1e-3
VR_arr   = np.array([1.42, 2.09, 2.93])

slope_R, se_slope_R, intcpt_R, se_intcpt_R, r2_R = ols_fit(I_arr_A, VR_arr)


# -----------------------------------------------------------------------------
# 4. Detailed Circuit Calculations & Error Propagation
# -----------------------------------------------------------------------------

# --- PART 1: RL Circuit ---
# 1.1 Resistance of R
R_1 = VR_1 / I_1
dR_1 = R_1 * np.sqrt((dV / VR_1)**2 + (dI / I_1)**2)

# 1.2 Angle delta between V_R and V_L: V_Z^2 = V_R^2 + V_L^2 + 2*V_R*V_L*cos(delta)
# cos(delta) = (V_Z^2 - V_R^2 - V_L^2) / (2 * V_R * V_L)
cos_d1 = (VZ_1**2 - VR_1**2 - VL_1**2) / (2.0 * VR_1 * VL_1)
d1 = np.arccos(cos_d1)
# Error propagation on cos_d1:
dcos_d1_dVZ = VZ_1 / (VR_1 * VL_1)
dcos_d1_dVR = -(VZ_1**2 + VR_1**2 - VL_1**2) / (2.0 * VR_1**2 * VL_1)
dcos_d1_dVL = -(VZ_1**2 + VL_1**2 - VR_1**2) / (2.0 * VR_1 * VL_1**2)
dcos_d1 = np.sqrt((dcos_d1_dVZ * dV)**2 + (dcos_d1_dVR * dV)**2 + (dcos_d1_dVL * dV)**2)
dd1 = dcos_d1 / np.sin(d1)  # radians

# 1.3 In-phase and quadrature components of V_L
VrL_1 = VL_1 * cos_d1
# VrL = (VZ^2 - VR^2 - VL^2) / (2 * VR)
dVrL_dVZ = VZ_1 / VR_1
dVrL_dVR = -(VZ_1**2 + VR_1**2 - VL_1**2) / (2.0 * VR_1**2)
dVrL_dVL = -VL_1 / VR_1
dVrL_1 = np.sqrt((dVrL_dVZ * dV)**2 + (dVrL_dVR * dV)**2 + (dVrL_dVL * dV)**2)

VXL_1 = VL_1 * np.sin(d1)
dVXL_1 = np.sqrt((VL_1 * dV / VXL_1)**2 + (VrL_1 * dVrL_1 / VXL_1)**2)

# 1.4 Inductor internal resistance r_L and inductive reactance X_L
rL_1 = VrL_1 / I_1
drL_1 = rL_1 * np.sqrt((dVrL_1 / VrL_1)**2 + (dI / I_1)**2)

XL_1 = VXL_1 / I_1
dXL_1 = XL_1 * np.sqrt((dVXL_1 / VXL_1)**2 + (dI / I_1)**2)

# 1.5 Inductance L
L_1 = XL_1 / omega
dL_1 = L_1 * np.sqrt((dXL_1 / XL_1)**2 + (df / f)**2)

# 1.6 Total resistance and impedances
Rtot_1 = R_1 + rL_1
dRtot_1 = np.sqrt(dR_1**2 + drL_1**2)

Z_1_meas = VZ_1 / I_1
dZ_1_meas = Z_1_meas * np.sqrt((dV / VZ_1)**2 + (dI / I_1)**2)

Z_1_th = np.sqrt(Rtot_1**2 + XL_1**2)
dZ_1_th = np.sqrt((Rtot_1 * dRtot_1 / Z_1_th)**2 + (XL_1 * dXL_1 / Z_1_th)**2)

# 1.7 Phase angle phi
phi_1_rad = np.arctan2(XL_1, Rtot_1)
phi_1_deg = np.degrees(phi_1_rad)
dphi_1_rad = np.sqrt((Rtot_1 * dXL_1)**2 + (XL_1 * dRtot_1)**2) / (Z_1_th**2)
dphi_1_deg = np.degrees(dphi_1_rad)


# --- PART 2: RC Circuit ---
# 2.1 Resistance of R
R_2 = VR_2 / I_2
dR_2 = R_2 * np.sqrt((dV / VR_2)**2 + (dI / I_2)**2)

# 2.2 Capacitive reactance X_C and Capacitance C
XC_2 = VC_2 / I_2
dXC_2 = XC_2 * np.sqrt((dV / VC_2)**2 + (dI / I_2)**2)

C_2 = 1.0 / (omega * XC_2)
dC_2 = C_2 * np.sqrt((dXC_2 / XC_2)**2 + (df / f)**2)

# 2.3 Orthogonality check between VR and VC: cos(delta_RC)
cos_d2 = (VZ_2**2 - VR_2**2 - VC_2**2) / (2.0 * VR_2 * VC_2)
d2 = np.arccos(cos_d2)
dcos_d2_dVZ = VZ_2 / (VR_2 * VC_2)
dcos_d2_dVR = -(VZ_2**2 + VR_2**2 - VC_2**2) / (2.0 * VR_2**2 * VC_2)
dcos_d2_dVC = -(VZ_2**2 + VC_2**2 - VR_2**2) / (2.0 * VR_2 * VC_2**2)
dcos_d2 = np.sqrt((dcos_d2_dVZ * dV)**2 + (dcos_d2_dVR * dV)**2 + (dcos_d2_dVC * dV)**2)
dd2 = dcos_d2 / np.sin(d2)

# 2.4 Measured and theoretical impedance
Z_2_meas = VZ_2 / I_2
dZ_2_meas = Z_2_meas * np.sqrt((dV / VZ_2)**2 + (dI / I_2)**2)

Z_2_th = np.sqrt(R_2**2 + XC_2**2)
dZ_2_th = np.sqrt((R_2 * dR_2 / Z_2_th)**2 + (XC_2 * dXC_2 / Z_2_th)**2)

# 2.5 Phase angle (voltage lags current, phi < 0)
phi_2_rad = -np.arctan2(XC_2, R_2)
phi_2_deg = np.degrees(phi_2_rad)
dphi_2_rad = np.sqrt((R_2 * dXC_2)**2 + (XC_2 * dR_2)**2) / (Z_2_th**2)
dphi_2_deg = np.degrees(dphi_2_rad)


# --- PART 3: RLC Circuit ---
# 3.1 Resistance of R
R_3 = VR_3 / I_3
dR_3 = R_3 * np.sqrt((dV / VR_3)**2 + (dI / I_3)**2)

# 3.2 Inductor internal angle from V_RL
cos_d3 = (VRL_3**2 - VR_3**2 - VL_3**2) / (2.0 * VR_3 * VL_3)
d3 = np.arccos(cos_d3)
dcos_d3_dVRL = VRL_3 / (VR_3 * VL_3)
dcos_d3_dVR  = -(VRL_3**2 + VR_3**2 - VL_3**2) / (2.0 * VR_3**2 * VL_3)
dcos_d3_dVL  = -(VRL_3**2 + VL_3**2 - VR_3**2) / (2.0 * VR_3 * VL_3**2)
dcos_d3 = np.sqrt((dcos_d3_dVRL * dV)**2 + (dcos_d3_dVR * dV)**2 + (dcos_d3_dVL * dV)**2)
dd3 = dcos_d3 / np.sin(d3)

# 3.3 Inductor components in RLC
VrL_3 = VL_3 * cos_d3
dVrL_3_dVRL = VRL_3 / VR_3
dVrL_3_dVR  = -(VRL_3**2 + VR_3**2 - VL_3**2) / (2.0 * VR_3**2)
dVrL_3_dVL  = -VL_3 / VR_3
dVrL_3 = np.sqrt((dVrL_3_dVRL * dV)**2 + (dVrL_3_dVR * dV)**2 + (dVrL_3_dVL * dV)**2)

VXL_3 = VL_3 * np.sin(d3)
dVXL_3 = np.sqrt((VL_3 * dV / VXL_3)**2 + (VrL_3 * dVrL_3 / VXL_3)**2)

rL_3 = VrL_3 / I_3
drL_3 = rL_3 * np.sqrt((dVrL_3 / VrL_3)**2 + (dI / I_3)**2)

XL_3 = VXL_3 / I_3
dXL_3 = XL_3 * np.sqrt((dVXL_3 / VXL_3)**2 + (dI / I_3)**2)

L_3 = XL_3 / omega
dL_3 = L_3 * np.sqrt((dXL_3 / XL_3)**2 + (df / f)**2)

# 3.4 Capacitor in RLC
XC_3 = VC_3 / I_3
dXC_3 = XC_3 * np.sqrt((dV / VC_3)**2 + (dI / I_3)**2)

C_3 = 1.0 / (omega * XC_3)
dC_3 = C_3 * np.sqrt((dXC_3 / XC_3)**2 + (df / f)**2)

# 3.5 Total resistance, net reactance, impedance
Rtot_3 = R_3 + rL_3
dRtot_3 = np.sqrt(dR_3**2 + drL_3**2)

Xnet_3 = XL_3 - XC_3
dXnet_3 = np.sqrt(dXL_3**2 + dXC_3**2)

Z_3_meas = VZ_3 / I_3
dZ_3_meas = Z_3_meas * np.sqrt((dV / VZ_3)**2 + (dI / I_3)**2)

Z_3_th = np.sqrt(Rtot_3**2 + Xnet_3**2)
dZ_3_th = np.sqrt((Rtot_3 * dRtot_3 / Z_3_th)**2 + (Xnet_3 * dXnet_3 / Z_3_th)**2)

# Cross-prediction of Z_RLC using parameters from RL (rL, L) and RC (C)
Rtot_cross = R_2 + rL_1
dRtot_cross = np.sqrt(dR_2**2 + drL_1**2)
Xnet_cross = XL_1 - XC_2
dXnet_cross = np.sqrt(dXL_1**2 + dXC_2**2)
Z_3_cross = np.sqrt(Rtot_cross**2 + Xnet_cross**2)
dZ_3_cross = np.sqrt((Rtot_cross * dRtot_cross / Z_3_cross)**2 + (Xnet_cross * dXnet_cross / Z_3_cross)**2)

# 3.6 Phase angle
phi_3_rad = np.arctan2(Xnet_3, Rtot_3)
phi_3_deg = np.degrees(phi_3_rad)
dphi_3_rad = np.sqrt((Rtot_3 * dXnet_3)**2 + (Xnet_3 * dRtot_3)**2) / (Z_3_th**2)
dphi_3_deg = np.degrees(dphi_3_rad)

# 3.7 Resonance frequency
f0_1 = 1.0 / (2.0 * np.pi * np.sqrt(L_1 * C_2))
df0_1 = f0_1 * 0.5 * np.sqrt((dL_1 / L_1)**2 + (dC_2 / C_2)**2)

f0_3 = 1.0 / (2.0 * np.pi * np.sqrt(L_3 * C_3))
df0_3 = f0_3 * 0.5 * np.sqrt((dL_3 / L_3)**2 + (dC_3 / C_3)**2)


# -----------------------------------------------------------------------------
# 5. Formatted Console Output
# -----------------------------------------------------------------------------
def print_header(title):
    print("=" * 72)
    print(f" {title}")
    print("=" * 72)

print_header("EXPERIMENT 6: AC CIRCUITS (RL, RC, RLC) ANALYSIS REPORT")
print(f"Supply frequency: f = {f:.1f} +/- {df:.1f} Hz (omega = {omega:.2f} rad/s)")
print()

print_header("1. RESISTOR LINEARITY FIT (Ohm's Law: V_R = R * I + offset)")
print(f"  Slope (Resistor R) : {slope_R:.3f} +/- {se_slope_R:.3f} Ohm")
print(f"  Intercept (Offset) : {intcpt_R*1e3:.3f} +/- {se_intcpt_R*1e3:.3f} mV")
print(f"  Coefficient R^2    : {r2_R:.6f}")
print()

print_header("2. SERIES RL CIRCUIT")
print(f"  Resistor R             : {R_1:.2f} +/- {dR_1:.2f} Ohm")
print(f"  Coil angle delta       : {np.degrees(d1):.2f} +/- {np.degrees(dd1):.2f} deg  (cos(delta) = {cos_d1:.4f} +/- {dcos_d1:.4f})")
print(f"  Coil resistive drop VrL: {VrL_1:.3f} +/- {dVrL_1:.3f} V")
print(f"  Coil reactive drop VXL : {VXL_1:.3f} +/- {dVXL_1:.3f} V")
print(f"  Coil int. resistance rL: {rL_1:.2f} +/- {drL_1:.2f} Ohm")
print(f"  Inductive reactance XL : {XL_1:.2f} +/- {dXL_1:.2f} Ohm")
print(f"  Inductance L           : {L_1:.4f} +/- {dL_1:.4f} H")
print(f"  Total resistance R_tot : {Rtot_1:.2f} +/- {dRtot_1:.2f} Ohm")
print(f"  Measured impedance |Z| : {Z_1_meas:.2f} +/- {dZ_1_meas:.2f} Ohm")
print(f"  Calculated imped. |Z|  : {Z_1_th:.2f} +/- {dZ_1_th:.2f} Ohm")
print(f"  Phase angle phi (I->V) : {phi_1_deg:.2f} +/- {dphi_1_deg:.2f} deg (inductive, V leads I)")
print()

print_header("3. SERIES RC CIRCUIT")
print(f"  Resistor R             : {R_2:.2f} +/- {dR_2:.2f} Ohm")
print(f"  Capacitive reactance XC: {XC_2:.2f} +/- {dXC_2:.2f} Ohm")
print(f"  Capacitance C          : {C_2*1e6:.2f} +/- {dC_2*1e6:.2f} uF")
print(f"  Orthogonality delta_RC : {np.degrees(d2):.2f} +/- {np.degrees(dd2):.2f} deg  (cos(delta) = {cos_d2:.5f})")
print(f"  Measured impedance |Z| : {Z_2_meas:.2f} +/- {dZ_2_meas:.2f} Ohm")
print(f"  Calculated imped. |Z|  : {Z_2_th:.2f} +/- {dZ_2_th:.2f} Ohm")
print(f"  Phase angle phi (I->V) : {phi_2_deg:.2f} +/- {dphi_2_deg:.2f} deg (capacitive, I leads V)")
print()

print_header("4. SERIES RLC CIRCUIT")
print(f"  Resistor R             : {R_3:.2f} +/- {dR_3:.2f} Ohm")
print(f"  Coil angle delta (RLC) : {np.degrees(d3):.2f} +/- {np.degrees(dd3):.2f} deg  (cos(delta) = {cos_d3:.4f})")
print(f"  Coil int. resistance rL: {rL_3:.2f} +/- {drL_3:.2f} Ohm")
print(f"  Inductive reactance XL : {XL_3:.2f} +/- {dXL_3:.2f} Ohm")
print(f"  Inductance L           : {L_3:.4f} +/- {dL_3:.4f} H")
print(f"  Capacitive reactance XC: {XC_3:.2f} +/- {dXC_3:.2f} Ohm")
print(f"  Capacitance C          : {C_3*1e6:.2f} +/- {dC_3*1e6:.2f} uF")
print(f"  Total resistance R_tot : {Rtot_3:.2f} +/- {dRtot_3:.2f} Ohm")
print(f"  Net reactance (XL - XC): {Xnet_3:.2f} +/- {dXnet_3:.2f} Ohm")
print(f"  Measured impedance |Z| : {Z_3_meas:.2f} +/- {dZ_3_meas:.2f} Ohm")
print(f"  Internal theoretical Z : {Z_3_th:.2f} +/- {dZ_3_th:.2f} Ohm")
print(f"  Cross-predicted Z      : {Z_3_cross:.2f} +/- {dZ_3_cross:.2f} Ohm  (diff: {abs(Z_3_meas - Z_3_cross)/Z_3_meas*100:.2f}%)")
print(f"  Phase angle phi (I->V) : {phi_3_deg:.2f} +/- {dphi_3_deg:.2f} deg (inductive, XL > XC)")
print(f"  Resonance frequency f0 : {f0_1:.2f} +/- {df0_1:.2f} Hz  (f = 50 Hz > f0 => inductive behavior)")
print()


# -----------------------------------------------------------------------------
# 6. Plot Generation
# -----------------------------------------------------------------------------
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['mathtext.fontset'] = 'cm'

# -----------------------------------------------------------------------------
# Plot 1: Resistor Linearity OLS Fit
# -----------------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(7, 5), dpi=300)
ax1.errorbar(I_arr_mA, VR_arr, xerr=dI*1e3, yerr=dV, fmt='o', color='#1f77b4',
             ecolor='black', elinewidth=1.2, capsize=4, capthick=1.2,
             label='Measured Data (RL, RLC, RC)', zorder=5)

I_line_mA = np.linspace(10, 32, 100)
I_line_A  = I_line_mA * 1e-3
VR_line   = slope_R * I_line_A + intcpt_R
ax1.plot(I_line_mA, VR_line, 'r--', lw=1.8,
         label=f'OLS Fit: $V_R = R \\cdot I + V_0$\n'
               f'$R = ({slope_R:.2f} \\pm {se_slope_R:.2f})\\ \\Omega$\n'
               f'$V_0 = ({intcpt_R*1e3:.1f} \\pm {se_intcpt_R*1e3:.1f})\\ \\mathrm{{mV}}$\n'
               f'$R^2 = {r2_R:.5f}$', zorder=4)

ax1.set_xlabel('Current $I$ [mA]', fontsize=12, labelpad=8)
ax1.set_ylabel('Resistor Voltage $V_R$ [V]', fontsize=12, labelpad=8)
ax1.set_title("Ohm's Law Verification across AC Circuits (RL, RC, RLC)", fontsize=13, pad=12, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='upper left', fontsize=10, frameon=True, framealpha=0.95)
ax1.set_xlim(10, 32)
ax1.set_ylim(1.0, 3.3)
fig1.tight_layout()
fig1.savefig('plots/VR_vs_I_resistor.pdf')
plt.close(fig1)
print("Saved: plots/VR_vs_I_resistor.pdf")


# -----------------------------------------------------------------------------
# Plot 2: Detailed Phasor Diagram — RL Circuit
# -----------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7, 6), dpi=300)
x_VR, y_VR = VR_1, 0.0
x_VrL, y_VrL = VR_1 + VrL_1, 0.0
x_top, y_top = VR_1 + VrL_1, VXL_1

draw_vector(ax2, 0, 0, x_VR, 0, '#2ca02c', label=f'$V_R = {VR_1:.2f}\\ \\mathrm{{V}}$', lw=2.2)
draw_vector(ax2, x_VR, y_VR, VrL_1, 0, '#bcbd22', linestyle='--',
            label=f'$V_{{r_L}} = {VrL_1:.2f}\\ \\mathrm{{V}}$ (Internal resistance)', lw=1.8)
draw_vector(ax2, x_VrL, y_VrL, 0, VXL_1, '#9467bd', linestyle='--',
            label=f'$V_{{X_L}} = {VXL_1:.2f}\\ \\mathrm{{V}}$ (Inductive reactance)', lw=1.8)
draw_vector(ax2, x_VR, y_VR, VrL_1, VXL_1, '#d62728',
            label=f'$V_L = {VL_1:.2f}\\ \\mathrm{{V}}$ (Total coil)', lw=2.2)
draw_vector(ax2, 0, 0, x_top, y_top, '#1f77b4',
            label=f'$V_Z = {VZ_1:.2f}\\ \\mathrm{{V}}$ (Total input)', lw=2.5)

arc_phi = Arc((0, 0), 1.6, 1.6, angle=0, theta1=0, theta2=phi_1_deg, color='black', lw=1.2)
ax2.add_patch(arc_phi)
ax2.text(1.0, 0.5, f'$\\varphi = {phi_1_deg:.1f}^\\circ$', fontsize=11, fontweight='bold')

arc_delta = Arc((VR_1, 0), 1.2, 1.2, angle=0, theta1=0, theta2=np.degrees(d1), color='#8c564b', lw=1.2)
ax2.add_patch(arc_delta)
ax2.text(VR_1 + 0.2, 0.8, f'$\\delta = {np.degrees(d1):.1f}^\\circ$', fontsize=10, color='#8c564b')

ax2.set_xlim(-0.5, 3.5)
ax2.set_ylim(-0.5, 5.8)
ax2.set_xlabel('In-phase Component (along Current $\\vec{I}$) [V]', fontsize=11)
ax2.set_ylabel('Quadrature Component ($+90^\\circ$ leading) [V]', fontsize=11)
ax2.set_title('Phasor Diagram — Series RL Circuit ($f = 50\\ \\mathrm{Hz}$)', fontsize=13, pad=12, fontweight='bold')
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper left', fontsize=9, framealpha=0.95)
ax2.set_aspect('equal')
fig2.tight_layout()
fig2.savefig('plots/phasor_RL.pdf')
plt.close(fig2)
print("Saved: plots/phasor_RL.pdf")


# -----------------------------------------------------------------------------
# Plot 3: Detailed Phasor Diagram — RC Circuit
# -----------------------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(7, 6), dpi=300)
draw_vector(ax3, 0, 0, VR_2, 0, '#2ca02c', label=f'$V_R = {VR_2:.2f}\\ \\mathrm{{V}}$', lw=2.2)
draw_vector(ax3, VR_2, 0, 0, -VC_2, '#ff7f0e', label=f'$V_C = {VC_2:.2f}\\ \\mathrm{{V}}$ (Capacitor)', lw=2.2)
draw_vector(ax3, 0, 0, VR_2, -VC_2, '#1f77b4', label=f'$V_Z = {VZ_2:.2f}\\ \\mathrm{{V}}$ (Total input)', lw=2.5)

arc_phi2 = Arc((0, 0), 1.8, 1.8, angle=0, theta1=phi_2_deg, theta2=0, color='black', lw=1.2)
ax3.add_patch(arc_phi2)
ax3.text(1.1, -0.6, f'$|\\varphi| = {abs(phi_2_deg):.1f}^\\circ$', fontsize=11, fontweight='bold')

rect_size = 0.25
ax3.plot([VR_2 - rect_size, VR_2 - rect_size, VR_2],
         [0, -rect_size, -rect_size], color='gray', lw=1.2)

ax3.set_xlim(-0.5, 4.0)
ax3.set_ylim(-5.2, 0.8)
ax3.set_xlabel('In-phase Component (along Current $\\vec{I}$) [V]', fontsize=11)
ax3.set_ylabel('Quadrature Component ($-90^\\circ$ lagging) [V]', fontsize=11)
ax3.set_title('Phasor Diagram — Series RC Circuit ($f = 50\\ \\mathrm{Hz}$)', fontsize=13, pad=12, fontweight='bold')
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.legend(loc='lower left', fontsize=9, framealpha=0.95)
ax3.set_aspect('equal')
fig3.tight_layout()
fig3.savefig('plots/phasor_RC.pdf')
plt.close(fig3)
print("Saved: plots/phasor_RC.pdf")


# -----------------------------------------------------------------------------
# Plot 4: Detailed Phasor Diagram — RLC Circuit
# -----------------------------------------------------------------------------
fig4, ax4 = plt.subplots(figsize=(8, 7), dpi=300)
x_VR3 = VR_3
x_Rtot3 = VR_3 + VrL_3
y_net3 = VXL_3 - VC_3

draw_vector(ax4, 0, 0, x_VR3, 0, '#2ca02c', label=f'$V_R = {VR_3:.2f}\\ \\mathrm{{V}}$', lw=2.0)
draw_vector(ax4, x_VR3, 0, 0, -VC_3, '#ff7f0e', label=f'$V_C = {VC_3:.2f}\\ \\mathrm{{V}}$ (Downwards)', lw=2.0)
draw_vector(ax4, x_VR3, -VC_3, VrL_3, VXL_3, '#d62728', label=f'$V_L = {VL_3:.2f}\\ \\mathrm{{V}}$ (Coil vector)', lw=2.2)
draw_vector(ax4, 0, 0, x_Rtot3, VXL_3, '#8c564b', linestyle=':', label=f'$V_{{RL}} = {VRL_3:.2f}\\ \\mathrm{{V}}$ (R + Coil)', lw=1.8)
draw_vector(ax4, 0, 0, x_Rtot3, y_net3, '#1f77b4', label=f'$V_Z = {VZ_3:.2f}\\ \\mathrm{{V}}$ (Total input)', lw=2.5)

# Net reactive line
ax4.plot([x_Rtot3, x_Rtot3], [0, y_net3], 'k--', lw=1.4, alpha=0.7,
         label=f'$V_{{net,X}} = V_{{X_L}} - V_C = {y_net3:.2f}\\ \\mathrm{{V}}$')

arc_phi3 = Arc((0, 0), 2.0, 2.0, angle=0, theta1=0, theta2=phi_3_deg, color='black', lw=1.2)
ax4.add_patch(arc_phi3)
ax4.text(1.3, 0.7, f'$\\varphi = {phi_3_deg:.1f}^\\circ$', fontsize=11, fontweight='bold')

ax4.set_xlim(-0.5, 4.5)
ax4.set_ylim(-4.0, 8.5)
ax4.set_xlabel('In-phase Component [V]', fontsize=11)
ax4.set_ylabel('Quadrature Component [V]', fontsize=11)
ax4.set_title('Phasor Diagram — Series RLC Circuit ($f = 50\\ \\mathrm{Hz}$)', fontsize=13, pad=12, fontweight='bold')
ax4.grid(True, linestyle=':', alpha=0.6)
ax4.legend(loc='upper left', fontsize=8.5, framealpha=0.95)
ax4.set_aspect('equal')
fig4.tight_layout()
fig4.savefig('plots/phasor_RLC.pdf')
plt.close(fig4)
print("Saved: plots/phasor_RLC.pdf")


# -----------------------------------------------------------------------------
# Plot 5: 3-Panel Comprehensive Phasor Comparison
# -----------------------------------------------------------------------------
fig5, (p_rl, p_rc, p_rlc) = plt.subplots(1, 3, figsize=(15, 5), dpi=300)

# Panel 1: RL
draw_vector(p_rl, 0, 0, VR_1, 0, '#2ca02c', label='$V_R$', lw=2.0)
draw_vector(p_rl, VR_1, 0, VrL_1, VXL_1, '#d62728', label='$V_L$', lw=2.0)
draw_vector(p_rl, 0, 0, VR_1 + VrL_1, VXL_1, '#1f77b4', label='$V_Z$', lw=2.2)
p_rl.add_patch(Arc((0, 0), 1.6, 1.6, angle=0, theta1=0, theta2=phi_1_deg, color='black', lw=1.2))
p_rl.text(0.9, 0.5, f'$\\varphi={phi_1_deg:.1f}^\\circ$', fontsize=10)
p_rl.set_xlim(-0.5, 3.0)
p_rl.set_ylim(-0.5, 5.5)
p_rl.set_title(f'RL Circuit ($\\varphi = {phi_1_deg:.1f}^\\circ$)', fontsize=11, fontweight='bold')
p_rl.set_xlabel('In-phase [V]', fontsize=10)
p_rl.set_ylabel('Quadrature [V]', fontsize=10)
p_rl.grid(True, linestyle=':', alpha=0.5)
p_rl.legend(loc='upper left', fontsize=9)
p_rl.set_aspect('equal')

# Panel 2: RC
draw_vector(p_rc, 0, 0, VR_2, 0, '#2ca02c', label='$V_R$', lw=2.0)
draw_vector(p_rc, VR_2, 0, 0, -VC_2, '#ff7f0e', label='$V_C$', lw=2.0)
draw_vector(p_rc, 0, 0, VR_2, -VC_2, '#1f77b4', label='$V_Z$', lw=2.2)
p_rc.add_patch(Arc((0, 0), 1.8, 1.8, angle=0, theta1=phi_2_deg, theta2=0, color='black', lw=1.2))
p_rc.text(1.1, -0.6, f'$|\\varphi|={abs(phi_2_deg):.1f}^\\circ$', fontsize=10)
p_rc.set_xlim(-0.5, 3.8)
p_rc.set_ylim(-5.2, 0.8)
p_rc.set_title(f'RC Circuit ($\\varphi = {phi_2_deg:.1f}^\\circ$)', fontsize=11, fontweight='bold')
p_rc.set_xlabel('In-phase [V]', fontsize=10)
p_rc.grid(True, linestyle=':', alpha=0.5)
p_rc.legend(loc='lower left', fontsize=9)
p_rc.set_aspect('equal')

# Panel 3: RLC
draw_vector(p_rlc, 0, 0, VR_3, 0, '#2ca02c', label='$V_R$', lw=2.0)
draw_vector(p_rlc, VR_3, 0, 0, -VC_3, '#ff7f0e', label='$V_C$', lw=2.0)
draw_vector(p_rlc, VR_3, -VC_3, VrL_3, VXL_3, '#d62728', label='$V_L$', lw=2.0)
draw_vector(p_rlc, 0, 0, VR_3 + VrL_3, VXL_3 - VC_3, '#1f77b4', label='$V_Z$', lw=2.2)
p_rlc.add_patch(Arc((0, 0), 2.0, 2.0, angle=0, theta1=0, theta2=phi_3_deg, color='black', lw=1.2))
p_rlc.text(1.2, 0.6, f'$\\varphi={phi_3_deg:.1f}^\\circ$', fontsize=10)
p_rlc.set_xlim(-0.5, 4.0)
p_rlc.set_ylim(-3.8, 5.5)
p_rlc.set_title(f'RLC Circuit ($\\varphi = {phi_3_deg:.1f}^\\circ$)', fontsize=11, fontweight='bold')
p_rlc.set_xlabel('In-phase [V]', fontsize=10)
p_rlc.grid(True, linestyle=':', alpha=0.5)
p_rlc.legend(loc='upper left', fontsize=9)
p_rlc.set_aspect('equal')

fig5.suptitle('Comparison of AC Voltage Phasor Diagrams at 50 Hz', fontsize=14, fontweight='bold', y=0.98)
fig5.tight_layout()
fig5.savefig('plots/phasor_summary.pdf')
plt.close(fig5)
print("Saved: plots/phasor_summary.pdf")


# -----------------------------------------------------------------------------
# Plot 6: Impedance Comparison (Measured vs Theoretical)
# -----------------------------------------------------------------------------
fig6, ax6 = plt.subplots(figsize=(7, 5), dpi=300)
categories = ['RL Circuit', 'RC Circuit', 'RLC Circuit']
Z_meas_vals = [Z_1_meas, Z_2_meas, Z_3_meas]
dZ_meas_vals = [dZ_1_meas, dZ_2_meas, dZ_3_meas]
Z_th_vals = [Z_1_th, Z_2_th, Z_3_cross]
dZ_th_vals = [dZ_1_th, dZ_2_th, dZ_3_cross]

x_pos = np.arange(len(categories))
bar_width = 0.35

ax6.bar(x_pos - bar_width/2, Z_meas_vals, bar_width, yerr=dZ_meas_vals,
        capsize=5, color='#1f77b4', edgecolor='black', alpha=0.85, label='Measured $|Z| = V_Z / I$')
ax6.bar(x_pos + bar_width/2, Z_th_vals, bar_width, yerr=dZ_th_vals,
        capsize=5, color='#ff7f0e', edgecolor='black', alpha=0.85, label='Theoretical $|Z|$ (with $r_L$)')

for i in range(len(categories)):
    diff_pct = abs(Z_meas_vals[i] - Z_th_vals[i]) / Z_meas_vals[i] * 100
    ax6.text(x_pos[i], max(Z_meas_vals[i], Z_th_vals[i]) + 20,
             f'$\\Delta = {diff_pct:.1f}\\%$', ha='center', fontsize=9, fontweight='bold')

ax6.set_ylabel('Impedance $|Z|\\ [\\Omega]$', fontsize=12, labelpad=8)
ax6.set_title('Comparison of Measured and Theoretical Impedances', fontsize=13, pad=12, fontweight='bold')
ax6.set_xticks(x_pos)
ax6.set_xticklabels(categories, fontsize=11)
ax6.set_ylim(0, 470)
ax6.grid(True, axis='y', linestyle=':', alpha=0.6)
ax6.legend(loc='upper right', fontsize=10, framealpha=0.95)
fig6.tight_layout()
fig6.savefig('plots/impedance_comparison.pdf')
plt.close(fig6)
print("Saved: plots/impedance_comparison.pdf")


# -----------------------------------------------------------------------------
# Plot 7: RLC Frequency Response & Resonance Curve
# -----------------------------------------------------------------------------
fig7, (ax7a, ax7b) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, dpi=300)

f_range = np.linspace(10, 100, 500)
omega_range = 2.0 * np.pi * f_range

# Use mean component values from experiment:
R_eff = R_2 + rL_1  # ~143.1 Ohm
L_eff = L_1          # ~1.176 H
C_eff = C_2          # ~20.2 uF
V_source = 5.34      # V

XL_f = omega_range * L_eff
XC_f = 1.0 / (omega_range * C_eff)
Z_f = np.sqrt(R_eff**2 + (XL_f - XC_f)**2)
I_f_mA = (V_source / Z_f) * 1e3
phi_f_deg = np.degrees(np.arctan2(XL_f - XC_f, R_eff))

# Subplot A: Current vs Frequency (Resonance Peak)
ax7a.plot(f_range, I_f_mA, 'b-', lw=2, label='Series Current $I(f)$')
ax7a.axvline(f0_1, color='red', linestyle='--', lw=1.5,
             label=f'Resonance $f_0 = {f0_1:.1f}\\ \\mathrm{{Hz}}$ ($I_{{\\max}} = {V_source/R_eff*1e3:.1f}\\ \\mathrm{{mA}}$)')
ax7a.axvline(50.0, color='darkgreen', linestyle=':', lw=1.5,
             label=f'Operating Point $f = 50\\ \\mathrm{{Hz}}$ ($I = {I_3*1e3:.1f}\\ \\mathrm{{mA}}$)')
ax7a.scatter([50.0], [I_3*1e3], color='darkgreen', s=60, zorder=5)
ax7a.scatter([f0_1], [V_source/R_eff*1e3], color='red', s=60, zorder=5)

ax7a.set_ylabel('Current $I$ [mA]', fontsize=11)
ax7a.set_title('Series RLC Circuit Frequency Response & Resonance Analysis', fontsize=13, pad=12, fontweight='bold')
ax7a.grid(True, linestyle=':', alpha=0.6)
ax7a.legend(loc='upper right', fontsize=9.5, framealpha=0.95)

# Subplot B: Impedance & Phase vs Frequency
ax7b.plot(f_range, Z_f, 'k-', lw=1.8, label='Impedance $|Z|(f)$')
ax7b.axvline(f0_1, color='red', linestyle='--', lw=1.5)
ax7b.axvline(50.0, color='darkgreen', linestyle=':', lw=1.5)
ax7b.set_xlabel('Frequency $f$ [Hz]', fontsize=11)
ax7b.set_ylabel('Impedance $|Z|\\ [\\Omega]$', fontsize=11)
ax7b.grid(True, linestyle=':', alpha=0.6)

# Add twin axis for phase angle
ax7b_twin = ax7b.twinx()
ax7b_twin.plot(f_range, phi_f_deg, 'purple', linestyle='-.', lw=1.6, label='Phase Angle $\\varphi(f)$')
ax7b_twin.axhline(0, color='gray', linestyle='-', lw=0.8, alpha=0.7)
ax7b_twin.set_ylabel('Phase $\\varphi\\ [^\\circ]$', color='purple', fontsize=11)
ax7b_twin.tick_params(axis='y', labelcolor='purple')

# Combine legends for lower plot
lines_b, labels_b = ax7b.get_legend_handles_labels()
lines_twin, labels_twin = ax7b_twin.get_legend_handles_labels()
ax7b.legend(lines_b + lines_twin, labels_b + labels_twin, loc='center right', fontsize=9.5, framealpha=0.95)

fig7.tight_layout()
fig7.savefig('plots/resonance_curve.pdf')
plt.close(fig7)
print("Saved: plots/resonance_curve.pdf")

print_header("ALL PLOTS GENERATED AND SAVED TO plots/ SUCCESSFULLY")
