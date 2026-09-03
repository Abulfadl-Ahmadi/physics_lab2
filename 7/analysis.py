"""
Experiment 7: Capacitor Characteristics & Dielectric Constant Measurement
              (with Oscilloscope Signal Analysis)
Sharif University of Technology - Department of Physics
Lab-Phy-II

This script performs:
1. Ordinary Least Squares (OLS) regression for I vs f (Plexiglass) -> C_plexi, K_plexi
2. OLS regression for I vs f (Air) -> C_air, epsilon_air, K_air
3. Geometry-independent ratio C_plexi / C_air
4. OLS regression for I vs 1/d (Distance dependence) -> slope, intercept, stray capacitance
5. Oscilloscope frequency error propagation and calibration analysis
6. High-quality publication-ready plots saved to plots/ directory
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Set plotting style
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'errorbar.capsize': 3,
    'grid.linestyle': '--',
    'grid.alpha': 0.6
})

# ---------------------------------------------------------
# OLS Fit Helper Function (matching Lab-Phys-IV standard)
# ---------------------------------------------------------
def ols_fit(x, y):
    """
    Perform Ordinary Least Squares linear regression y = m*x + c.
    Returns: slope (m), se_slope, intercept (c), se_intercept, r_squared
    """
    n = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 1.0
    
    s_yx = np.sqrt(ss_res / (n - 2))
    ss_x = np.sum((x - np.mean(x))**2)
    se_m = s_yx / np.sqrt(ss_x)
    se_c = s_yx * np.sqrt(1.0 / n + (np.mean(x)**2) / ss_x)
    
    return m, se_m, c, se_c, r2

# ---------------------------------------------------------
# Physical Constants & Apparatus Geometry
# ---------------------------------------------------------
epsilon_0 = 8.8541878128e-12  # F/m (Permittivity of vacuum)
r_plate = 0.10                # m (radius = 10 cm)
delta_r = 0.001               # m (1 mm uncertainty in radius)
A_plate = np.pi * r_plate**2   # m^2 (Area)
delta_A = 2 * A_plate * (delta_r / r_plate)

d_fixed = 2.8e-3              # m (2.8 mm spacing for Plexiglass & Air)
delta_d = 0.05e-3             # m (0.05 mm uncertainty)

# Voltage settings
Vm = 4.0                      # V (Peak amplitude)
delta_Vm = 0.05               # V
V_rms = Vm / np.sqrt(2)       # V_rms = 2.8284 V
delta_V_rms = delta_Vm / np.sqrt(2)

print("="*65)
print("EXPERIMENT 7: CAPACITANCE & DIELECTRIC CONSTANT MEASUREMENT")
print("="*65)
print(f"Plate radius: r = {r_plate*100:.1f} +/- {delta_r*100:.1f} cm")
print(f"Plate area:   A = {A_plate*1e4:.2f} +/- {delta_A*1e4:.2f} cm^2")
print(f"Plate gap:    d = {d_fixed*1e3:.2f} +/- {delta_d*1e3:.2f} mm")
print(f"Peak voltage: Vm = {Vm:.2f} V -> V_rms = {V_rms:.4f} V")
print()

# Ensure plots/ directory exists
plots_dir = Path('plots')
plots_dir.mkdir(exist_ok=True)

# ---------------------------------------------------------
# Part 1: Plexiglass Dielectric (Table 1 from 7.tex)
# ---------------------------------------------------------
f_plexi_kHz = np.array([1.0, 5.0, 9.0, 10.0, 13.0, 17.0, 21.0, 25.0])
I_plexi_uA  = np.array([4.53, 21.95, 39.06, 43.34, 56.28, 72.95, 89.70, 106.61])

# Convert to SI units (Hz and A)
f_plexi_Hz = f_plexi_kHz * 1e3
I_plexi_A  = I_plexi_uA * 1e-6

# Fit in native plot units: I [uA] vs f [kHz]
m1_plot, se_m1_plot, c1_plot, se_c1_plot, r2_1 = ols_fit(f_plexi_kHz, I_plexi_uA)

# Fit in SI units: I [A] vs f [Hz]
m1, se_m1, c1, se_c1, _ = ols_fit(f_plexi_Hz, I_plexi_A)

# Capacitance: C = slope / (2*pi*V_rms)
C_plexi = m1 / (2 * np.pi * V_rms)
# Uncertainty propagation: partial derivatives w.r.t slope and V_rms
se_C_plexi = C_plexi * np.sqrt((se_m1 / m1)**2 + (delta_V_rms / V_rms)**2)
se_C_plexi_stat = m1 / (2 * np.pi * V_rms) * (se_m1 / m1)

# Dielectric constant: K = C * d / (epsilon_0 * A)
K_plexi = (C_plexi * d_fixed) / (epsilon_0 * A_plate)
se_K_plexi = K_plexi * np.sqrt(
    (se_C_plexi / C_plexi)**2 +
    (delta_d / d_fixed)**2 +
    (delta_A / A_plate)**2
)
se_K_plexi_stat = K_plexi * (se_C_plexi_stat / C_plexi)

print("--- 1. PLEXIGLASS DIELECTRIC MEASUREMENT ---")
print(f"Slope (uA/kHz): {m1_plot:.4f} +/- {se_m1_plot:.4f}")
print(f"Intercept (uA): {c1_plot:.4f} +/- {se_c1_plot:.4f}")
print(f"R^2:            {r2_1:.6f}")
print(f"C_plexi:        {C_plexi*1e12:.2f} +/- {se_C_plexi_stat*1e12:.2f} (stat) +/- {se_C_plexi*1e12:.2f} (total) pF")
print(f"K_plexi:        {K_plexi:.3f} +/- {se_K_plexi_stat:.3f} (stat) +/- {se_K_plexi:.3f} (total)")
print()

# ---------------------------------------------------------
# Part 2: Air Dielectric (Table 2 from 7.tex)
# ---------------------------------------------------------
f_air_kHz = np.array([1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0])
I_air_uA  = np.array([1.02, 5.25, 9.50, 13.73, 17.18, 22.05, 26.14])

f_air_Hz = f_air_kHz * 1e3
I_air_A  = I_air_uA * 1e-6

m2_plot, se_m2_plot, c2_plot, se_c2_plot, r2_2 = ols_fit(f_air_kHz, I_air_uA)
m2, se_m2, c2, se_c2, _ = ols_fit(f_air_Hz, I_air_A)

C_air = m2 / (2 * np.pi * V_rms)
se_C_air = C_air * np.sqrt((se_m2 / m2)**2 + (delta_V_rms / V_rms)**2)
se_C_air_stat = C_air * (se_m2 / m2)

# Permittivity of Air: epsilon_air = C_air * d / A
eps_air = (C_air * d_fixed) / A_plate
se_eps_air = eps_air * np.sqrt(
    (se_C_air / C_air)**2 +
    (delta_d / d_fixed)**2 +
    (delta_A / A_plate)**2
)
se_eps_air_stat = eps_air * (se_C_air_stat / C_air)

K_air = eps_air / epsilon_0
se_K_air = se_eps_air / epsilon_0
se_K_air_stat = se_eps_air_stat / epsilon_0

# Relative error compared to vacuum permittivity
rel_diff_eps0 = np.abs(eps_air - epsilon_0) / epsilon_0 * 100.0

# Geometry-independent capacitance ratio
cap_ratio = m1 / m2
se_cap_ratio = cap_ratio * np.sqrt((se_m1 / m1)**2 + (se_m2 / m2)**2)

print("--- 2. AIR DIELECTRIC MEASUREMENT ---")
print(f"Slope (uA/kHz): {m2_plot:.4f} +/- {se_m2_plot:.4f}")
print(f"Intercept (uA): {c2_plot:.4f} +/- {se_c2_plot:.4f}")
print(f"R^2:            {r2_2:.6f}")
print(f"C_air:          {C_air*1e12:.2f} +/- {se_C_air_stat*1e12:.2f} (stat) +/- {se_C_air*1e12:.2f} (total) pF")
print(f"epsilon_air:    ({eps_air*1e12:.3f} +/- {se_eps_air_stat*1e12:.3f}) x 10^-12 F/m")
print(f"K_air:          {K_air:.3f} +/- {se_K_air_stat:.3f} (stat) +/- {se_K_air:.3f} (total)")
print(f"Discrepancy vs epsilon_0: {rel_diff_eps0:.1f}%")
print(f"Capacitance Ratio C_plexi / C_air: {cap_ratio:.3f} +/- {se_cap_ratio:.3f}")
print()

# ---------------------------------------------------------
# Part 3: Distance Dependence at f = 14 kHz (Table 3 from 7.tex)
# ---------------------------------------------------------
d_mm = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
I_d_uA = np.array([14.22, 12.34, 10.96, 9.90, 9.05, 8.37, 7.78])
inv_d_mm = 1.0 / d_mm  # in mm^-1

f_const = 14.0e3  # Hz
m3, se_m3, c3, se_c3, r2_3 = ols_fit(inv_d_mm, I_d_uA)

# Theoretical slope using epsilon_0:
# I = (2*pi*f*V*eps_0*A) * (1/d)
# Convert d in mm: 1/d[m] = 1000 * 1/d[mm]
slope_theory_uA_mm = (2 * np.pi * f_const * V_rms * epsilon_0 * A_plate) * 1e3 * 1e6
# With measured effective air permittivity:
slope_eff_uA_mm = (2 * np.pi * f_const * V_rms * eps_air * A_plate) * 1e3 * 1e6

# Stray capacitance from intercept c3:
# c3 = 2*pi*f*C_stray*V -> C_stray = c3 / (2*pi*f*V)
C_stray = (c3 * 1e-6) / (2 * np.pi * f_const * V_rms)
se_C_stray = (se_c3 * 1e-6) / (2 * np.pi * f_const * V_rms)

print("--- 3. DISTANCE DEPENDENCE AT f = 14 kHz ---")
print(f"Slope a3:       {m3:.3f} +/- {se_m3:.3f} uA*mm")
print(f"Intercept b3:   {c3:.3f} +/- {se_c3:.3f} uA")
print(f"R^2:            {r2_3:.6f}")
print(f"Theoretical slope (vacuum): {slope_theory_uA_mm:.2f} uA*mm")
print(f"Effective slope (eps_air):  {slope_eff_uA_mm:.2f} uA*mm")
print(f"Stray capacitance (C_stray): {C_stray*1e12:.2f} +/- {se_C_stray*1e12:.2f} pF")
print()

# ---------------------------------------------------------
# Part 4: Oscilloscope Calibration & Signal Analysis
# ---------------------------------------------------------
# Oscilloscope setting from Section 3 of 7.tex:
# Generator at 10 kHz, Time/Div = 50 us
time_per_div = 50e-6  # s
N_div = 2.0           # divisions per period
delta_N = 0.1         # visual reading uncertainty
T_scope = N_div * time_per_div
f_scope = 1.0 / T_scope
rel_err_f = delta_N / N_div
delta_f_scope = f_scope * rel_err_f

print("--- 4. OSCILLOSCOPE FREQUENCY MEASUREMENT ---")
print(f"Time/Div:        {time_per_div*1e6:.1f} us/div")
print(f"Observed period: N = {N_div:.1f} +/- {delta_N:.1f} div -> T = {T_scope*1e6:.1f} us")
print(f"Measured frequency: f = {f_scope*1e-3:.2f} +/- {delta_f_scope*1e-3:.2f} kHz")
print(f"Relative uncertainty: {rel_err_f*100:.1f}%")
print("="*65)

# ---------------------------------------------------------
# Plot 1: Current vs. Frequency (Plexiglass)
# ---------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(7, 4.8))
f_grid_plexi = np.linspace(0, 27, 200)
fit_line_plexi = m1_plot * f_grid_plexi + c1_plot

# Estimate instrumental uncertainty for error bars (0.5 uA or 1%)
err_I_plexi = np.maximum(0.3, 0.01 * I_plexi_uA)

ax1.errorbar(f_plexi_kHz, I_plexi_uA, yerr=err_I_plexi, fmt='s', color='#1f77b4',
             ecolor='#1f77b4', elinewidth=1.2, capsize=3, label='Experimental Data', zorder=3)
ax1.plot(f_grid_plexi, fit_line_plexi, color='#d62728', linestyle='-',
         label=f'OLS Fit: $I = ({m1_plot:.3f}\\pm{se_m1_plot:.3f})f + ({c1_plot:.2f}\\pm{se_c1_plot:.2f})$\n$R^2 = {r2_1:.5f}$')

ax1.set_xlabel('Frequency $f$ [kHz]')
ax1.set_ylabel('Current $I$ [$\\mu$A]')
ax1.set_title('Capacitor Current vs. Frequency (Plexiglass Dielectric)')
ax1.set_xlim(0, 27)
ax1.set_ylim(0, 115)
ax1.grid(True)
ax1.legend(loc='upper left', frameon=True)
plt.tight_layout()
p1_path = plots_dir / 'I_vs_f_plexi.pdf'
fig1.savefig(p1_path, dpi=300)
plt.close(fig1)
print(f"Saved: {p1_path}")

# ---------------------------------------------------------
# Plot 2: Current vs. Frequency (Air)
# ---------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7, 4.8))
f_grid_air = np.linspace(0, 27, 200)
fit_line_air = m2_plot * f_grid_air + c2_plot
err_I_air = np.maximum(0.1, 0.01 * I_air_uA)

ax2.errorbar(f_air_kHz, I_air_uA, yerr=err_I_air, fmt='o', color='#2ca02c',
             ecolor='#2ca02c', elinewidth=1.2, capsize=3, label='Experimental Data', zorder=3)
ax2.plot(f_grid_air, fit_line_air, color='#ff7f0e', linestyle='-',
         label=f'OLS Fit: $I = ({m2_plot:.3f}\\pm{se_m2_plot:.3f})f + ({c2_plot:.3f}\\pm{se_c2_plot:.3f})$\n$R^2 = {r2_2:.5f}$')

ax2.set_xlabel('Frequency $f$ [kHz]')
ax2.set_ylabel('Current $I$ [$\\mu$A]')
ax2.set_title('Capacitor Current vs. Frequency (Air Dielectric)')
ax2.set_xlim(0, 27)
ax2.set_ylim(0, 29)
ax2.grid(True)
ax2.legend(loc='upper left', frameon=True)
plt.tight_layout()
p2_path = plots_dir / 'I_vs_f_air.pdf'
fig2.savefig(p2_path, dpi=300)
plt.close(fig2)
print(f"Saved: {p2_path}")

# ---------------------------------------------------------
# Plot 3: Current vs. Inverse Distance (1/d)
# ---------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(7, 4.8))
inv_d_grid = np.linspace(0.09, 0.36, 200)
fit_line_d = m3 * inv_d_grid + c3
err_I_d = np.maximum(0.1, 0.01 * I_d_uA)

ax3.errorbar(inv_d_mm, I_d_uA, yerr=err_I_d, fmt='^', color='#9467bd',
             ecolor='#9467bd', elinewidth=1.2, capsize=3, label='Experimental Data', zorder=3)
ax3.plot(inv_d_grid, fit_line_d, color='#8c564b', linestyle='-',
         label=f'OLS Fit: $I = ({m3:.2f}\\pm{se_m3:.2f})(1/d) + ({c3:.2f}\\pm{se_c3:.2f})$\n$R^2 = {r2_3:.5f}$')

ax3.set_xlabel('Inverse Plate Distance $1/d$ [$\\mathrm{mm^{-1}}$]')
ax3.set_ylabel('Current $I$ [$\\mu$A]')
ax3.set_title('Current vs. Inverse Distance ($f = 14$ kHz, $V_m = 4$ V)')
ax3.set_xlim(0.09, 0.36)
ax3.set_ylim(6, 16)
ax3.grid(True)
ax3.legend(loc='lower right', frameon=True)
plt.tight_layout()
p3_path = plots_dir / 'I_vs_inv_d.pdf'
fig3.savefig(p3_path, dpi=300)
plt.close(fig3)
print(f"Saved: {p3_path}")

# ---------------------------------------------------------
# Plot 4: Oscilloscope Screen Display (Time-Domain Waveform)
# ---------------------------------------------------------
fig4, ax4 = plt.subplots(figsize=(7, 5))
# Simulate 8x10 graticule (standard CRT screen: 10 horizontal, 8 vertical div)
t_us = np.linspace(-250, 250, 1000)
V_signal = 10.0 * np.sin(2 * np.pi * 10e3 * (t_us * 1e-6))  # 10 kHz, 10 V amplitude

# Plot screen grid: horizontal divisions from -5 to +5 (50 us/div), vertical -4 to +4 (5 V/div)
ax4.set_facecolor('#0f1b0f')  # dark phosphor green background
ax4.plot(t_us / 50.0, V_signal / 5.0, color='#00ff66', lw=2.2, label='Ch 1: Sine Wave (10 kHz, $V_m=10$ V)')

# Grid lines
for h in range(-5, 6):
    ax4.axvline(h, color='#2a5a2a', lw=0.7, ls=':')
for v in range(-4, 5):
    ax4.axhline(v, color='#2a5a2a', lw=0.7, ls=':')
ax4.axvline(0, color='#3a7a3a', lw=1.2)
ax4.axhline(0, color='#3a7a3a', lw=1.2)

# Annotations showing Period measurement
ax4.annotate('', xy=(0, 0), xytext=(2.0, 0),
             arrowprops=dict(arrowstyle='<->', color='#ffff33', lw=1.8))
ax4.text(1.0, 0.4, '$T = 2.0\\,\\mathrm{div} \\times 50\\,\\mu\\mathrm{s/div} = 100\\,\\mu\\mathrm{s}$\n$f = 1/T = 10\\,\\mathrm{kHz}$',
         color='#ffff33', fontsize=10, ha='center',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#112211', edgecolor='#ffff33', alpha=0.9))

# Annotations showing Peak-to-Peak measurement
ax4.annotate('', xy=(-2.5, -2.0), xytext=(-2.5, 2.0),
             arrowprops=dict(arrowstyle='<->', color='#00ffff', lw=1.8))
ax4.text(-2.55, 0.0, '$V_{pp} = 4.0\\,\\mathrm{div} \\times 5\\,\\mathrm{V/div} = 20\\,\\mathrm{V}$',
         color='#00ffff', fontsize=10, ha='right', va='center', rotation=90,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#112211', edgecolor='#00ffff', alpha=0.9))

ax4.set_xlim(-5, 5)
ax4.set_ylim(-4, 4)
ax4.set_xlabel('Horizontal Scale: Time/Div = 50 $\\mu$s [divisions]', color='#00ff66')
ax4.set_ylabel('Vertical Scale: Volts/Div = 5 V [divisions]', color='#00ff66')
ax4.tick_params(colors='#00ff66')
ax4.set_title('Oscilloscope Time-Domain Measurement ($f=10$ kHz, $V_m=10$ V)', color='#00ff66', pad=12)
ax4.legend(loc='upper right', facecolor='#112211', edgecolor='#2a5a2a', labelcolor='#00ff66')
plt.tight_layout()
p4_path = plots_dir / 'oscilloscope_waveform.pdf'
fig4.savefig(p4_path, dpi=300)
plt.close(fig4)
print(f"Saved: {p4_path}")

# ---------------------------------------------------------
# Plot 5: Lissajous Figures for Phase Measurement
# ---------------------------------------------------------
fig5, axes = plt.subplots(1, 4, figsize=(13, 3.4))
phases = [0, np.pi/6, np.pi/4, np.pi/2]
phase_labels = ['0^\\circ', '30^\\circ', '45^\\circ', '90^\\circ']
t_lis = np.linspace(0, 2*np.pi, 500)

for ax, phi, plabel in zip(axes, phases, phase_labels):
    x_sig = np.sin(t_lis)
    y_sig = np.sin(t_lis + phi)
    ax.plot(x_sig, y_sig, color='#1f77b4', lw=2.0)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.axhline(0, color='black', lw=0.6)
    ax.axvline(0, color='black', lw=0.6)
    
    # Intercept at x = 0 is y0 = sin(phi)
    y0 = np.sin(phi)
    if phi > 0 and phi < np.pi/2:
        ax.plot([0, 0], [0, y0], 'ro-', markersize=4)
        ax.text(0.08, y0/2, '$y_0$', color='red', fontsize=11, va='center')
        ax.plot([0, 0], [0, 1.0], 'k:', lw=1)
        ax.text(-0.08, 0.9, '$y_{\\max}$', color='black', fontsize=9, ha='right')
    ax.set_title(f'$\\Delta\\phi = {plabel}$\\n$\\sin\\Delta\\phi = y_0/y_{{\\max}}$', fontsize=11)

plt.suptitle('Lissajous Patterns (X-Y Mode) for Phase Shift Determination', fontsize=13, y=1.05)
plt.tight_layout()
p5_path = plots_dir / 'lissajous_figures.pdf'
fig5.savefig(p5_path, dpi=300, bbox_inches='tight')
plt.close(fig5)
print(f"Saved: {p5_path}")

print("="*65)
print("ALL ANALYSES & PLOTS SUCCESSFULLY GENERATED!")
print("="*65)
