#!/usr/bin/env python3
"""
Analysis of Experiment 8: AC Circuits, Phase Difference via Lissajous Figures,
and Determination of Inductance (L), Capacitance (C), and Resonance Frequency (f_res).

Physics Laboratory II - Sharif University of Technology
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
    'errorbar.capsize': 3.5,
    'grid.alpha': 0.4,
    'grid.linestyle': '--'
})

# Ensure plots directory exists
plots_dir = Path('plots')
plots_dir.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# Regression Helpers
# -----------------------------------------------------------------------------
def ols_fit(x, y):
    """
    Ordinary Least Squares (OLS) linear regression: y = m * x + c
    Returns: m, se_m, c, se_c, r2
    """
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

def ols_origin_fit(x, y):
    """
    Linear regression through origin: y = a * x
    Returns: a, se_a, r2, s_yx
    """
    N = len(x)
    a = np.sum(x * y) / np.sum(x**2)
    y_pred = a * x
    residuals = y - y_pred
    sse = np.sum(residuals**2)
    s2 = sse / (N - 1) if N > 1 else 0.0
    se_a = np.sqrt(s2 / np.sum(x**2))
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (sse / ss_tot) if ss_tot > 0 else 1.0
    return a, se_a, r2, np.sqrt(s2)

def propagate_tan_error(y0, B, delta_y0=0.04, delta_B=0.04):
    """
    Propagate uncertainties on y0 and B to sin(phi) and tan(phi):
    sin(phi) = y0 / B
    tan(phi) = sin(phi) / sqrt(1 - sin^2(phi))
    sigma_tan = sigma_sin / (1 - sin^2(phi))^(3/2)
    """
    sin_phi = y0 / B
    sigma_sin = np.sqrt((delta_y0 / B)**2 + (y0 * delta_B / (B**2))**2)
    sigma_tan = sigma_sin / ((1.0 - sin_phi**2)**1.5)
    return sin_phi, sigma_sin, sigma_tan

# =============================================================================
# 1. Unknown Frequency via Lissajous Figures
# =============================================================================
print("=" * 70)
print("1. DETERMINATION OF UNKNOWN FREQUENCY VIA LISSAJOUS FIGURES")
print("=" * 70)

f_x = 100.0  # Reference frequency in Hz
delta_fx = 0.5  # Uncertainty in Hz
Nx = 2
Ny = 4
fy = f_x * (Nx / Ny)
delta_fy = fy * (delta_fx / f_x)  # Assuming contact counts are exact integers

print(f"Reference frequency (f_x): {f_x:.1f} +/- {delta_fx:.1f} Hz")
print(f"Tangency ratio (Nx / Ny): {Nx} / {Ny} = {Nx/Ny:.2f}")
print(f"Unknown frequency (f_y):  {fy:.1f} +/- {delta_fy:.2f} Hz\n")

# =============================================================================
# 2. RL Series Circuit: Determination of Self-Inductance (L)
# =============================================================================
print("=" * 70)
print("2. RL SERIES CIRCUIT: DETERMINATION OF INDUCTANCE (L)")
print("=" * 70)

R_nominal = 300.0  # Ohms
delta_R = 3.0      # 1% resistor box tolerance (Ohms)

f_RL = np.array([30.0, 60.0, 90.0, 120.0])  # Hz
y0_RL = np.array([0.88, 1.48, 1.76, 1.88])   # div
B_RL = np.array([1.92, 2.08, 2.13, 2.26])    # div

sin_phi_RL, d_sin_phi_RL, d_tan_phi_RL = propagate_tan_error(y0_RL, B_RL)
phi_RL_rad = np.arcsin(sin_phi_RL)
tan_phi_RL = np.tan(phi_RL_rad)

print(f"{'f (Hz)':>8} {'y0':>6} {'B':>6} {'sin(phi)':>10} {'phi (rad)':>10} {'tan(phi)':>10} {'d(tan(phi))':>12}")
for i in range(len(f_RL)):
    print(f"{f_RL[i]:8.1f} {y0_RL[i]:6.2f} {B_RL[i]:6.2f} {sin_phi_RL[i]:10.4f} "
          f"{phi_RL_rad[i]:10.4f} {tan_phi_RL[i]:10.4f} {d_tan_phi_RL[i]:12.4f}")

# Origin-constrained fit: tan(phi) = a * f, where a = 2*pi*L / R
a_RL, se_a_RL, r2_RL_origin, s_RL = ols_origin_fit(f_RL, tan_phi_RL)
L_RL = (a_RL * R_nominal) / (2.0 * np.pi)
se_L_RL = L_RL * np.sqrt((se_a_RL / a_RL)**2 + (delta_R / R_nominal)**2)

# Unconstrained OLS fit for comparison
m_RL, se_m_RL, c_RL, se_c_RL, r2_RL_ols = ols_fit(f_RL, tan_phi_RL)
L_RL_ols = (m_RL * R_nominal) / (2.0 * np.pi)
se_L_RL_ols = L_RL_ols * np.sqrt((se_m_RL / m_RL)**2 + (delta_R / R_nominal)**2)

print("\n--- RL Regression Results ---")
print(f"Origin fit:        slope a = {a_RL:.5e} +/- {se_a_RL:.5e} Hz^-1, R^2 = {r2_RL_origin:.4f}")
print(f"Derived L:         L = {L_RL:.4f} +/- {se_L_RL:.4f} H  [(0.69 +/- 0.06) H in report]")
print(f"Unconstrained fit: slope m = {m_RL:.5e} +/- {se_m_RL:.5e} Hz^-1, c = {c_RL:.4f} +/- {se_c_RL:.4f}, R^2 = {r2_RL_ols:.4f}")
print(f"Unconstrained L:   L = {L_RL_ols:.4f} +/- {se_L_RL_ols:.4f} H\n")

# Residuals for origin fit
res_RL = tan_phi_RL - a_RL * f_RL
sse_RL = np.sum(res_RL**2)
print(f"Sum of Squared Errors (SSE): {sse_RL:.5e}")
for i in range(len(f_RL)):
    print(f"  Point {i+1} (f={f_RL[i]:3.0f} Hz): y_meas={tan_phi_RL[i]:.4f}, y_fit={a_RL*f_RL[i]:.4f}, residual={res_RL[i]:+.4f}")
print()

# =============================================================================
# 3. RC Series Circuit: Determination of Capacitance (C)
# =============================================================================
print("=" * 70)
print("3. RC SERIES CIRCUIT: DETERMINATION OF CAPACITANCE (C)")
print("=" * 70)

C_nominal = 10.0e-6  # 10 uF

f_RC = np.array([30.0, 60.0, 90.0, 120.0])  # Hz
inv_f_RC = 1.0 / f_RC                       # s
y0_RC = np.array([1.72, 1.28, 0.96, 0.76])   # div
B_RC = np.array([2.06, 2.04, 2.09, 1.96])    # div

sin_phi_RC, d_sin_phi_RC, d_tan_phi_RC = propagate_tan_error(y0_RC, B_RC)
abs_tan_phi_RC = sin_phi_RC / np.sqrt(1.0 - sin_phi_RC**2)
tan_phi_RC = -abs_tan_phi_RC  # Physical sign: capacitive circuit has phi < 0

print(f"{'f (Hz)':>8} {'1/f (s)':>10} {'y0':>6} {'B':>6} {'sin(phi)':>10} {'tan(phi)':>10} {'d(tan(phi))':>12}")
for i in range(len(f_RC)):
    print(f"{f_RC[i]:8.1f} {inv_f_RC[i]:10.5f} {y0_RC[i]:6.2f} {B_RC[i]:6.2f} "
          f"{sin_phi_RC[i]:10.4f} {tan_phi_RC[i]:10.4f} {d_tan_phi_RC[i]:12.4f}")

# Origin-constrained fit: tan(phi) = m * (1/f), where m = -1 / (2*pi*R*C)
m_origin_RC, se_m_origin_RC, r2_RC_origin, s_RC = ols_origin_fit(inv_f_RC, tan_phi_RC)
C_origin_RC = -1.0 / (2.0 * np.pi * R_nominal * m_origin_RC)
se_C_origin_RC = C_origin_RC * np.sqrt((se_m_origin_RC / abs(m_origin_RC))**2 + (delta_R / R_nominal)**2)

# Unconstrained OLS fit
m_RC, se_m_RC, c_RC, se_c_RC, r2_RC_ols = ols_fit(inv_f_RC, tan_phi_RC)
C_ols_RC = -1.0 / (2.0 * np.pi * R_nominal * m_RC)
se_C_ols_RC = C_ols_RC * np.sqrt((se_m_RC / abs(m_RC))**2 + (delta_R / R_nominal)**2)

print("\n--- RC Regression Results ---")
print(f"Origin fit:        slope m = {m_origin_RC:.4f} +/- {se_m_origin_RC:.4f} s^-1, R^2 = {r2_RC_origin:.4f}")
print(f"Derived C:         C = {C_origin_RC*1e6:.2f} +/- {se_C_origin_RC*1e6:.2f} uF  [(11.5 +/- 0.2) uF in report]")
print(f"Nominal C diff:    Relative difference = {abs(C_origin_RC - C_nominal)/C_nominal * 100:.1f}%")
print(f"Unconstrained fit: slope m = {m_RC:.4f} +/- {se_m_RC:.4f} s^-1, c = {c_RC:.4f} +/- {se_c_RC:.4f}, R^2 = {r2_RC_ols:.4f}")
print(f"Unconstrained C:   C = {C_ols_RC*1e6:.2f} +/- {se_C_ols_RC*1e6:.2f} uF\n")

# Residuals for origin fit
res_RC = tan_phi_RC - m_origin_RC * inv_f_RC
sse_RC = np.sum(res_RC**2)
print(f"Sum of Squared Errors (SSE): {sse_RC:.5e}")
for i in range(len(f_RC)):
    print(f"  Point {i+1} (1/f={inv_f_RC[i]:.5f} s): y_meas={tan_phi_RC[i]:.4f}, y_fit={m_origin_RC*inv_f_RC[i]:.4f}, residual={res_RC[i]:+.4f}")
print()

# =============================================================================
# 4. RLC Series Circuit: Resonance and Comparison
# =============================================================================
print("=" * 70)
print("4. RLC SERIES CIRCUIT: RESONANCE FREQUENCY & L DETERMINATION")
print("=" * 70)

f_RLC = np.array([20.0, 35.0, 48.8, 65.0, 80.0])  # Hz
y0_RLC = np.array([1.72, 1.04, 0.0, 0.88, 1.32])  # div
B_RLC = np.array([2.00, 1.94, 2.00, 2.00, 2.08])   # div (B at res set to 2.0 for reference)

# Signs: capacitive (negative) below resonance, inductive (positive) above
signs_RLC = np.array([-1.0, -1.0, 0.0, 1.0, 1.0])
sin_phi_RLC = y0_RLC / B_RLC
abs_tan_phi_RLC = np.where(sin_phi_RLC < 1.0, sin_phi_RLC / np.sqrt(np.maximum(1e-9, 1.0 - sin_phi_RLC**2)), 0.0)
tan_phi_RLC = signs_RLC * abs_tan_phi_RLC

d_tan_phi_RLC = np.zeros_like(tan_phi_RLC)
for i in range(len(f_RLC)):
    if y0_RLC[i] == 0.0:
        d_tan_phi_RLC[i] = 0.02
    else:
        _, _, dt = propagate_tan_error(y0_RLC[i], B_RLC[i])
        d_tan_phi_RLC[i] = dt

print(f"{'f (Hz)':>8} {'y0':>6} {'B':>6} {'|sin(phi)|':>12} {'tan(phi)':>12} {'d(tan(phi))':>12}")
for i in range(len(f_RLC)):
    print(f"{f_RLC[i]:8.1f} {y0_RLC[i]:6.2f} {B_RLC[i]:6.2f} {sin_phi_RLC[i]:12.4f} "
          f"{tan_phi_RLC[i]:12.4f} {d_tan_phi_RLC[i]:12.4f}")

# Resonance frequency and derived inductance
f_res = 48.8
delta_f_res = 0.5  # Hz

# Using nominal C (10 uF)
L_res_nom = 1.0 / ((2.0 * np.pi * f_res)**2 * C_nominal)
se_L_res_nom = L_res_nom * 2.0 * (delta_f_res / f_res)

# Using measured C (11.45 uF)
L_res_meas = 1.0 / ((2.0 * np.pi * f_res)**2 * C_origin_RC)
se_L_res_meas = L_res_meas * np.sqrt((2.0 * delta_f_res / f_res)**2 + (se_C_origin_RC / C_origin_RC)**2)

# Theoretical resonance from RL and RC measurements
f0_theory = 1.0 / (2.0 * np.pi * np.sqrt(L_RL * C_origin_RC))
rel_err_f0 = 0.5 * np.sqrt((se_L_RL / L_RL)**2 + (se_C_origin_RC / C_origin_RC)**2)
se_f0_theory = f0_theory * rel_err_f0

print("\n--- Resonance Results & Cross-Validation ---")
print(f"Observed resonance frequency f_res:      {f_res:.1f} +/- {delta_f_res:.1f} Hz")
print(f"L derived from f_res (with C_nom=10uF):  {L_res_nom:.3f} +/- {se_L_res_nom:.3f} H  [(1.06 +/- 0.02) H in report]")
print(f"L derived from f_res (with C_meas):      {L_res_meas:.3f} +/- {se_L_res_meas:.3f} H")
print(f"L derived from RL phase slope:           {L_RL:.3f} +/- {se_L_RL:.3f} H")
print(f"Predicted resonance f0 from (L_RL, C_RC): {f0_theory:.1f} +/- {se_f0_theory:.1f} Hz  [(56.8 +/- 2.4) Hz in report]")
print(f"Discrepancy (f0_theory - f_res):         {f0_theory - f_res:.1f} Hz (attributed to coil internal resistance R_L & capacitor ESR)")
print("=" * 70 + "\n")

# =============================================================================
# 5. Figure 1: RL Circuit Plot
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 4.8))
f_grid = np.linspace(0, 130, 200)

ax.errorbar(f_RL, tan_phi_RL, yerr=d_tan_phi_RL, fmt='s', color='#1f77b4',
            ecolor='#1f77b4', elinewidth=1.5, capsize=4, label='Experimental Data', zorder=4)

ax.plot(f_grid, a_RL * f_grid, color='#d62728', linestyle='-',
        label=f'Origin Fit: $\\tan\\phi = ({a_RL*1e2:.3f}\\times 10^{{-2}})\\,f$\n($L = {L_RL:.2f} \\pm {se_L_RL:.2f}$ H, $R^2 = {r2_RL_origin:.4f}$)')

ax.plot(f_grid, m_RL * f_grid + c_RL, color='#2ca02c', linestyle='--', alpha=0.7,
        label=f'OLS Fit: $\\tan\\phi = ({m_RL*1e2:.3f}\\times 10^{{-2}})\\,f + {c_RL:.3f}$\n($R^2 = {r2_RL_ols:.4f}$)')

ax.set_xlabel('Frequency $f$ (Hz)')
ax.set_ylabel(r'$\tan\phi$')
ax.set_title(r'RL Series Circuit: $\tan\phi$ vs. Frequency ($R = 300\ \Omega$)')
ax.set_xlim(0, 130)
ax.set_ylim(0, 1.8)
ax.grid(True)
ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=False)
plt.tight_layout()
fig.savefig(plots_dir / 'RL_tan_vs_f.pdf')
plt.close(fig)
print("Saved plots/RL_tan_vs_f.pdf")

# =============================================================================
# 6. Figure 2: RC Circuit Plot
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 4.8))
inv_f_grid = np.linspace(0, 0.04, 200)

ax.errorbar(inv_f_RC, tan_phi_RC, yerr=d_tan_phi_RC, fmt='s', color='#ff7f0e',
            ecolor='#ff7f0e', elinewidth=1.5, capsize=4, label='Experimental Data', zorder=4)

ax.plot(inv_f_grid, m_origin_RC * inv_f_grid, color='#1f77b4', linestyle='-',
        label=f'Origin Fit: $\\tan\\phi = ({m_origin_RC:.2f})\\,(1/f)$\n($C = {C_origin_RC*1e6:.1f} \\pm {se_C_origin_RC*1e6:.1f}\\ \\mu$F, $R^2 = {r2_RC_origin:.4f}$)')

ax.plot(inv_f_grid, m_RC * inv_f_grid + c_RC, color='#2ca02c', linestyle='--', alpha=0.7,
        label=f'OLS Fit: $\\tan\\phi = ({m_RC:.2f})\\,(1/f) + {c_RC:.3f}$\n($R^2 = {r2_RC_ols:.4f}$)')

ax.axhline(0, color='gray', linestyle=':', linewidth=0.9)
ax.set_xlabel(r'Inverse Frequency $1/f$ (s)')
ax.set_ylabel(r'$\tan\phi$')
ax.set_title(r'RC Series Circuit: $\tan\phi$ vs. $1/f$ ($R = 300\ \Omega$)')
ax.set_xlim(0, 0.038)
ax.set_ylim(-1.8, 0.1)
ax.grid(True)
ax.legend(loc='lower left', frameon=True, fancybox=True, shadow=False)
plt.tight_layout()
fig.savefig(plots_dir / 'RC_tan_vs_inv_f.pdf')
plt.close(fig)
print("Saved plots/RC_tan_vs_inv_f.pdf")

# =============================================================================
# 7. Figure 3: RLC Circuit Plot (Resonance Curve)
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 4.8))
ax.errorbar(f_RLC, tan_phi_RLC, yerr=d_tan_phi_RLC, fmt='s', color='#d62728',
            ecolor='#d62728', elinewidth=1.5, capsize=4, label='Experimental Data', zorder=4)

# Theoretical RLC curve using extracted L and C
f_curve = np.linspace(15, 85, 300)
omega_curve = 2.0 * np.pi * f_curve
tan_phi_theory = (omega_curve * L_RL - 1.0 / (omega_curve * C_origin_RC)) / R_nominal
ax.plot(f_curve, tan_phi_theory, color='#1f77b4', linestyle='--', alpha=0.6,
        label=r'Theory from $(L_{RL}, C_{RC})$: $\tan\phi = \frac{\omega L - 1/(\omega C)}{R}$')

# Empirical smooth spline through experimental points for visualization
ax.axhline(0, color='black', linestyle='-', linewidth=1.0)
ax.axvline(f_res, color='purple', linestyle=':', linewidth=1.4,
           label=f'Resonance: $f_{{res}} = {f_res:.1f}$ Hz ($\\phi = 0$)')

ax.annotate(r'Capacitive ($\phi < 0$)', xy=(25, -1.2), xytext=(22, -1.5),
            color='#333333', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='goldenrod'))
ax.annotate(r'Inductive ($\phi > 0$)', xy=(70, 0.6), xytext=(65, 0.9),
            color='#333333', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='goldenrod'))

ax.set_xlabel('Frequency $f$ (Hz)')
ax.set_ylabel(r'$\tan\phi$')
ax.set_title(r'RLC Series Circuit: Phase Transition across Resonance ($f_{\mathrm{res}} = 48.8$ Hz)')
ax.set_xlim(15, 85)
ax.set_ylim(-2.0, 1.2)
ax.grid(True)
ax.legend(loc='lower right', frameon=True, fancybox=True, shadow=False)
plt.tight_layout()
fig.savefig(plots_dir / 'RLC_tan_vs_f.pdf')
plt.close(fig)
print("Saved plots/RLC_tan_vs_f.pdf")

# =============================================================================
# 8. Figure 4: Lissajous Figures Simulation
# =============================================================================
t = np.linspace(0, 2.0 * np.pi, 500)
fig, axs = plt.subplots(2, 3, figsize=(10.5, 7.0))
phases = [0, np.pi/6, np.pi/2, 5*np.pi/6, np.pi]
phase_titles = [
    r'$\phi = 0^\circ$ (In-phase / Resonance)',
    r'$\phi = 30^\circ$ (Ellipse)',
    r'$\phi = 90^\circ$ (Perpendicular axes)',
    r'$\phi = 150^\circ$ (Oblique Ellipse)',
    r'$\phi = 180^\circ$ (Anti-phase Line)'
]

for idx in range(5):
    row, col = idx // 3, idx % 3
    ax = axs[row, col]
    phi_val = phases[idx]
    x_val = np.sin(t)
    y_val = np.sin(t + phi_val)
    ax.plot(x_val, y_val, color='#1f77b4', linewidth=2.0)
    ax.set_title(phase_titles[idx], fontsize=10)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axhline(0, color='gray', linewidth=0.6, linestyle=':')
    ax.axvline(0, color='gray', linewidth=0.6, linestyle=':')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

# 6th panel: Frequency ratio fx : fy = 2 : 1 (Nx : Ny = 2 : 4 = 1 : 2)
ax_freq = axs[1, 2]
t_freq = np.linspace(0, 2.0 * np.pi, 800)
x_freq = np.sin(2.0 * t_freq)
y_freq = np.sin(1.0 * t_freq + np.pi/4)
ax_freq.plot(x_freq, y_freq, color='#d62728', linewidth=2.0)
ax_freq.set_title(r'$f_x : f_y = 2 : 1$ ($f_x = 100$ Hz, $f_y = 50$ Hz)', fontsize=10)
ax_freq.set_xlim(-1.2, 1.2)
ax_freq.set_ylim(-1.2, 1.2)
ax_freq.axhline(0, color='gray', linewidth=0.6, linestyle=':')
ax_freq.axvline(0, color='gray', linewidth=0.6, linestyle=':')
ax_freq.set_aspect('equal')
ax_freq.grid(True, alpha=0.3)

fig.suptitle('Lissajous Figure Geometries in XY Mode ($x(t) = A\\sin(\\omega_x t), y(t) = B\\sin(\\omega_y t + \\phi)$)', fontsize=13)
plt.tight_layout()
fig.savefig(plots_dir / 'Lissajous_simulation.pdf')
plt.close(fig)
print("Saved plots/Lissajous_simulation.pdf")

print("\nAnalysis complete. All figures successfully generated.")
