"""
Experiment 4: RC Circuits — Charging, Discharging, Time Constant, and Capacitance
Lab-Phy-II, Sharif University of Technology

Analyses:
  Part 1 - Charging capacitor C1 (20 uF) through voltmeter Rv1
  Part 2 - Discharging capacitor C2 (4 uF) through voltmeter Rv2
  Part 3 - Discharging series combination (C1 and C2) through voltmeter Rv
  Part 4 - Discharging parallel combination (C1 and C2) through voltmeter Rv
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Create plots directory
plots_dir = Path('plots')
plots_dir.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# Instrument Uncertainties
# -----------------------------------------------------------------------------
delta_V = 0.02       # Voltmeter reading uncertainty (V)
delta_t = 0.5        # Timing uncertainty (s)
delta_C1_nom = 1.0e-6 # 5% nominal tolerance on C1 (F)
delta_C2_nom = 0.2e-6 # 5% nominal tolerance on C2 (F)

# -----------------------------------------------------------------------------
# OLS Regression Helper
# -----------------------------------------------------------------------------
def ols_fit(x, y):
    """
    Ordinary Least Squares (OLS) linear regression helper.
    Returns: slope, se_slope, intercept, se_intercept, r2
    """
    N = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals ** 2)
    s_yx = np.sqrt(ss_res / (N - 2))
    ss_x = np.sum((x - np.mean(x)) ** 2)
    se_m = s_yx / np.sqrt(ss_x)
    se_c = s_yx * np.sqrt(1.0 / N + (np.mean(x) ** 2) / ss_x)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    return m, se_m, c, se_c, r2

# -----------------------------------------------------------------------------
# 1. Raw Experimental Data
# -----------------------------------------------------------------------------
# Part 1: C1 = 20 uF charging via DC source (V0 = 10 V) and voltmeter Rv
t1 = np.array([0, 15, 30, 45, 60, 75, 90, 105, 120, 135], dtype=float)
V1 = np.array([10.00, 9.19, 8.46, 7.79, 7.17, 6.59, 6.07, 5.59, 5.14, 4.73])
V0_1 = 10.00
C1_nom = 20.0e-6  # 20 uF

# Part 2: C2 = 4 uF discharging through voltmeter Rv
t2 = np.array([0, 15, 30, 45, 60, 75, 90, 105, 120, 135], dtype=float)
V2 = np.array([10.00, 7.36, 5.37, 3.92, 2.87, 2.11, 1.54, 1.13, 0.83, 0.60])
V0_2 = 10.00
C2_nom = 4.0e-6   # 4 uF

# Part 3: C1 and C2 in series discharging through voltmeter Rv
t3 = np.array([0, 15, 30, 45, 60, 75, 90, 105, 120, 135], dtype=float)
V3 = np.array([10.00, 6.95, 4.63, 3.12, 2.11, 1.42, 0.96, 0.65, 0.48, 0.30])
V0_3 = 10.00
C_series_th = (C1_nom * C2_nom) / (C1_nom + C2_nom)  # 3.33 uF

# Part 4: C1 and C2 in parallel discharging through voltmeter Rv
t4 = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300], dtype=float)
V4 = np.array([9.83, 8.65, 7.59, 6.65, 5.83, 5.11, 4.48, 3.93, 3.46, 3.02, 2.65])
V0_4 = 9.83
C_parallel_th = C1_nom + C2_nom  # 24.0 uF

# -----------------------------------------------------------------------------
# 2. OLS Fits on ln(V / V0) vs t
# -----------------------------------------------------------------------------
# For exponential decay V(t) = V0 * exp(-t / tau):
# y = ln(V/V0) = -t / tau = m * t + c  ==>  tau = -1 / m
# Also log10(V/V0) = -t / (tau * ln 10) = b * t + a ==> tau = -1 / (b * ln 10)

def analyze_decay(t_arr, V_arr, V0_val):
    ratio = V_arr / V0_val
    y_ln = np.log(ratio)
    y_log10 = np.log10(ratio)
    
    # Error in ln(V/V0): delta_y = delta_V / V
    delta_y_ln = delta_V / V_arr
    
    # Fits
    m_ln, se_m_ln, c_ln, se_c_ln, r2_ln = ols_fit(t_arr, y_ln)
    m_log10, se_m_log10, c_log10, se_c_log10, r2_log10 = ols_fit(t_arr, y_log10)
    
    tau = -1.0 / m_ln
    # Error propagation: delta_tau = |d(tau)/dm| * se_m = (1 / m^2) * se_m = tau * (se_m / |m|)
    delta_tau = tau * (se_m_ln / abs(m_ln))
    
    return {
        'ratio': ratio,
        'y_ln': y_ln,
        'y_log10': y_log10,
        'delta_y_ln': delta_y_ln,
        'm_ln': m_ln,
        'se_m_ln': se_m_ln,
        'c_ln': c_ln,
        'se_c_ln': se_c_ln,
        'r2_ln': r2_ln,
        'm_log10': m_log10,
        'se_m_log10': se_m_log10,
        'c_log10': c_log10,
        'se_c_log10': se_c_log10,
        'r2_log10': r2_log10,
        'tau': tau,
        'delta_tau': delta_tau
    }

res1 = analyze_decay(t1, V1, V0_1)
res2 = analyze_decay(t2, V2, V0_2)
res3 = analyze_decay(t3, V3, V0_3)
res4 = analyze_decay(t4, V4, V0_4)

# -----------------------------------------------------------------------------
# 3. Voltmeter Internal Resistance (Rv)
# -----------------------------------------------------------------------------
# Rv1 = tau1 / C1
Rv1 = res1['tau'] / C1_nom
dRv1_fit = res1['delta_tau'] / C1_nom
# Total error including nominal capacitor tolerance:
dRv1_tot = Rv1 * np.sqrt((res1['delta_tau'] / res1['tau']) ** 2 + (delta_C1_nom / C1_nom) ** 2)

# Rv2 = tau2 / C2
Rv2 = res2['tau'] / C2_nom
dRv2_fit = res2['delta_tau'] / C2_nom
dRv2_tot = Rv2 * np.sqrt((res2['delta_tau'] / res2['tau']) ** 2 + (delta_C2_nom / C2_nom) ** 2)

# Mean Rv
Rv_mean = 0.5 * (Rv1 + Rv2)
dRv_mean_fit = 0.5 * np.sqrt(dRv1_fit ** 2 + dRv2_fit ** 2)
dRv_mean_tot = 0.5 * np.sqrt(dRv1_tot ** 2 + dRv2_tot ** 2)
# Difference between the two measurements
Rv_diff = abs(Rv2 - Rv1)

# -----------------------------------------------------------------------------
# 4. Equivalent Capacitances (Series & Parallel)
# -----------------------------------------------------------------------------
# C_eq = tau / Rv_mean
# Series
C_series_exp = res3['tau'] / Rv_mean
dC_series_fit = C_series_exp * np.sqrt((res3['delta_tau'] / res3['tau']) ** 2 + (dRv_mean_fit / Rv_mean) ** 2)
dC_series_tot = C_series_exp * np.sqrt((res3['delta_tau'] / res3['tau']) ** 2 + (dRv_mean_tot / Rv_mean) ** 2)
err_series_pct = abs(C_series_exp - C_series_th) / C_series_th * 100.0

# Parallel
C_parallel_exp = res4['tau'] / Rv_mean
dC_parallel_fit = C_parallel_exp * np.sqrt((res4['delta_tau'] / res4['tau']) ** 2 + (dRv_mean_fit / Rv_mean) ** 2)
dC_parallel_tot = C_parallel_exp * np.sqrt((res4['delta_tau'] / res4['tau']) ** 2 + (dRv_mean_tot / Rv_mean) ** 2)
err_parallel_pct = abs(C_parallel_exp - C_parallel_th) / C_parallel_th * 100.0

# Linear interpolation method for tau (as done in report text at V = V0/e)
def find_tau_interp(t_arr, V_arr, V0_val):
    target_V = V0_val / np.e
    # Find points bounding target_V
    idx = np.where(V_arr <= target_V)[0][0]
    t_a, t_b = t_arr[idx - 1], t_arr[idx]
    V_a, V_b = V_arr[idx - 1], V_arr[idx]
    slope = (V_b - V_a) / (t_b - t_a)
    tau_est = t_a + (target_V - V_a) / slope
    return target_V, tau_est

v_e3, tau_interp3 = find_tau_interp(t3, V3, V0_3)
v_e4, tau_interp4 = find_tau_interp(t4, V4, V0_4)

# -----------------------------------------------------------------------------
# 5. Print Results
# -----------------------------------------------------------------------------
print("=" * 75)
print("EXPERIMENT 4: RC CIRCUITS — REGRESSION & UNCERTAINTY ANALYSIS")
print("Sharif University of Technology - Physics Lab II")
print("=" * 75)

print("\n--- PART 1: Charging Capacitor C1 (20 uF) ---")
print(f"  Linear fit (ln V/V0 vs t):")
print(f"    Slope m       = ({res1['m_ln']:.6e} +/- {res1['se_m_ln']:.6e}) s^-1")
print(f"    Intercept c   = ({res1['c_ln']:.6e} +/- {res1['se_c_ln']:.6e})")
print(f"    R^2           = {res1['r2_ln']:.6f}")
print(f"  Base-10 fit (log10 V/V0 vs t):")
print(f"    Slope b       = ({res1['m_log10']:.6e} +/- {res1['se_m_log10']:.6e}) s^-1")
print(f"    Intercept a   = ({res1['c_log10']:.6e} +/- {res1['se_c_log10']:.6e})")
print(f"  Time constant tau_1  = {res1['tau']:.3f} +/- {res1['delta_tau']:.3f} s")
print(f"  Voltmeter Rv1        = {Rv1*1e-6:.4f} +/- {dRv1_fit*1e-6:.4f} (fit) +/- {dRv1_tot*1e-6:.4f} (tot) MOhm")

print("\n--- PART 2: Discharging Capacitor C2 (4 uF) ---")
print(f"  Linear fit (ln V/V0 vs t):")
print(f"    Slope m       = ({res2['m_ln']:.6e} +/- {res2['se_m_ln']:.6e}) s^-1")
print(f"    Intercept c   = ({res2['c_ln']:.6e} +/- {res2['se_c_ln']:.6e})")
print(f"    R^2           = {res2['r2_ln']:.6f}")
print(f"  Base-10 fit (log10 V/V0 vs t):")
print(f"    Slope b       = ({res2['m_log10']:.6e} +/- {res2['se_m_log10']:.6e}) s^-1")
print(f"    Intercept a   = ({res2['c_log10']:.6e} +/- {res2['se_c_log10']:.6e})")
print(f"  Time constant tau_2  = {res2['tau']:.3f} +/- {res2['delta_tau']:.3f} s")
print(f"  Voltmeter Rv2        = {Rv2*1e-6:.4f} +/- {dRv2_fit*1e-6:.4f} (fit) +/- {dRv2_tot*1e-6:.4f} (tot) MOhm")

print("\n--- VOLTMETER RESISTANCE SUMMARY ---")
print(f"  Rv1 (from C1)        = {Rv1*1e-6:.3f} MOhm")
print(f"  Rv2 (from C2)        = {Rv2*1e-6:.3f} MOhm")
print(f"  Mean Rv              = {Rv_mean*1e-6:.3f} +/- {dRv_mean_fit*1e-6:.3f} (fit) +/- {dRv_mean_tot*1e-6:.3f} (tot) MOhm")
print(f"  Absolute difference  = {Rv_diff*1e-6:.3f} MOhm (relative: {Rv_diff/Rv_mean*100:.1f}%)")

print("\n--- PART 3: Series Combination (C1 = 20 uF, C2 = 4 uF) ---")
print(f"  Theoretical C_eq     = {C_series_th*1e6:.3f} uF")
print(f"  Linear fit (ln V/V0 vs t):")
print(f"    Slope m       = ({res3['m_ln']:.6e} +/- {res3['se_m_ln']:.6e}) s^-1")
print(f"    Intercept c   = ({res3['c_ln']:.6e} +/- {res3['se_c_ln']:.6e})")
print(f"    R^2           = {res3['r2_ln']:.6f}")
print(f"  Time constant tau_s  = {res3['tau']:.3f} +/- {res3['delta_tau']:.3f} s  (Interpolated: {tau_interp3:.2f} s)")
print(f"  Experimental C_eq    = {C_series_exp*1e6:.3f} +/- {dC_series_fit*1e6:.3f} (fit) +/- {dC_series_tot*1e6:.3f} (tot) uF")
print(f"  Relative error       = {err_series_pct:.2f}%")

print("\n--- PART 4: Parallel Combination (C1 = 20 uF, C2 = 4 uF) ---")
print(f"  Theoretical C_eq     = {C_parallel_th*1e6:.3f} uF")
print(f"  Linear fit (ln V/V0 vs t):")
print(f"    Slope m       = ({res4['m_ln']:.6e} +/- {res4['se_m_ln']:.6e}) s^-1")
print(f"    Intercept c   = ({res4['c_ln']:.6e} +/- {res4['se_c_ln']:.6e})")
print(f"    R^2           = {res4['r2_ln']:.6f}")
print(f"  Time constant tau_p  = {res4['tau']:.3f} +/- {res4['delta_tau']:.3f} s  (Interpolated: {tau_interp4:.2f} s)")
print(f"  Experimental C_eq    = {C_parallel_exp*1e6:.3f} +/- {dC_parallel_fit*1e6:.3f} (fit) +/- {dC_parallel_tot*1e6:.3f} (tot) uF")
print(f"  Relative error       = {err_parallel_pct:.2f}%")
print("=" * 75)

# -----------------------------------------------------------------------------
# 6. Plotting Functions (PDF Output)
# -----------------------------------------------------------------------------
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.family': 'serif'
})

# --- Plot 1: Part 1 Charging C1 ---
fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(13, 5))

# Semi-log / Linear ln(V/V0)
ax1a.errorbar(t1, res1['ratio'], yerr=delta_V / V0_1, fmt='s', color='#1f77b4',
              capsize=4, label='Data: $V(t)/V_0$')
t_fine = np.linspace(0, 140, 200)
ax1a.plot(t_fine, np.exp(-t_fine / res1['tau']), 'r--',
          label=f'Fit: $\\tau_1 = ({res1["tau"]:.1f} \\pm {res1["delta_tau"]:.1f})$ s\n$R^2 = {res1["r2_ln"]:.5f}$')
ax1a.set_yscale('log')
ax1a.set_xlabel('Time $t$ [s]')
ax1a.set_ylabel('$V(t) / V_0$ (Log Scale)')
ax1a.set_title('(a) Semilogarithmic Decay of Voltmeter Voltage')
ax1a.grid(True, which='both', linestyle=':', alpha=0.6)
ax1a.legend()

# Voltage vs time
ax1b.errorbar(t1, V1, yerr=delta_V, fmt='o', color='#2ca02c', capsize=4, label='Measured $V(t)$')
ax1b.plot(t_fine, V0_1 * np.exp(-t_fine / res1['tau']), 'b-',
          label=f'$V(t) = {V0_1:.1f} e^{{-t / {res1["tau"]:.1f}}}$ V')
# Also show capacitor voltage V_C(t) = V0 * (1 - exp(-t/tau))
ax1b.plot(t_fine, V0_1 * (1.0 - np.exp(-t_fine / res1['tau'])), 'm:',
          label='Capacitor Voltage $V_C(t) = V_0(1 - e^{-t/\\tau})$')
ax1b.set_xlabel('Time $t$ [s]')
ax1b.set_ylabel('Voltage [V]')
ax1b.set_title('(b) Voltmeter Reading & Capacitor Voltage')
ax1b.grid(True, linestyle='--', alpha=0.6)
ax1b.legend()

fig1.suptitle('Part 1: Charging of Capacitor $C_1 = 20\\,\\mu\\mathrm{F}$ through Voltmeter $R_v$', y=1.02)
fig1.tight_layout()
fig1.savefig(plots_dir / 'part1_C1_charging.pdf', bbox_inches='tight')
plt.close(fig1)
print("  -> Saved plots/part1_C1_charging.pdf")

# --- Plot 2: Part 2 Discharging C2 ---
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(13, 5))

ax2a.errorbar(t2, res2['ratio'], yerr=delta_V / V0_2, fmt='s', color='#ff7f0e',
              capsize=4, label='Data: $V(t)/V_0$')
ax2a.plot(t_fine, np.exp(-t_fine / res2['tau']), 'r--',
          label=f'Fit: $\\tau_2 = ({res2["tau"]:.2f} \\pm {res2["delta_tau"]:.2f})$ s\n$R^2 = {res2["r2_ln"]:.5f}$')
ax2a.set_yscale('log')
ax2a.set_xlabel('Time $t$ [s]')
ax2a.set_ylabel('$V(t) / V_0$ (Log Scale)')
ax2a.set_title('(a) Semilogarithmic Discharge')
ax2a.grid(True, which='both', linestyle=':', alpha=0.6)
ax2a.legend()

ax2b.errorbar(t2, V2, yerr=delta_V, fmt='o', color='#d62728', capsize=4, label='Measured $V(t)$')
ax2b.plot(t_fine, V0_2 * np.exp(-t_fine / res2['tau']), 'k-',
          label=f'$V(t) = {V0_2:.1f} e^{{-t / {res2["tau"]:.1f}}}$ V')
ax2b.axhline(V0_2 / np.e, color='blue', linestyle='--', alpha=0.7, label=f'$V_0/e = {V0_2/np.e:.2f}$ V')
ax2b.axvline(res2['tau'], color='blue', linestyle=':', alpha=0.7, label=f'$\\tau_2 = {res2["tau"]:.1f}$ s')
ax2b.set_xlabel('Time $t$ [s]')
ax2b.set_ylabel('Voltage $V$ [V]')
ax2b.set_title('(b) Exponential Discharge Curve')
ax2b.grid(True, linestyle='--', alpha=0.6)
ax2b.legend()

fig2.suptitle('Part 2: Discharging of Capacitor $C_2 = 4\\,\\mu\\mathrm{F}$ through Voltmeter $R_v$', y=1.02)
fig2.tight_layout()
fig2.savefig(plots_dir / 'part2_C2_discharging.pdf', bbox_inches='tight')
plt.close(fig2)
print("  -> Saved plots/part2_C2_discharging.pdf")

# --- Plot 3: Part 3 Series Discharging ---
fig3, ax3 = plt.subplots(figsize=(7, 5))
ax3.errorbar(t3, V3, yerr=delta_V, fmt='s', color='#9467bd', capsize=4, label='Experimental data')
ax3.plot(t_fine, V0_3 * np.exp(-t_fine / res3['tau']), 'r-',
         label=f'OLS Fit: $\\tau_s = ({res3["tau"]:.2f} \\pm {res3["delta_tau"]:.2f})$ s\n$C_{{eq}} = ({C_series_exp*1e6:.2f} \\pm {dC_series_tot*1e6:.2f})\\,\\mu$F')
ax3.axhline(V0_3 / np.e, color='grey', linestyle='--', label=f'$V_0/e = {V0_3/np.e:.2f}$ V')
ax3.axvline(res3['tau'], color='grey', linestyle=':', label=f'$\\tau = {res3["tau"]:.1f}$ s')
ax3.set_xlabel('Time $t$ [s]')
ax3.set_ylabel('Voltage $V(t)$ [V]')
ax3.set_title('Part 3: Series Combination Discharge ($C_1$ & $C_2$ in Series)')
ax3.grid(True, linestyle='--', alpha=0.6)
ax3.legend()
fig3.tight_layout()
fig3.savefig(plots_dir / 'part3_series_discharging.pdf', bbox_inches='tight')
plt.close(fig3)
print("  -> Saved plots/part3_series_discharging.pdf")

# --- Plot 4: Part 4 Parallel Discharging ---
fig4, ax4 = plt.subplots(figsize=(7, 5))
t4_fine = np.linspace(0, 310, 200)
ax4.errorbar(t4, V4, yerr=delta_V, fmt='s', color='#17becf', capsize=4, label='Experimental data')
ax4.plot(t4_fine, V0_4 * np.exp(-t4_fine / res4['tau']), 'r-',
         label=f'OLS Fit: $\\tau_p = ({res4["tau"]:.1f} \\pm {res4["delta_tau"]:.1f})$ s\n$C_{{eq}} = ({C_parallel_exp*1e6:.2f} \\pm {dC_parallel_tot*1e6:.2f})\\,\\mu$F')
ax4.axhline(V0_4 / np.e, color='grey', linestyle='--', label=f'$V_0/e = {V0_4/np.e:.2f}$ V')
ax4.axvline(res4['tau'], color='grey', linestyle=':', label=f'$\\tau = {res4["tau"]:.1f}$ s')
ax4.set_xlabel('Time $t$ [s]')
ax4.set_ylabel('Voltage $V(t)$ [V]')
ax4.set_title('Part 4: Parallel Combination Discharge ($C_1$ & $C_2$ in Parallel)')
ax4.grid(True, linestyle='--', alpha=0.6)
ax4.legend()
fig4.tight_layout()
fig4.savefig(plots_dir / 'part4_parallel_discharging.pdf', bbox_inches='tight')
plt.close(fig4)
print("  -> Saved plots/part4_parallel_discharging.pdf")

# --- Plot 5: All 4 Regressions Comparison (log10(V/V0) vs t) ---
fig5, ax5 = plt.subplots(figsize=(8, 5.5))
ax5.plot(t1, res1['y_log10'], 's-', color='#1f77b4',
         label=f'$C_1$ ($20\\,\\mu$F): slope = {res1["m_log10"]*1e3:.2f} $\\times 10^{{-3}}\\,\\mathrm{{s}}^{{-1}}$, $R^2={res1["r2_log10"]:.5f}$')
ax5.plot(t2, res2['y_log10'], 'o-', color='#ff7f0e',
         label=f'$C_2$ ($4\\,\\mu$F): slope = {res2["m_log10"]*1e3:.2f} $\\times 10^{{-3}}\\,\\mathrm{{s}}^{{-1}}$, $R^2={res2["r2_log10"]:.5f}$')
ax5.plot(t3, res3['y_log10'], '^-', color='#9467bd',
         label=f'Series: slope = {res3["m_log10"]*1e3:.2f} $\\times 10^{{-3}}\\,\\mathrm{{s}}^{{-1}}$, $R^2={res3["r2_log10"]:.5f}$')
ax5.plot(t4, res4['y_log10'], 'd-', color='#2ca02c',
         label=f'Parallel: slope = {res4["m_log10"]*1e3:.2f} $\\times 10^{{-3}}\\,\\mathrm{{s}}^{{-1}}$, $R^2={res4["r2_log10"]:.5f}$')

ax5.set_xlabel('Time $t$ [s]')
ax5.set_ylabel('$\\log_{10}(V / V_0)$')
ax5.set_title('Linear Regression Comparison: $\\log_{10}(V/V_0)$ vs. $t$ for All Circuits')
ax5.grid(True, linestyle='--', alpha=0.6)
ax5.legend(loc='lower left')
fig5.tight_layout()
fig5.savefig(plots_dir / 'all_circuits_log_comparison.pdf', bbox_inches='tight')
plt.close(fig5)
print("  -> Saved plots/all_circuits_log_comparison.pdf")

# --- Plot 6: Equivalent Capacitance Comparison (Theory vs Experiment) ---
fig6, ax6 = plt.subplots(figsize=(7, 5))
labels = ['Series ($C_1 \\| C_2$)', 'Parallel ($C_1 + C_2$)']
th_vals = [C_series_th * 1e6, C_parallel_th * 1e6]
exp_vals = [C_series_exp * 1e6, C_parallel_exp * 1e6]
exp_errs = [dC_series_tot * 1e6, dC_parallel_tot * 1e6]

x_pos = np.arange(len(labels))
width = 0.35

rects1 = ax6.bar(x_pos - width/2, th_vals, width, label='Theoretical', color='#3498db', edgecolor='black')
rects2 = ax6.bar(x_pos + width/2, exp_vals, width, yerr=exp_errs, capsize=6,
                 label='Experimental (from $\\tau / \\bar{R}_v$)', color='#e74c3c', edgecolor='black')

ax6.set_ylabel('Equivalent Capacitance [$\\mu$F]')
ax6.set_title('Comparison of Theoretical vs. Experimental Capacitances')
ax6.set_xticks(x_pos)
ax6.set_xticklabels(labels)
ax6.legend()
ax6.grid(True, axis='y', linestyle='--', alpha=0.6)

# Add values above bars
for rect in rects1:
    h = rect.get_height()
    ax6.annotate(f'{h:.2f} $\\mu$F', xy=(rect.get_x() + rect.get_width()/2, h),
                 xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
for rect in rects2:
    h = rect.get_height()
    ax6.annotate(f'{h:.2f} $\\mu$F', xy=(rect.get_x() + rect.get_width()/2, h),
                 xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

fig6.tight_layout()
fig6.savefig(plots_dir / 'capacitance_comparison.pdf', bbox_inches='tight')
plt.close(fig6)
print("  -> Saved plots/capacitance_comparison.pdf")

print("\nAll analyses and plots successfully generated!")
