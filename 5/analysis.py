"""
Experiment 5 - Force on a Current-Carrying Wire in a Magnetic Field
Lab-Phy-II, Sharif University of Technology

Analyses:
  Part A  - DeltaF vs L  (i=4 A, Im=2 A)   -> B from slope
  Part B  - DeltaF vs i  (L=10 cm, Im=2 A) -> B from slope
  Part C  - DeltaF vs Im (L=10 cm, i=4 A)  -> alpha = B/Im conversion factor
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

Path('plots').mkdir(exist_ok=True)

delta_F_mN  = 0.05
delta_L_cm  = 0.05
delta_i_A   = 0.05
delta_Im_A  = 0.05


def ols_fit(x, y):
    N = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals ** 2)
    s_yx   = np.sqrt(ss_res / (N - 2))
    ss_x   = np.sum((x - np.mean(x)) ** 2)
    se_m   = s_yx / np.sqrt(ss_x)
    se_c   = s_yx * np.sqrt(1 / N + np.mean(x) ** 2 / ss_x)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2     = 1 - ss_res / ss_tot
    return m, se_m, c, se_c, r2


# PART A - DeltaF vs L
L_cm   = np.array([1.25, 2.50, 5.00, 10.00])
F_A_mN = np.array([0.95, 1.68, 2.74,  5.16])
L_m    = L_cm  * 1e-2
F_A_N  = F_A_mN * 1e-3
i_A_part = 4.0

slope_A, se_slope_A, intcpt_A, se_intcpt_A, r2_A = ols_fit(L_m, F_A_N)
B_A  = slope_A / i_A_part
dB_A = B_A * np.sqrt((se_slope_A / slope_A) ** 2 + (delta_i_A / i_A_part) ** 2)

print("=" * 60)
print("PART A - DeltaF vs L  (i = 4 A, Im = 2 A)")
print("=" * 60)
print(f"  slope     = {slope_A:.6f} +/- {se_slope_A:.6f}  N/m")
print(f"  intercept = {intcpt_A*1e3:.4f} +/- {se_intcpt_A*1e3:.4f}  mN")
print(f"  R2        = {r2_A:.6f}")
print(f"  B         = {B_A*1e3:.3f} +/- {dB_A*1e3:.3f}  mT")
print()

fig_A, ax_A = plt.subplots(figsize=(7, 5))
ax_A.errorbar(L_cm, F_A_mN, xerr=delta_L_cm, yerr=delta_F_mN,
              fmt='s', color='steelblue', capsize=4, label='Experimental data')
L_fit = np.linspace(0, 11, 200)
ax_A.plot(L_fit, (slope_A * L_fit * 1e-2 + intcpt_A) * 1e3, 'r-',
          label=f'OLS fit ($R^2={r2_A:.4f}$)\n$B = ({B_A*1e3:.2f} \\pm {dB_A*1e3:.2f})$ mT')
ax_A.set_xlabel('Wire length  L  [cm]', fontsize=12)
ax_A.set_ylabel('Magnetic force  dF  [mN]', fontsize=12)
ax_A.set_title('DeltaF vs L  (i=4A, Im=2A)', fontsize=13)
ax_A.legend(fontsize=10)
ax_A.grid(True, linestyle='--', alpha=0.5)
ax_A.set_xlim(0, 11); ax_A.set_ylim(0, 6)
fig_A.tight_layout()
fig_A.savefig('plots/F_vs_L.pdf')
plt.close(fig_A)
print("  -> plots/F_vs_L.pdf saved")
print()


# PART B - DeltaF vs i
i_B_A  = np.array([1.00, 2.00, 3.00, 4.00])
F_B_mN = np.array([0.23, 0.46, 1.91, 3.16])
F_B_N  = F_B_mN * 1e-3
L_B_m  = 0.10
delta_L_m = delta_L_cm * 1e-2

slope_B, se_slope_B, intcpt_B, se_intcpt_B, r2_B = ols_fit(i_B_A, F_B_N)
B_B  = slope_B / L_B_m
dB_B = B_B * np.sqrt((se_slope_B / slope_B) ** 2 + (delta_L_m / L_B_m) ** 2)

print("=" * 60)
print("PART B - DeltaF vs i  (L = 10 cm, Im = 2 A)")
print("=" * 60)
print(f"  slope     = {slope_B*1e3:.6f} +/- {se_slope_B*1e3:.6f}  mN/A")
print(f"  intercept = {intcpt_B*1e3:.4f} +/- {se_intcpt_B*1e3:.4f}  mN")
print(f"  R2        = {r2_B:.6f}")
print(f"  B         = {B_B*1e3:.3f} +/- {dB_B*1e3:.3f}  mT")
print()

fig_B, ax_B = plt.subplots(figsize=(7, 5))
ax_B.errorbar(i_B_A, F_B_mN, xerr=delta_i_A, yerr=delta_F_mN,
              fmt='s', color='darkorange', capsize=4, label='Experimental data')
i_fit = np.linspace(0, 4.5, 200)
ax_B.plot(i_fit, (slope_B * i_fit + intcpt_B) * 1e3, 'r-',
          label=f'OLS fit ($R^2={r2_B:.4f}$)\n$B = ({B_B*1e3:.2f} \\pm {dB_B*1e3:.2f})$ mT')
ax_B.set_xlabel('Loop current  i  [A]', fontsize=12)
ax_B.set_ylabel('Magnetic force  dF  [mN]', fontsize=12)
ax_B.set_title('DeltaF vs i  (L=10cm, Im=2A)', fontsize=13)
ax_B.legend(fontsize=10)
ax_B.grid(True, linestyle='--', alpha=0.5)
ax_B.set_xlim(0, 4.5)
fig_B.tight_layout()
fig_B.savefig('plots/F_vs_i.pdf')
plt.close(fig_B)
print("  -> plots/F_vs_i.pdf saved")
print()


# PART C - DeltaF vs Im
Im_C_A = np.array([0.50, 1.00, 1.50, 2.00])
F_C_mN = np.array([0.32, 0.66, 2.00, 3.16])
F_C_N  = F_C_mN * 1e-3
i_C_A  = 4.0
L_C_m  = 0.10

slope_C, se_slope_C, intcpt_C, se_intcpt_C, r2_C = ols_fit(Im_C_A, F_C_N)
alpha_C  = slope_C / (i_C_A * L_C_m)
dalpha_C = alpha_C * np.sqrt((se_slope_C / slope_C) ** 2
                             + (delta_i_A / i_C_A) ** 2
                             + (delta_L_m / L_C_m) ** 2)
B_C_at2  = alpha_C  * 2.0
dB_C_at2 = dalpha_C * 2.0

print("=" * 60)
print("PART C - DeltaF vs Im  (L = 10 cm, i = 4 A)")
print("=" * 60)
print(f"  slope     = {slope_C*1e3:.6f} +/- {se_slope_C*1e3:.6f}  mN/A")
print(f"  intercept = {intcpt_C*1e3:.4f} +/- {se_intcpt_C*1e3:.4f}  mN")
print(f"  R2        = {r2_C:.6f}")
print(f"  alpha     = {alpha_C*1e3:.4f} +/- {dalpha_C*1e3:.4f}  mT/A")
print(f"  B(Im=2 A) = {B_C_at2*1e3:.3f} +/- {dB_C_at2*1e3:.3f}  mT")
print()

fig_C, ax_C = plt.subplots(figsize=(7, 5))
ax_C.errorbar(Im_C_A, F_C_mN, xerr=delta_Im_A, yerr=delta_F_mN,
              fmt='s', color='forestgreen', capsize=4, label='Experimental data')
Im_fit = np.linspace(0, 2.2, 200)
ax_C.plot(Im_fit, (slope_C * Im_fit + intcpt_C) * 1e3, 'r-',
          label=f'OLS fit ($R^2={r2_C:.4f}$)\n$\\alpha = ({alpha_C*1e3:.2f} \\pm {dalpha_C*1e3:.2f})$ mT/A')
ax_C.set_xlabel('Magnet coil current  Im  [A]', fontsize=12)
ax_C.set_ylabel('Magnetic force  dF  [mN]', fontsize=12)
ax_C.set_title('DeltaF vs Im  (L=10cm, i=4A)', fontsize=13)
ax_C.legend(fontsize=10)
ax_C.grid(True, linestyle='--', alpha=0.5)
ax_C.set_xlim(0, 2.2)
fig_C.tight_layout()
fig_C.savefig('plots/F_vs_Im.pdf')
plt.close(fig_C)
print("  -> plots/F_vs_Im.pdf saved")
print()


# SUMMARY
print("=" * 60)
print("SUMMARY - Magnetic Field B at Im = 2 A")
print("=" * 60)
print(f"  Part A (F vs L):  B = {B_A*1e3:.2f} +/- {dB_A*1e3:.2f}  mT")
print(f"  Part B (F vs i):  B = {B_B*1e3:.2f} +/- {dB_B*1e3:.2f}  mT")
print(f"  Part C (F vs Im): B = {B_C_at2*1e3:.2f} +/- {dB_C_at2*1e3:.2f}  mT")
print(f"  alpha (B/Im)    = {alpha_C*1e3:.4f} +/- {dalpha_C*1e3:.4f}  mT/A")
B_vals  = np.array([B_A, B_B, B_C_at2])
dB_vals = np.array([dB_A, dB_B, dB_C_at2])
weights = 1.0 / dB_vals ** 2
B_avg   = np.sum(weights * B_vals) / np.sum(weights)
dB_avg  = 1.0 / np.sqrt(np.sum(weights))
print(f"  Weighted average: B = {B_avg*1e3:.2f} +/- {dB_avg*1e3:.2f}  mT")
print("=" * 60)
