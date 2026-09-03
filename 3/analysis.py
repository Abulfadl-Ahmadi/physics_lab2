"""
Experiment 3: Magnetic Field Analysis & Earth's Magnetic Field
Physics Laboratory II - Sharif University of Technology

Analyses:
1. Axial magnetic field profile of Helmholtz coils: experimental vs theoretical.
2. Linear regression of B_H vs I to determine geometric/coil constant K and vacuum permeability mu_0.
3. Linear regression of tangent galvanometer data (I*K vs tan(alpha)) to measure Earth's horizontal magnetic field B_E^h.
4. Calculation of vertical component B_E^v and total Earth's magnetic field |B_E| with rigorous error propagation.
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ensure plots directory exists
plots_dir = Path(__file__).parent / 'plots'
plots_dir.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# OLS Regression Helper (Standard Lab-IV Pattern)
# -----------------------------------------------------------------------------
def ols_fit(x, y):
    """
    Ordinary Least Squares linear regression y = m*x + c.
    Returns:
        m      : slope
        se_m   : standard error of slope
        c      : intercept
        se_c   : standard error of intercept
        r2     : coefficient of determination R^2
    """
    N = len(x)
    m, c = np.polyfit(x, y, 1)
    y_pred = m * x + c
    residuals = y - y_pred
    ss_res = np.sum(residuals**2)
    s_yx = np.sqrt(ss_res / (N - 2))
    ss_x = np.sum((x - np.mean(x))**2)
    se_m = s_yx / np.sqrt(ss_x)
    se_c = s_yx * np.sqrt(1.0 / N + np.mean(x)**2 / ss_x)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (ss_res / ss_tot)
    return m, se_m, c, se_c, r2

# =============================================================================
# PART 1: Axial Magnetic Field Profile (B vs D)
# =============================================================================
# Table 1: Axial displacement D [cm] and measured magnetic field B_H [mT]
D_cm = np.array([0.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0])
B_H_mT = np.array([1.49, 1.42, 1.26, 1.00, 0.73, 0.50, 0.35, 0.24, 0.18, 0.13])

delta_D_cm = 0.2     # uncertainty in position reading [cm]
delta_B_mT = 0.01    # uncertainty in teslameter reading [mT]

# Theoretical model comparison:
# Coils: radius R = 0.20 m (20 cm).
# Case A: Helmholtz separation a = R = 0.20 m
# Case B: Narrow separation a = R/2 = 0.10 m
R_val = 0.20
D_cont = np.linspace(0, 0.55, 300)

def B_axial_norm(D, a, R):
    """Normalized axial magnetic field for two parallel circular coils of radius R separated by distance a."""
    term1 = 1.0 / (R**2 + (D - a / 2.0)**2)**1.5
    term2 = 1.0 / (R**2 + (D + a / 2.0)**2)**1.5
    return term1 + term2

# Theoretical shapes scaled to match central field B_H(0) = 1.49 mT
B_theo_Helm = B_axial_norm(D_cont, a=R_val, R=R_val)
B_theo_Helm = B_theo_Helm / B_axial_norm(0.0, a=R_val, R=R_val) * B_H_mT[0]

B_theo_narrow = B_axial_norm(D_cont, a=R_val/2.0, R=R_val)
B_theo_narrow = B_theo_narrow / B_axial_norm(0.0, a=R_val/2.0, R=R_val) * B_H_mT[0]

# Plot 1: B vs D axial profile
plt.figure(figsize=(7, 4.8), dpi=300)
plt.errorbar(D_cm, B_H_mT, xerr=delta_D_cm, yerr=delta_B_mT, fmt='s',
             color='#1f77b4', ecolor='#1f77b4', elinewidth=1.2, capsize=3,
             label=r'Experimental Data ($B_H \pm 0.01\,\mathrm{mT}$)')
plt.plot(D_cont * 100, B_theo_Helm, '--', color='#2ca02c', lw=1.8,
         label=r'Theory Helmholtz ($a = R = 20\,\mathrm{cm}$)')
plt.plot(D_cont * 100, B_theo_narrow, ':', color='#d62728', lw=1.8,
         label=r'Theory Narrow ($a = R/2 = 10\,\mathrm{cm}$)')

plt.title('Axial Magnetic Field Distribution of Helmholtz Coils', fontsize=12, pad=10)
plt.xlabel('Displacement along Axis $D$ (cm)', fontsize=11)
plt.ylabel('Magnetic Field $B_H$ (mT)', fontsize=11)
plt.xlim(-1, 55)
plt.ylim(0, 1.7)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
plt.tight_layout()
plt.savefig(plots_dir / 'B_vs_D_axial_profile.pdf')
plt.close()

# =============================================================================
# PART 2: Determination of Permeability of Free Space mu_0 (B vs I)
# =============================================================================
# Table 2: Current I [A] and Magnetic field B_H [mT]
I_A = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.23, 1.43, 1.65, 1.77])
B_mT = np.array([0.13, 0.28, 0.45, 0.61, 0.79, 1.01, 1.18, 1.35, 1.44])

delta_I_A = 0.01     # uncertainty in current meter [A]
delta_B_mu0 = 0.01   # uncertainty in field [mT]

# Coil parameters
n_turns = 154         # turns per coil
R_m = 0.20           # coil radius in meters (20 cm)
delta_R_m = 0.002     # uncertainty in radius (2 mm)

# Perform OLS fit: B_mT = m * I_A + c
m_B, se_m_B, c_B, se_c_B, r2_B = ols_fit(I_A, B_mT)

# Slope in SI units [T/A]
K_exp = m_B * 1e-3
se_K = se_m_B * 1e-3

# Determination of mu_0:
# Theoretical relation for Helmholtz center field:
# B = (8 / (5 * sqrt(5))) * (mu_0 * n / R) * I = K * I
# => mu_0 = K * R * 5 * sqrt(5) / (8 * n)
geom_factor = (5.0 * np.sqrt(5.0)) / (8.0 * n_turns)
mu_0_exp = K_exp * R_m * geom_factor

# Error propagation for mu_0:
# mu_0 = f(K, R) => (d_mu0 / mu0)^2 = (se_K / K)^2 + (delta_R / R)^2
rel_err_K = se_K / K_exp
rel_err_R = delta_R_m / R_m
se_mu_0 = mu_0_exp * np.sqrt(rel_err_K**2 + rel_err_R**2)

mu_0_theory = 4.0 * np.pi * 1e-7
rel_discrepancy = abs(mu_0_exp - mu_0_theory) / mu_0_theory * 100.0

# Plot 2: B vs I linear regression
plt.figure(figsize=(7, 4.8), dpi=300)
plt.errorbar(I_A, B_mT, xerr=delta_I_A, yerr=delta_B_mu0, fmt='o',
             color='#d62728', ecolor='#d62728', elinewidth=1.2, capsize=3,
             label=r'Experimental Data ($B_H \pm 0.01\,\mathrm{mT}$)')

I_fit = np.linspace(0, 2.0, 200)
B_fit = m_B * I_fit + c_B
plt.plot(I_fit, B_fit, '-', color='#1f77b4', lw=1.8,
         label=f'OLS Fit: $B_H = ({m_B:.3f} \\pm {se_m_B:.3f})I + ({c_B:.3f} \\pm {se_c_B:.3f})$\n$R^2 = {r2_B:.4f}$')

plt.title(r'Magnetic Field vs Coil Current ($B_H$ vs $I$)', fontsize=12, pad=10)
plt.xlabel('Current $I$ (A)', fontsize=11)
plt.ylabel('Magnetic Field $B_H$ (mT)', fontsize=11)
plt.xlim(0, 2.0)
plt.ylim(0, 1.7)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5, loc='upper left')
plt.tight_layout()
plt.savefig(plots_dir / 'B_vs_I_linear_fit.pdf')
plt.close()

# =============================================================================
# PART 3: Earth's Horizontal Magnetic Field B_E^h (Tangent Galvanometer)
# =============================================================================
# Table 3: Current I [mA], Deflection angle alpha [deg], tan(alpha)
I_galv_mA = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
alpha_deg = np.array([14.0, 31.0, 42.0, 51.0, 60.0, 64.0, 66.0, 69.0, 70.0, 74.0])
tan_alpha = np.tan(np.radians(alpha_deg))

delta_alpha_deg = 1.0  # uncertainty in compass scale angle [deg]
# tan(alpha) uncertainty: d(tan a)/da = sec^2(a) * da_rad
delta_alpha_rad = np.radians(delta_alpha_deg)
delta_tan_alpha = (1.0 / np.cos(np.radians(alpha_deg))**2) * delta_alpha_rad

# Field generated by coil: B_H = K * I [Tesla]
# Using experimental K from Part 2
I_galv_A = I_galv_mA * 1e-3
y_IK = I_galv_A * K_exp  # [Tesla]

# Uncertainty in y_IK: d(IK) = IK * sqrt((dI/I)^2 + (dK/K)^2)
delta_I_galv_A = 1.0 * 1e-3  # 1 mA
delta_y_IK = y_IK * np.sqrt((delta_I_galv_A / I_galv_A)**2 + (se_K / K_exp)**2)

# OLS Fit: y_IK = m_E * tan_alpha + c_E
# Theory: y_IK = B_E^h * tan(alpha) => slope m_E is horizontal Earth field B_E^h
m_E, se_m_E, c_E, se_c_E, r2_E = ols_fit(tan_alpha, y_IK)
B_Eh_exp = m_E
se_B_Eh = se_m_E

# Plot 3: IK vs tan(alpha)
plt.figure(figsize=(7, 4.8), dpi=300)
plt.errorbar(tan_alpha, y_IK * 1e6, xerr=delta_tan_alpha, yerr=delta_y_IK * 1e6, fmt='s',
             color='#2ca02c', ecolor='#2ca02c', elinewidth=1.2, capsize=3,
             label=r'Experimental Data ($I\cdot K$ vs $\tan\alpha$)')

tan_fit = np.linspace(0, 3.8, 200)
IK_fit = (m_E * tan_fit + c_E) * 1e6
plt.plot(tan_fit, IK_fit, '-', color='#9467bd', lw=1.8,
         label=f'OLS Fit: $IK = ({m_E*1e6:.2f} \\pm {se_m_E*1e6:.2f})\\tan\\alpha + ({c_E*1e6:.2f} \\pm {se_c_E*1e6:.2f})$\n$R^2 = {r2_E:.4f}$')

plt.title(r'Tangent Galvanometer: $I\cdot K$ vs $\tan\alpha$', fontsize=12, pad=10)
plt.xlabel(r'$\tan\alpha$', fontsize=11)
plt.ylabel(r'$I\cdot K$ ($\mu\mathrm{T}$)', fontsize=11)
plt.xlim(0, 3.8)
plt.ylim(0, 100)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5, loc='upper left')
plt.tight_layout()
plt.savefig(plots_dir / 'earth_field_tangent_galvanometer.pdf')
plt.close()

# =============================================================================
# PART 4: Vertical Component and Total Earth's Field
# =============================================================================
# Dip angle nu
nu_deg = 55.0
delta_nu_deg = 1.0
nu_rad = np.radians(nu_deg)
delta_nu_rad = np.radians(delta_nu_deg)

# Vertical component: B_E^v = B_E^h * tan(nu)
B_Ev_exp = B_Eh_exp * np.tan(nu_rad)

# Error propagation for B_E^v:
# d(B_Ev) = sqrt( (tan(nu) * d_BEh)^2 + (B_Eh * sec^2(nu) * d_nu)^2 )
dB_Ev_dBEh = np.tan(nu_rad)
dB_Ev_dnu = B_Eh_exp * (1.0 / np.cos(nu_rad)**2)
se_B_Ev = np.sqrt((dB_Ev_dBEh * se_B_Eh)**2 + (dB_Ev_dnu * delta_nu_rad)**2)

# Total Earth Field: |B_E| = sqrt((B_Eh)^2 + (B_Ev)^2)
B_E_total = np.sqrt(B_Eh_exp**2 + B_Ev_exp**2)

# Error propagation for |B_E|:
# d|B_E| = (1 / |B_E|) * sqrt( (B_Eh * d_BEh)^2 + (B_Ev * d_BEv)^2 )
se_B_E_total = (1.0 / B_E_total) * np.sqrt((B_Eh_exp * se_B_Eh)**2 + (B_Ev_exp * se_B_Ev)**2)

# Print Summary of Results
print("=" * 65)
print("EXPERIMENT 3: MAGNETIC FIELD & HELMHOLTZ COIL ANALYSIS RESULTS")
print("=" * 65)
print("\n--- Part 2: Vacuum Permeability (mu_0) ---")
print(f"OLS Slope m (K)      = {m_B:.4f} +/- {se_m_B:.4f} mT/A  ({K_exp:.4e} +/- {se_K:.4e} T/A)")
print(f"OLS Intercept c      = {c_B:.4f} +/- {se_c_B:.4f} mT")
print(f"Fit R^2              = {r2_B:.6f}")
print(f"Experimental mu_0    = ({mu_0_exp*1e6:.3f} +/- {se_mu_0*1e6:.3f}) x 10^-6 H/m")
print(f"Theoretical mu_0     = {mu_0_theory*1e6:.3f} x 10^-6 H/m")
print(f"Relative Discrepancy = {rel_discrepancy:.2f}%")

print("\n--- Part 3: Horizontal Component of Earth's Field (B_E^h) ---")
print(f"OLS Slope (B_E^h)    = ({B_Eh_exp*1e6:.2f} +/- {se_B_Eh*1e6:.2f}) uT")
print(f"OLS Intercept c      = ({c_E*1e6:.2f} +/- {se_c_E*1e6:.2f}) uT")
print(f"Fit R^2              = {r2_E:.6f}")

print("\n--- Part 4: Vertical Component & Total Earth Magnetic Field ---")
print(f"Dip Angle (nu)       = {nu_deg:.1f} +/- {delta_nu_deg:.1f} deg")
print(f"Vertical Field B_E^v = ({B_Ev_exp*1e6:.2f} +/- {se_B_Ev*1e6:.2f}) uT")
print(f"Total Field |B_E|    = ({B_E_total*1e6:.2f} +/- {se_B_E_total*1e6:.2f}) uT")
print("=" * 65)
