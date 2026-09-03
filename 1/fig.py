import numpy as np
import matplotlib.pyplot as plt

# enabling LaTeX for mtplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})

# exp. data
I = np.array([100, 200, 300, 400, 500])
V = np.array([0.055, 0.125, 0.156, 0.233, 0.265])

# error bars
# dI = np.full_like(I, 5)
# dV = np.full_like(V, 0.005)

# fit line
a, b = np.polyfit(I, V, 1)

# ploting
plt.figure(figsize=(6, 4))
# plt.errorbar(I, V, xerr=dI, yerr=dV, fmt="bo", label="Experiment Data", capsize=4)

I_fit = np.linspace(0, 500, 100)
V_fit = a * I_fit + b

plt.plot(I, V, "bo", label="Experiment Data")
plt.plot(I_fit, V_fit, "r-", label=f"V = {a:.4f} I + {b:.4f}")
plt.xlabel("Current [mA]")
plt.ylabel("Voltage [V]")
plt.title("Graph of Voltage as a Functions of Current")
plt.legend()
plt.grid(True)

# saving as PDF
plt.savefig("ohm_law_plot.pdf", format="pdf")

# plt.show()
