import numpy as np
import matplotlib.pyplot as plt

R = 0.2
a = R/2
D = np.linspace(0,0.5,200)

def B_func(D):
    return ( 1/((R**2 + (D - a/2)**2)**1.5) +
             1/((R**2 + (D + a/2)**2)**1.5) )

plt.figure(figsize=(5,3.5))
plt.plot(D*100, B_func(D))
plt.xlabel("D (cm)")
plt.ylabel("B (arb. units)")
plt.title("Theoretical B(D) for a = R/2")
plt.grid()
plt.savefig("fig-01.pdf")
# plt.show()
