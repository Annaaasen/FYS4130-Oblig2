import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import plot_utils

#Setting T = T/J and k_B=1
mag8 = pd.read_csv("8 magnetisations_vs_T_2D.csv", header=None, names=["T", "m", "m1", "m2", "m4"])
gamma8 = mag8["m4"] / mag8["m2"]**2
plt.plot(mag8["T"], gamma8, label="L=8")

mag16 = pd.read_csv("16 magnetisations_vs_T_2D.csv", header=None, names=["T", "m", "m1", "m2", "m4"])
gamma16 = mag16["m4"] / mag16["m2"]**2
plt.plot(mag16["T"], gamma16, label="L=16")

# mag24 = pd.read_csv("24 magnetisations_vs_T.csv", header=None, names=["T", "m", "m1", "m2", "m4"])
# gamma24 = mag24["m4"] / mag24["m2"]**2
# plt.plot(mag24["T"], gamma24)

mag32 = pd.read_csv("32 magnetisations_vs_T_2D.csv", header=None, names=["T", "m", "m1", "m2", "m4"])
gamma32 = mag32["m4"] / mag32["m2"]**2
plt.plot(mag32["T"], gamma32, label="L=32")

plt.xlabel("T/J")
plt.ylabel(r"$\Gamma$")
plt.legend()
plt.savefig("../tex/figs/gamma_T.pdf")
plt.show()


plt.plot(mag8["T"], gamma8, label="L=8")
plt.plot(mag16["T"], gamma16, label="L=16")
plt.plot(mag32["T"], gamma32, label="L=32")
plt.xlim((0.985, 0.99))
plt.ylim((1.12, 1.13))
plt.xlabel("T/J")
plt.ylabel(r"$\Gamma$")
plt.legend()
plt.savefig("../tex/figs/gamma_T_zoom.pdf")
plt.show()


# print(gamma8[idx])
# print(idx, mag16["T"][idx], gamma16[idx])
