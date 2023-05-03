import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

#Setting T = T/J and k_B=1

mag = pd.read_csv("16 magnetisations_vs_T.csv")
mag.plot(x="T", y="m")
plt.show()


mag = pd.read_csv("16 magnetisations_vs_T_20.csv", header=None, names=["T", "m", "m1", "m2", "m4"])
mag.plot(x="T", y="m")
plt.show()