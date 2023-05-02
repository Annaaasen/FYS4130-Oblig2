import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

#Setting T = T/J and k_B=1

N = 16
Ts = np.array([0.25, 0.5])
r = np.linspace(0, 16, 100)

#####################
### Analytic part ###
####################
def analytic_C(N, r, T):
    beta = 1/T
    a = np.exp(beta) + 2
    b = np.exp(beta) - 1
    nom = a**r * b**(N-r) + b**r * a**(N-r) + b**N
    denom = 2*b**N + a**N
    return nom / denom 


for T in Ts:
    plt.plot(r, analytic_C(N, r, T))
plt.show()



