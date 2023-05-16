import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import plot_utils

#Setting T = T/J and k_B=1

N = 16
Ts = np.array([0.25, 0.5])
r = np.linspace(0, N, 100)

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


#####################
### Numeric part ###
####################
corr_025 = pd.read_csv("correlation,3,16,0.25.csv", header=None, names=["r", "corr_func"])
corr_05 =  pd.read_csv("correlation,3,16,0.50.csv", header=None, names=["r", "corr_func"])
corrs = np.array([corr_025["corr_func"], corr_05["corr_func"]])
num_rs = np.array([corr_025["r"], corr_05["r"]])



#####################
### Plotting part ###
####################
ana_labels = ["Exact T=0.25", "Exact T=0.5"]
num_labels = ["Num T=0.25", "Num T=0.5"]
for T, corr, num_r, ana_label, num_label in zip(Ts, corrs, num_rs, ana_labels, num_labels):
    plt.plot(r, analytic_C(N, r, T), label=ana_label)
    plt.scatter(num_r, corr, label=num_label)
plt.xlabel("r")
plt.ylabel("C(r)")
plt.xlim((-0.5, 27))
plt.legend()
# plt.title(r"Comparing analytical and numerical solution of $C(r)$")
plt.savefig("../tex/figs/plot_b.pdf")
plt.show()



