## Include any code for analyzing or visualizing the simulation outputs.

import matplotlib.pyplot as plt
import numpy as np
import math

rangeR = 40
recordInterval = 20
dT = .01

pi_array = np.genfromtxt("EKG_simulation_pi.csv", delimiter =",")
xi_array = np.genfromtxt("EKG_simulation_xi.csv", delimiter =",")
alpha_array = np.genfromtxt("EKG_simulation_alpha.csv", delimiter =",")
beta_array = np.genfromtxt("EKG_simulation_beta.csv", delimiter =",")
psi_array = np.genfromtxt("EKG_simulation_psi.csv", delimiter =",")
expansion_array = np.genfromtxt("EKG_simulation_expansion.csv", delimiter =",")


nT = pi_array.shape[0]
nR = pi_array.shape[1]



r_array = np.linspace(0,rangeR, nR)


for i in range(nT):
    plt.plot(r_array, xi_array[i,:])
    plt.title(i* recordInterval*dT)
    #plt.plot(r_array[3700:4000], pi_array[i,3700:4000])
    plt.draw()
    plt.pause(.1)
    plt.clf()
