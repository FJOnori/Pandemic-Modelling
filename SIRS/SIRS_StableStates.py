import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def StableStates():
    S = []
    I = []
    R = []

    for b in np.linspace(0,1,100):
        S.append((0.1/b))
        I.append((1-(0.1/b))/(1+(210*0.1)))
        R.append((210*0.1)*((1-(0.1/b))/(1+(210*0.1))))

    plt.grid(ls=":",c='grey')
    plt.plot(np.linspace(0,1,100), S, c='#1e68b3', label="Susceptible")
    plt.plot(np.linspace(0,1,100), I, c='#b81111', label="Infected")
    plt.plot(np.linspace(0,1,100), R, c='#aaacad', label="Recovered")
    plt.ylabel("Proportion of Population")
    plt.xlabel("Contact Rate")
    plt.title("SIRS COVID-19 stable states")
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.savefig("SIRS/SIRSstablestates.png", dpi=227)


StableStates()