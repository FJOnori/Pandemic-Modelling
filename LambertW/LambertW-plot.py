import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

def LambertW_plot():

    x = []
    y = []
    for r in np.arange(0,100,0.05):
        y.append(lambertw(r))
        x.append(r)

    plt.grid(ls=":",c='grey')
    plt.plot(np.arange(0,100,0.05),y, c="green")
    plt.title("Lambert W function")
    plt.xlabel("X-input")
    plt.ylabel("Y-output")
    plt.xlim(0,100)
    plt.ylim(0,3.5)
    plt.savefig("LambertW/lambertW.png", dpi=227)
    plt.show()


def stable_states():
    N = 1
    s0 = N - 0.000001
    R0 = []
    Sw = []

    for r in np.arange(1,3.5,0.05):
        Wsub = -(r/N)*s0*np.exp(-r)
        Wmul = -(N/r)
        Sw.append(Wmul*lambertw(Wsub))
        R0.append(r)
        
    plt.grid(ls=":",c='grey')
    plt.plot(R0, np.array(Sw), c="#1e68b3", label="Susceptible")
    plt.plot(R0, (N - np.array(Sw)), c="#aaacad", label="Recovered")
    plt.title("SIR Stable States")
    plt.xlabel("R\u2080 value")
    plt.ylabel("Proportion of Population")
    plt.xlim(1,3)
    plt.ylim(0,1)
    plt.legend()
    plt.savefig("LambertW/Stable-States.png",dpi=227)


LambertW_plot()