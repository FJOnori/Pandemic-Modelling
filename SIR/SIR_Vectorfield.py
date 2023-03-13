import numpy as np
import matplotlib.pyplot as plt

def SIR_VectorField(S_initial=1-0.0001, I_inital=0.0001, beta=0.3, gamma=0.4):
    
    population = S_initial + I_inital
    x,y = np.meshgrid(np.linspace(0.001,1,100),np.linspace(0.001,1,100))
    dS = (-beta*(y/population)*x)
    dI = (beta*(y/population)*x - gamma*y)

    color = 2 * np.log(np.hypot(dS, dI))
    plt.streamplot(x,y,dS,dI, cmap="Greys", color=color, density=1.5)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIR Vector Field")
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.savefig("SIR/SIRvectorfield.png", dpi=227)

SIR_VectorField()