import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIR_Model():

    def __init__(self) -> None:
        self.S_initial = 1 - 0.000005
        self.I_inital = 0.000005
        self.population = self.S_initial + self.I_inital
        self.RecoveryRate = 0.05
        self.ContactRate = 0.15
        self.ImmunityLossRate = 0.0035

    def NumInt(self):

        S ,I, R, = [self.S_initial], [self.I_inital], [0]
        IP = [0]

        for n in np.arange(0,1000):

               dS = - self.ContactRate*(I[n]/self.population)*S[n] + self.ImmunityLossRate*R[n]
               dI = self.ContactRate*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n]
               dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n]

               self.population += (dI + dR + dS)*0.01

               S.append(S[n] + dS*0.01)
               I.append(I[n] + dI*0.01)
               R.append(R[n] + dR*0.01)
    

        return np.array(S), np.array(I), np.array(R)
    
    def Lineplot(self):
        S,I,R = self.NumInt()
        plt.plot(S, c='#1e68b3', label="Susceptible")
        plt.plot(I, c='#b81111', label="Infected")
        plt.plot(R, c='#aaacad', label="Recovered")
        plt.title("SIRS Numerical Integration")
        plt.ylabel("Population Proportion")
        plt.xlabel("Time (days)")
        plt.ylim(0,1)
        plt.legend()
        plt.savefig("SIRS/SIRSnumint.png", dpi=227)

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.Lineplot()