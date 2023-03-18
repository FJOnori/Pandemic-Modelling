import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIR_Model():

    def progress_bar(self, current, total, bar_length=100):
        fraction = current / total
        arrow = int(fraction * bar_length - 1) * '-' + '>'
        padding = int(bar_length - len(arrow)) * ' '
        ending = '\n' if current == total else '\r'
        print(f'Progress: [{arrow}{padding}] {round((fraction*100),5)}%', end=ending)


    def __init__(self, FinalTime=1000, dt=0.01, I_initial=0.00005, RecoveryRate=(1/10), GraphRounding=4, ContactFactor=0.98) -> None:

        COVIDdf             = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        self.dt             = dt
        self.FinalTime      = FinalTime
        self.iterations     = int(self.FinalTime/self.dt)

        self.population     = 1
        self.I_inital       = I_initial
        self.S_initial      = self.population - I_initial
        
        self.RecoveryRate   = np.round(RecoveryRate,3)
        self.ContactRate    = np.array(COVIDdf['reproduction_rate'])[:1000]*self.RecoveryRate*ContactFactor
        self.GraphRounding  = GraphRounding

    def NumInt(self):

        S ,I, R, NI, T = [self.S_initial], [self.I_inital], [0], [0],[0]
        beta = np.repeat(self.ContactRate, self.dt**(-1))
        NIRT = 0

        for n in np.arange(0, self.iterations):

               dS = (-beta[n]*(I[n]/self.population)*S[n])*self.dt
               dI = (beta[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n])*self.dt
               dR = (self.RecoveryRate*I[n])*self.dt

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)
               T.append(self.dt*n)
               
               self.progress_bar(n, self.iterations)
            
               NIRT += (beta[n]*(I[n]/self.population)*S[n])*self.dt
               if n%(self.dt**(-1)) == 0:
                    NI.append(NIRT); NIRT = 0

        return np.array(S), np.array(I), np.array(R), np.array(NI), np.array(T)
    
    def Lineplot(self):
        S,I,R,NI,T = self.NumInt()
        plt.plot(T,I*67000000, c='#b81111', label="Infected")
        plt.title("SIR Model Infections $\gamma = "+ str(self.RecoveryRate) + "$")
        plt.ylabel("Infected Population")
        plt.xlabel("Time (days)")
        plt.xlim(0,self.FinalTime)
        plt.grid(ls=":", c="grey", axis='y')
        plt.savefig("SIR/SIRFinalInfectionCurve-" + str(self.RecoveryRate) + ".png", dpi=227)

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.Lineplot()