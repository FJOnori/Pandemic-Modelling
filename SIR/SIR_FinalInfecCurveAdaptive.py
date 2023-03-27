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


    def __init__(self, FinalTime=1000, dt=0.01, I_initial=0.0001, RecoveryRate=(1/10), ContactFactor=1) -> None:

        COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)

        self.co

        self.dt             = dt
        self.FinalTime      = FinalTime
        self.iterations     = int(self.FinalTime/self.dt)
        self.ContactFactor  = ContactFactor
        self.population     = 1
        self.I_inital       = I_initial
        self.S_initial      = self.population - I_initial
        self.RecoveryRate   = np.round(RecoveryRate,3)
        self.popcut         = 0.01
        self.ContactRateH   = self.RecoveryRate*np.array(COVIDdf['reproduction_rate'])[:52]
        self.ContactRateL   = self.RecoveryRate*(ContactFactor)*np.array(COVIDdf['reproduction_rate'])[52:1000]

    def NumInt(self):

        IA = []
        for b in np.linspace(self.ContactRateL,self.ContactRateH,5):
            
            S ,I, R, T = [self.S_initial], [self.I_inital], [0], [0]

            beta = np.concatenate(self.ContactRateH, self.ContactRateL)

            for n in np.arange(0, self.iterations):                
                dS = (-beta[n]*(I[n]/self.population)*S[n])*self.dt
                dI = (beta[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n])*self.dt
                dR = (self.RecoveryRate*I[n])*self.dt

                S.append(S[n] + dS)
                I.append(I[n] + dI)
                R.append(R[n] + dR)
                T.append(self.dt*n)
                
                self.progress_bar(n, self.iterations)
            
            print()
            IA.append(np.array(I)*67000000)

        return np.array(IA), np.array(T)
    
    def Lineplot(self):
        IA, T = self.NumInt()
        labels = np.round(np.linspace(self.ContactRateL,self.ContactRateH,5),4)
        for j in range(0, len(IA)):
            plt.plot(T, IA[j], label=str(np.round(labels[j],2)), ls="--", alpha=0.8)

        plt.title("SIR Adaptive Model, Population Cutoff = " + str(self.popcut))
        plt.ylabel("Infected Population")
        plt.xlabel("Time (days)")
        plt.ylim(0)
        plt.xlim(0,self.FinalTime)
        plt.grid(ls=":", c="grey", axis='y')
        plt.legend()
        plt.savefig("SIRFinalInfectionCurveLockdown - " + str(self.popcut) + ".png", dpi=227)

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.Lineplot()