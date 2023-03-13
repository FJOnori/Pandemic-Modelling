import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIR_Model():

    def __init__(self) -> None:

        COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        infectioncurve = np.array(list(np.array(COVIDdf['new_cases_smoothed']))[:1000])/67000000
        kernel_size = 5
        kernel = np.ones(kernel_size) / kernel_size
        self.TrueInfectionCurve = np.convolve(infectioncurve, kernel, mode='same')

        self.S_initial = 1 - 0.0001
        self.I_inital = 0.0001
        self.population = self.S_initial + self.I_inital
        self.RecoveryRate = 0.1
        self.ContactRate = np.array(COVIDdf['reproduction_rate'])[:1000]*self.RecoveryRate
        self.ImmunityLossRate = 1/365

    def NumInt(self):

        S ,I, R, = [self.S_initial], [self.I_inital], [0]

        for n in np.arange(0,1000):

               dS = - self.ContactRate[n]*(I[n]/self.population)*S[n] + self.ImmunityLossRate*R[n]
               dI = self.ContactRate[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n]
               dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n]

               self.population += dI + dR + dS

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)

        return np.array(S), np.array(I), np.array(R)
    
    def Lineplot(self):
        S,I,R = self.NumInt()
        plt.plot(self.TrueInfectionCurve, c='#bfacac', label="Infected", ls=":")
        plt.plot(I, c='#b81111', label="Infected")
        plt.title("SIRS Numerical Integration")
        plt.ylabel("Population Proportion")
        plt.xlabel("Time (days)")
        plt.xlim(0,1000)
        plt.ylim(0,0.003)
        plt.savefig("SIRS/SIRSFinalInfectionCurve.png", dpi=227)

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.Lineplot()