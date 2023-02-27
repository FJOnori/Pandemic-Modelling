import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIR_Model():

    def __init__(self) -> None:
        COVIDdf = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        self.S_initial = 1 - 0.0000056
        self.I_inital = 0.0000056
        self.population = self.S_initial + self.I_inital
        self.RecoveryRate = 1/10
        self.InfectionRate = 0.19
        self.TDInfectionRate = (self.RecoveryRate)*np.array(COVIDdf['reproduction_rate'][30:1030])

    def NumInt(self, TD=True):

        S ,I, R, = [self.S_initial], [self.I_inital], [0]

        if TD: beta = self.TDInfectionRate 
        else: beta = np.array([self.InfectionRate]*len(self.TDInfectionRate))

        for n in np.arange(0, len(self.TDInfectionRate)):

               dS = - beta[n]*(I[n]/self.population)*S[n] 
               dI = beta[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n]
               dR = self.RecoveryRate*I[n]

               self.population += dI + dR + dS

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)

        return np.array(S), np.array(I), np.array(R)
    
    def InfectionCurve(self):
        S,I,R = self.NumInt()
        plt.plot(I*67000000)
        plt.show()

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.InfectionCurve()