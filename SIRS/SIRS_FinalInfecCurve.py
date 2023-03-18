import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIR_Model():

    def __init__(self) -> None:
        
        self.FinalTime = 1000
        self.dt = 0.01
        self.iterations = int(self.FinalTime/self.dt)
        COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        infectioncurve = np.array(list(np.array(COVIDdf['new_cases_smoothed']))[:1000])/67000000
        kernel_size = 5
        kernel = np.ones(kernel_size) / kernel_size
        self.TrueInfectionCurve = np.convolve(infectioncurve, kernel, mode='same')

        self.S_initial = 1 - 0.00002
        self.I_inital = 0.00002
        self.population = self.S_initial + self.I_inital
        self.RecoveryRate = 1/7
        self.ContactRate = np.array(COVIDdf['reproduction_rate'])[:1000]*self.RecoveryRate
        self.ImmunityLossRate = 1/365
    
    def progress_bar(self, current, total, bar_length=100):
        fraction = current / total
        arrow = int(fraction * bar_length - 1) * '-' + '>'
        padding = int(bar_length - len(arrow)) * ' '
        ending = '\n' if current == total else '\r'
        print(f'Progress: [{arrow}{padding}] {round((fraction*100),5)}%', end=ending)

    def NumInt(self):

        S ,I, R, T = [self.S_initial], [self.I_inital], [0], [0]
        beta = np.repeat(self.ContactRate, self.dt**(-1))

        for n in np.arange(0,self.iterations):

               dS = (-beta[n]*(I[n]/self.population)*S[n] + self.ImmunityLossRate*R[n])*self.dt
               dI = (beta[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n])*self.dt
               dR = (self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n])*self.dt

               self.progress_bar(n, self.iterations)

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)
               T.append(self.dt*n)

        return np.array(S), np.array(I), np.array(R), np.array(T)
    
    def Lineplot(self):
        S,I,R,T = self.NumInt()
        plt.plot(T,I*67000000, c='#b81111', label="Infected")
        plt.title("SIRS Numerical Integration")
        plt.ylabel("Infected Population")
        plt.xlabel("Time (days)")
        plt.xlim(0,self.FinalTime)
        plt.ylim(0,180000)
        plt.savefig("SIRS/SIRSFinalInfectionCurve.png", dpi=227)

if __name__ == "__main__":
    SIR = SIR_Model()
    SIR.Lineplot()