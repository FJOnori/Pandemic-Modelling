import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import random, randint
from scipy.ndimage.interpolation import shift


def deathCurve():
    
    COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
    dr = np.mean( shift(np.array(COVIDdf['new_cases_smoothed'])[109:313], -10,cval=0)  / np.array(COVIDdf['new_deaths_smoothed'])[109:313] )
    infec = shift(dr*np.array(COVIDdf['new_deaths_smoothed']), -10, cval=0)
    recinfec = np.array(COVIDdf['new_cases_smoothed'])
 
    plt.grid(ls=":",c='grey',axis='y')
    plt.plot(infec, c='#b81111')
    plt.plot(recinfec, c='#b81111', ls="--")
    plt.xlim(0,200)
    plt.ylim(0,120000)
    plt.title("Estimated Infection Curve")
    plt.xlabel("Time (days)")
    plt.ylabel("Total Infections")
    plt.savefig("TrueInfectionsEstimate.png", dpi=227)
    plt.show()

deathCurve()