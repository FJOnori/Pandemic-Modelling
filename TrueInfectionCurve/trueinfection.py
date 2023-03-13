import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import random, randint

def deathCurve():
    COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)

    cutpoint =180

    dr = np.mean(np.array(COVIDdf['total_cases']) / np.array(COVIDdf['total_deaths']))
    start = list(dr*np.array(COVIDdf['new_deaths_smoothed']))[:cutpoint]
    end = list(np.array(COVIDdf['new_cases_smoothed']))[cutpoint:]
    infectioncurve = np.concatenate((start,end))

    kernel_size = 5
    kernel = np.ones(kernel_size) / kernel_size
    infectioncurve_convolved = np.convolve(infectioncurve, kernel, mode='same')

    plt.grid(ls=":",c='grey')
    plt.plot(infectioncurve_convolved, c='#b81111')
    plt.xlim(0,1000)
    plt.ylim(0,200000)
    plt.title("Estimated Infection Curve")
    plt.xlabel("Time (days)")
    plt.ylabel("Total Infections")
    plt.savefig("TrueInfectionCurveEstimate.png", dpi=227)
    plt.show()

deathCurve()