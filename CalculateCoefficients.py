import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint


COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)


def ContactRateH():
    betah = (1/10) * (np.array(COVIDdf['reproduction_rate']))[30:53]
    betal = (1/10) * (np.array(COVIDdf['reproduction_rate']))[53:113]
    print(np.mean(betah))
    print(np.mean(betal))


def COVIDdeathrate():
    d = np.array(COVIDdf['total_deaths'])
    c = np.array(COVIDdf['total_cases'])
    dr = d[-1]/c[-1]
    print(dr**-1)
    return dr

def deathrate():
    pass

def birthrate():
    pass

ContactRateH()