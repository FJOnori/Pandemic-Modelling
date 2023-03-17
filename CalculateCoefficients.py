import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint


COVIDdf = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)

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

COVIDdeathrate()