import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def read_uk_csv():
    path = "/Users/finnjohnonori/Documents/GitHubR/owid-covid-data.csv"
    df = pd.read_csv(path)
    df = df[df['location'] == 'United Kingdom']
    df.to_csv("/Users/finnjohnonori/Documents/GitHubR/owid-covid-data-uk.csv")



def print_df():
    df = pd.read_csv("owid-covid-data-uk.csv")
    print(df)

print_df()