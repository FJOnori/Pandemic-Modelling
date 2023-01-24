import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class SEIRSV2_Model():

    def __init__(self):

        self.S_initial = 67570000
        self.E_inital = 1
        self.population = self.S_initial + self.E_inital

        self.T_final = 1100
        self.dt = 0.04

        self.beta = 0.211 #exposure rate, S to E
        self.gamma = 1/12 #recovery rate, I to R
        self.epsilon = 1/365 #immunity loss rate, R to S
        self.sigma = 1/4 #infection rate, E to I
        self.upsilon = 0.001 #vaccination rate S,R to V
        self.mu = 1/28092 #Birth rate
        self.nu = 1/29565 #Death rate

    def NumInt(self):

        S ,I, E ,R, V1, V2, T, n = [self.S_initial], [0], [self.E_inital], [0], [0], [0], [0], 0

        for t in np.arange(0, self.T_final, self.dt):
            
            dS = (self.mu*self.population - self.beta*((I[n]+ E[n])/self.population)*S[n] + self.epsilon*R[n] - self.nu*S[n] - self.upsilon*S[n]) * self.dt
            dE = (self.beta*((I[n] + E[n])/self.population)*S[n] - self.sigma*E[n] - self.nu*E[n]) * self.dt
            dI = (self.sigma*E[n] - self.gamma*I[n] - self.nu*I[n]) * self.dt
            dR = (self.gamma*I[n] - self.epsilon*R[n] - self.nu*R[n] - self.upsilon*R[n]) * self.dt
            dV1 = (self.upsilon*S[n] + self.upsilon*R[n] - self.upsilon*V1[n] - self.nu*V1[n]) * self.dt
            dV2 = (self.upsilon*V1[n] - self.nu*V2[n]) * self.dt
            
            S.append(S[n] + dS)
            E.append(E[n] + dE)
            I.append(I[n] + dI)
            R.append(R[n] + dR)
            V1.append(V1[n] + dV1)
            V2.append(V2[n] + dV2)
            T.append(T[n] + self.dt)

            n += 1

        plt.plot(T, S, label="Susceptible")
        plt.plot(T, E, label="Exposed")
        plt.plot(T, I, label="Infected")
        plt.plot(T, R, label="Recovered")
        plt.plot(T, V1, label="1 Vaccinated")
        plt.plot(T, V2, label="2 Vaccinated")
        plt.xlabel("Time (Days)")
        plt.ylabel("Population (Percentage)")
        plt.title("SEVIRS Pandemic model")
        plt.legend()
        plt.show()

    def calculate_vaccination_rate(self):
        df = pd.read_csv('owid-covid-data-uk.csv')
        v1_list = list(df['people_vaccinated_per_hundred'].fillna(0))
        v2_list = list(df['people_fully_vaccinated_per_hundred'].fillna(0))
        plt.plot(range(0,len(v1_list)), v1_list)
        plt.plot(range(0,len(v2_list)), v2_list)
        plt.ylim(0,100)
        plt.show()

    def calculate_beta(self):
        currentCases = 0
        newCases = 0
        suseptablePopulation = self.population - currentCases
        beta = - (newCases + self.nu*suseptablePopulation + self.upsilon*suseptablePopulation - self.mu*self.population)/(((currentCases)/self.population)*suseptablePopulation)


SEIRSV2 = SEIRSV2_Model()
SEIRSV2.calculate_vaccination_rate()