import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import sleep
from random import random

class SEIRSV2_Model():

    def __init__(self):
        
        self.T_final = 2000
        self.S_initial = 65000000
        self.E_inital = 1
        self.population = self.S_initial + self.E_inital

        self.ExposureRate = 0.211           #S to E
        self.ExposureRateV1 = 0.1             #V1 to E
        self.ExposureRateV2 = 0.05             #V2 to E
        self.RecoveryRate = 1/12            #I to R
        self.RecoveryRateV1 = 0.05             #I to V1 
        self.RecoveryRateV2 = 0.05             #I to V2
        self.ImmunityLossRate = 1/365       #R to S
        self.InfectionRate = 1/4            #E to I
        self.BirthRate = 0.001              #units
        self.DeathRate = 0.0001             #units
        self.COVIDDeathRate = 0.001             #units
        self.VaccinationRate1 = 0.0001           #units
        self.VaccinationRate2 = 0.0001           #units

    def NumInt(self):

        S ,I, E ,R, V1, V2, T, D = [self.S_initial], [0], [self.E_inital], [0], [0], [0], [0], [0]

        for n in np.arange(0, self.T_final):
            
            dS = (self.BirthRate*self.population) - self.ExposureRate*((I[n] + E[n])/self.population)*S[n] \
                + self.ImmunityLossRate*R[n] - self.DeathRate*S[n] - self.VaccinationRate1*S[n]

            dE = self.ExposureRate*((I[n] + E[n])/self.population)*S[n] - self.InfectionRate*E[n] - self.DeathRate*E[n] \
                + self.ExposureRateV1*((I[n] + E[n])/self.population)*V1[n] + self.ExposureRateV2*((I[n] + E[n])/self.population)*V2[n]

            dI = self.InfectionRate*E[n] - self.RecoveryRate*I[n] - self.COVIDDeathRate*I[n] \
                - self.RecoveryRateV1*I[n] - self.RecoveryRateV2*I[n]

            dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n] \
                 - self.VaccinationRate1*R[n]

            dV1 = self.VaccinationRate1*S[n] + self.VaccinationRate1*R[n] - self.VaccinationRate2*V1[n] \
                 - self.ExposureRateV1*((I[n] + E[n])/self.population)*V1[n] - self.DeathRate*V1[n] + self.RecoveryRateV1*I[n]

            dV2 = self.VaccinationRate2*V1[n] - self.DeathRate*V2[n] - self.ExposureRateV2*((I[n] + E[n])/self.population)*V2[n] \
                + self.RecoveryRateV2*I[n]

            COVIDDeaths = self.COVIDDeathRate*I[n]
            self.population += dS + dE + dI + dR + dV1 + dV2

            S.append(S[n] + dS)
            E.append(E[n] + dE)
            I.append(I[n] + dI)
            R.append(R[n] + dR)
            V1.append(V1[n] + dV1)
            V2.append(V2[n] + dV2)

            D.append(COVIDDeaths)     
            T.append(n)

        print("Final Population: ", int(self.population))
        print("Percentage Vaccinated:",str(100*V2[-1]/self.population) + "%")
        print("Total Deaths: ", int(D[-1]))
        return S, E, I , R, V1, V2, D, T

    def StackPlot(self):
        S, E, I, R, V1, V2, D, T = self.NumInt()
        plt.stackplot(T, S, E, I, R, V1, V2,  labels=["Susceptible","Exposed","Infected","Recovered","1st Dose", "2nd Dose"])
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.xlim(0,self.T_final)
        plt.ylim(0,self.population)
        plt.title("SEIRSV2 Pandemic model StackPlot")
        plt.legend()
        plt.show()

    def LinePlot(self):
        S, E, I, R, V1, V2, D, T = self.NumInt()
        plt.plot(T, S, label="Susceptible")
        plt.plot(T, E, label="Exposed")
        plt.plot(T, I, label="Infected")
        plt.plot(T, R, label="Recovered")
        plt.plot(T, V1, label="1st Dose")
        plt.plot(T, V2, label="2nd Dose")
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.title("SEIRSV2 Pandemic model Lineplot")
        plt.legend()
        plt.show()

    def InfectionsPlot(self):
        S, E, I, R, V1, V2, D, T = self.NumInt()
        plt.plot(T, np.array(E) + np.array(I), label="Infections")
        plt.xlabel("Time (Days)")
        plt.ylabel("Infections from COVID-19")
        plt.title("SEIRSV2 Pandemic model Deaths")
        plt.show()

    def DeathsPlot(self):
        S, E, I, R, V1, V2, D, T = self.NumInt()
        plt.plot(T, D, label="Deaths", c='k')
        plt.xlabel("Time (Days)")
        plt.ylabel("Deaths from COVID-19")
        plt.title("SEIRSV2 Pandemic model Deaths")
        plt.show()

    def PieChart(self):
        S, E, I, R, V1, V2, D, T = self.NumInt()
        data = np.array([S[-1], E[-1], I[-1], R[-1], V1[-1], V2[-1]])
        l = ["Susceptible","Exposed", "Infected", "Recovered", "1st Dose", "2nd Dose"]
        plt.pie(data, labels=l)
        plt.title("SEIRSV2 Pandemic model Final State")
        plt.show()

    def VectorField(self):
        s,i,e,r = np.meshgrid(np.linspace(1,self.population,25),\
                              np.linspace(1,self.population,25),\
                              np.linspace(1,self.population,25),\
                              np.linspace(1,self.population,25))

        dS = (self.BirthRate*self.population) - self.ExposureRate*((i + e)/self.population)*s \
                + self.ImmunityLossRate*r - self.DeathRate*s - self.VaccinationRate1*s

        dI = self.InfectionRate*e - self.RecoveryRate*i - self.COVIDDeathRate*i \
                - self.RecoveryRateV1*i - self.RecoveryRateV2*i
        
        plt.quiver(s,i,dS,dI)
        plt.xlabel("Susceptible")
        plt.ylabel("Infected")
        plt.title("SEIRSV2 vector field")
        plt.show()


    def ASCII_line(self):
        simtime = 10
        sim = np.array(["*","S","S","S","S","S","E","S","S","S","S","*"])
        
        for j in range(0,simtime):
            print(str(j), sim)

            newsim = sim.copy()
            for n in range(len(sim)-1):
                if sim[n] == "S":
                    if sim[n+1] == "E" or sim[n-1] == "E":
                        newsim[n] = "E"

            sim = newsim
            
    

    def animation_grid_frame(self, sim):
        
        newsim = sim.copy()
        for x in range(len(sim[0])-1):
            for y in range(len(sim[0])-1):

                if sim[x,y] == "S":
                    
                    if (sim[y,x+1] == "E" or sim[y,x-1] == "E" \
                    or sim[y+1,x] == "E" or sim[y-1,x] == "E" \
                    or sim[y,x+1] == "I" or sim[y,x-1] == "I" \
                    or sim[y+1,x] == "I" or sim[y-1,x] == "I") \
                    and random() < 0.9:
                        newsim[x,y] = "E"
                    
                    elif random() < 0.3:
                        sim[x,y] == "1"
                    
                elif sim[x,y] == "E":
                    
                    if random() < 0.5:
                        newsim[x,y] = "I"

                elif sim[x,y] == "E1":
                    
                    if random() < 0.5:
                        newsim[x,y] = "I1"
                
                elif sim[x,y] == "E2":
                    
                    if random() < 0.5:
                        newsim[x,y] = "I2"

                elif sim[x,y] == "I":
                    
                    if random() < 0.2:
                        newsim[x,y] = "R"

                elif sim[x,y] == "I1":
                    
                    if random() < 0.2:
                        newsim[x,y] = "1"
                
                elif sim[x,y] == "I2":
                    
                    if random() < 0.2:
                        newsim[x,y] = "2"

                elif sim[x,y] == "R":

                    if random() < 0.05:
                        newsim[x,y] = "S"
                    elif random() < 0.3:
                        newsim[x,y] = "1"
                
                elif sim[x,y] == "1":

                    if random() < 0.6:
                        newsim[x,y] = "2"
                    elif random() < 0.1:
                        newsim[x,y] = "E1"
                
                elif sim[x,y] == "2":

                    if random() < 0.1:
                        newsim[x,y] = "E2"


        return newsim
    
    def animation_grid(self):
        simtime = 100
        sim = np.array([["*","*","*","*","*","*","*","*","*","*","*","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","E","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","S","S","S","S","S","S","S","S","S","S","*"],
                        ["*","*","*","*","*","*","*","*","*","*","*","*"],])
        
        for j in range(0,simtime):
            print()
            print(str(j))
            print(sim)
            sim = self.animation_grid_frame(sim)

        

if __name__ == "__main__":
    SEIRSV2 = SEIRSV2_Model()
    SEIRSV2.animation_grid()