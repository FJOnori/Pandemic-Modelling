import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import random

class SEIRSV2_Model():

    def __init__(self):
    
        self.T_final = len(np.array(pd.read_csv("CSVfiles/beta.csv")['beta']))
        self.S_initial = 67508935
        self.E_inital = 1
        self.population = self.S_initial + self.E_inital
        self.V1rate = np.average(np.trim_zeros(np.array(pd.read_csv("CSVfiles/Vaccine1rate.csv")['vaccination_rate'])))
        self.V2rate = np.average(np.trim_zeros(np.array(pd.read_csv("CSVfiles/Vaccine2rate.csv")['vaccination_rate'])))
        self.vaccineDay = 600
        
        self.ExposureRate = np.array(pd.read_csv("CSVfiles/beta.csv")['beta'])
        self.ExposureRateV1 = self.ExposureRate * 0.19
        self.ExposureRateV2 = self.ExposureRate * 0.09
        self.RecoveryRate = 1/10
        self.RecoveryRateV1 = 1/10
        self.RecoveryRateV2 = 1/10
        self.ImmunityLossRate = 1/365
        self.InfectionRate = 0.17857
        self.BirthRate = 2.78363e-5
        self.DeathRate = 2.48321e-5
        self.COVIDDeathRate = 8.95027e-3
        self.COVIDDeathRate1 = self.COVIDDeathRate * 0.19
        self.COVIDDeathRate2 = self.COVIDDeathRate * 0.09
        self.VaccinationRate1 = np.array([0]*(self.vaccineDay) + [self.V1rate]*(self.T_final - self.vaccineDay))
        self.VaccinationRate2 = np.array([0]*(self.vaccineDay) + [self.V2rate]*(self.T_final - self.vaccineDay))

    def NumInt(self):

        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = [self.S_initial], [0], [self.E_inital], [0], [0], [0], [0], [0], [0], [0], [0], [0]
        totalcases = 0

        for n in np.arange(0, self.T_final):
            
            infec = I[n] + E[n] + I1[n] + E1[n] + I2[n] + E2[n]

            dS = (self.BirthRate*self.population) - self.ExposureRate[n]*((infec)/self.population)*S[n] \
                + self.ImmunityLossRate*R[n] - self.DeathRate*S[n] - self.VaccinationRate1[n]*S[n]

            dE = self.ExposureRate[n]*((infec)/self.population)*S[n] - self.InfectionRate*E[n] - self.DeathRate*E[n] 

            dI = self.InfectionRate*E[n] - self.RecoveryRate*I[n] - self.COVIDDeathRate*I[n]

            dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n] \
                 - self.VaccinationRate1[n]*R[n]

            dE1 = self.ExposureRateV1[n]*((infec)/self.population)*V1[n] - self.InfectionRate*E1[n] - self.DeathRate*E1[n]

            dI1 = - self.RecoveryRateV1*I1[n] + self.InfectionRate*E1[n] - self.COVIDDeathRate1*I1[n]

            dE2 = self.ExposureRateV2[n]*((infec)/self.population)*V2[n] - self.InfectionRate*E2[n] - self.DeathRate*E2[n]

            dI2 =  - self.RecoveryRateV2*I2[n] + self.InfectionRate*E2[n] - self.COVIDDeathRate2*I2[n]

            dV1 = self.VaccinationRate1[n]*S[n] + self.VaccinationRate1[n]*R[n] - self.VaccinationRate2[n]*V1[n] \
                 - self.ExposureRateV1[n]*((infec)/self.population)*V1[n] - self.DeathRate*V1[n] + self.RecoveryRateV1*I1[n]

            dV2 = self.VaccinationRate2[n]*V1[n] - self.DeathRate*V2[n] - self.ExposureRateV2[n]*((infec)/self.population)*V2[n] \
                + self.RecoveryRateV2*I2[n]

            COVIDDeaths = self.COVIDDeathRate*I[n] + self.COVIDDeathRate1*I1[n] + self.COVIDDeathRate2*I2[n]
            totalcases +=  self.ExposureRateV2[n]*((infec)/self.population)*V2[n] + self.ExposureRateV1[n]*((infec)/self.population)*V1[n] + self.ExposureRate[n]*((infec)/self.population)*S[n]
            self.population += dS + dE + dE1 + dE2 + dI + dI1 + dI2 + dR + dV1 + dV2

            S.append(S[n] + dS)
            E.append(E[n] + dE)
            E1.append(E1[n] + dE1)
            E2.append(E2[n] + dE2)
            I.append(I[n] + dI)
            I1.append(I1[n] + dI1)
            I2.append(I2[n] + dI2)
            R.append(R[n] + dR)
            V1.append(V1[n] + dV1)
            V2.append(V2[n] + dV2)
            D.append(COVIDDeaths)     
            T.append(n)

        print("total cases: ", int(totalcases))
        print("total deaths: ", int(sum(D)))
        print("vaccinated: ", int((np.array(I1) + np.array(E1) + np.array(V1))[-1]))
        print("fully vaccinated: ", int((np.array(I2) + np.array(E2) + np.array(V2))[-1]))
        return  S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D 

    def LinePlot(self):
        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = self.NumInt()
        plt.plot(T, S, label="Susceptible")
        plt.plot(T, np.array(I)+np.array(I1)+np.array(I2)+np.array(E)+np.array(E1)+np.array(E2), label="Infected")
        plt.plot(T, R, label="Recovered")
        plt.plot(T, V2, label="Fully Vaccinated")
        plt.plot(T, V1, label="Vaccinated")
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.title("SEIRSV2 Pandemic model Lineplot")
        plt.legend()
        plt.show()

    def StackPlot(self):
        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = self.NumInt()
        plt.stackplot(T, S, E, I, R, V1, V2,  labels=["Susceptible","Exposed","Infected","Recovered","1st Dose", "2nd Dose"])
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.xlim(0,self.T_final)
        plt.ylim(0,self.population)
        plt.title("SEIRSV2 Pandemic model StackPlot")
        plt.legend()
        plt.show()

    def InfectionsPlot(self):
        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = self.NumInt()
        plt.plot(T, np.array(E) + np.array(I), label="Infections")
        plt.xlabel("Time (Days)")
        plt.ylabel("Infections from COVID-19")
        plt.title("SEIRSV2 Pandemic model Deaths")
        plt.show()

    def DeathsPlot(self):
        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = self.NumInt()
        plt.plot(T, D, label="Deaths", c='k')
        plt.xlabel("Time (Days)")
        plt.ylabel("Deaths from COVID-19")
        plt.title("SEIRSV2 Pandemic model Deaths")
        plt.show()

    def PieChart(self):
        S ,I, E ,R, V1, V2, T, E2, E1, I1, I2, D = self.NumInt()
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
        simtime = 10
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
        
        colourMap = {"S":"blue",
                      "E":"orange",
                      "E1":"orange",
                      "E2":"orange",
                      "I":"red",
                      "I1":"red",
                      "I2":"red",
                      "R": "green",
                      "1": "grey",
                      "2":"dimgrey"}

        for j in range(0,simtime):
            print()
            print(str(j))
            print(sim)
            sim = self.animation_grid_frame(sim)

        

if __name__ == "__main__":
    SEIRSV2 = SEIRSV2_Model()
    SEIRSV2.LinePlot()