import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import random, randint

class SIRS_Model():

    def __init__(self) -> None:
        self.T_final = 1000
        self.S_initial = 67508935
        self.I_inital = 1
        self.population = self.S_initial + self.I_inital
        
        self.BirthRate = 2.78363e-5
        self.DeathRate = 2.48321e-5
        self.COVIDDeathRate = 8.95027e-3
        self.InfectionRate = 0.5
        self.ImmunityLossRate = 1/365
        self.RecoveryRate = 1/4


    def NumInt(self):

        S ,I, R, T, D = [self.S_initial], [self.I_inital], [0], [0], [0]
        
        for n in np.arange(0, self.T_final):

               dS = (self.BirthRate*self.population) - self.InfectionRate*(I[n]/self.population)*S[n] \
                + self.ImmunityLossRate*R[n] - self.DeathRate*S[n]

               dI = self.InfectionRate*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n] - self.COVIDDeathRate*I[n] 

               dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n]

               dD = self.COVIDDeathRate*I[n] 

               self.population += dI + dR + dS

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)
               D.append(dD)
               T.append(n)

        return S, I, R, D, T

    
    def LinePlot(self):
        S ,I, R, D, T = self.NumInt()
        plt.plot(T, S, label="Susceptible")
        plt.plot(T, I, label="Infected")
        plt.plot(T, R, label="Recovered")
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.title("SIRS Pandemic Model Lineplot")
        plt.legend()
        plt.show()

    def StackPlot(self):
        S ,I, R, D, T = self.NumInt()
        plt.stackplot(T, S, I, R, labels=["Susceptible","Infected","Recovered"])
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.xlim(0,self.T_final)
        plt.ylim(0,self.population)
        plt.title("SIRS Pandemic model StackPlot")
        plt.legend()
        plt.show()

    def InfectionsPlot(self):
        S ,I, R, D, T = self.NumInt()
        plt.plot(T, I, label="Infections")
        plt.xlabel("Time (Days)")
        plt.ylabel("Infections from COVID-19")
        plt.title("SIRS Pandemic model Infections")
        plt.show()

    def DeathsPlot(self):
        S ,I, R, D, T = self.NumInt()
        plt.plot(T, D, label="Deaths", c='k')
        plt.xlabel("Time (Days)")
        plt.ylabel("Deaths from COVID-19")
        plt.title("SIRS Pandemic model Deaths")
        plt.show()

    def PieChart(self):
        S ,I, R, D, T = self.NumInt()
        data = np.array([S[-1], I[-1], R[-1]])
        l = ["Susceptible", "Infected", "Recovered"]
        plt.pie(data, labels=l)
        plt.title("SIRS Pandemic Final State")
        plt.show()


    def VectorFieldSR(self):
        
        spacing = 25
        infected_Slice = 400000
        s, r = np.meshgrid(np.linspace(1,self.population,spacing),
                           np.linspace(1,self.population,spacing))

        dS = (self.BirthRate*self.population) - self.InfectionRate*(infected_Slice/self.population)*s \
                + self.ImmunityLossRate*r- self.DeathRate*s

        dR = self.RecoveryRate*infected_Slice - self.ImmunityLossRate*r - self.DeathRate*r

        plt.quiver(s,r,dS,dR)
        plt.xlabel("Suseptable")
        plt.ylabel("Recovered")
        plt.title("SIRS Model Vector field")
        plt.show()

    def VectorFieldSI(self):
        
        spacing = 25
        Recovered_Slice = 400000
        s, i = np.meshgrid(np.linspace(1,self.population,spacing),
                            np.linspace(1,self.population,spacing))

        dS = (self.BirthRate*self.population) - self.InfectionRate*(i/self.population)*s \
            + self.ImmunityLossRate*Recovered_Slice - self.DeathRate*s

        dI = self.InfectionRate*(i/self.population)*s - self.RecoveryRate*i - self.COVIDDeathRate*i

        dR = self.RecoveryRate*i - self.ImmunityLossRate*Recovered_Slice - self.DeathRate*Recovered_Slice

        plt.quiver(s,i,dS,dI)
        plt.xlabel("Suseptable")
        plt.ylabel("Infected")
        plt.title("SIRS Model Vector field")
        plt.show()

    def VectorFieldIR(self):
        
        spacing = 25
        Sus_Slice = 400000
        i, r = np.meshgrid(np.linspace(1,self.population,spacing),
                            np.linspace(1,self.population,spacing))

        dI = self.InfectionRate*(i/self.population)*Sus_Slice - self.RecoveryRate*i - self.COVIDDeathRate*i

        dR = self.RecoveryRate*i - self.ImmunityLossRate*r - self.DeathRate*r

        plt.quiver(r,i,dR,dI)
        plt.xlabel("Suseptable")
        plt.ylabel("Infected")
        plt.title("SIRS Model Vector field")
        plt.show()

    def animation_grid_frame(self, sim, mtr):
        
        InfectionRate = 0.1
        RecoveryRate = 0.1
        ImmunitylossRate = 0.05

        newsim = sim.copy()
        for x in range(len(sim[0])-mtr):
            for y in range(len(sim[0])-mtr):                
                if sim[x,y] == "S":
                    if (sim[y,x+randint(1,mtr)] == "I" or sim[y,x-randint(1,mtr)] == "I" \
                    or sim[y+randint(1,mtr),x] == "I" or sim[y-randint(1,mtr),x] == "I" \
                    or sim[y-randint(1,mtr),x+randint(1,mtr)] == "I" or sim[y+randint(1,mtr),x-randint(1,mtr)] == "I"
                    or sim[y-randint(1,mtr),x-randint(1,mtr)] == "I" or sim[y+randint(1,mtr),x+randint(1,mtr)] == "I") \
                    and random() < InfectionRate:
                        newsim[x,y] = "I"
                elif sim[x,y] == "I":
                    if random() < RecoveryRate:
                        newsim[x,y] = "R"
                elif sim[x,y] == "R":
                    if random() < ImmunitylossRate:
                        newsim[x,y] = "S"

        return newsim

    def animation_grid(self):
        simtime = 500
        S = list()
        I = list()
        R = list()
        T = list()
        gridsize = 100
        mtr = 2
        row = ["*"]*mtr + ["S"]*gridsize + ["*"]*mtr
        edge = ["*"]*(gridsize+(mtr*2))
        sim = np.array([edge]*mtr + [row]*gridsize + [edge]*mtr)
        sim[int(gridsize/2), int(gridsize/2)] = "I"

        for j in range(0,simtime):
            print(j+1)
            sim = self.animation_grid_frame(sim,mtr)
            S.append(np.sum(sim == "S"))
            I.append(np.sum(sim == "I"))
            R.append(np.sum(sim == "R"))
            T.append(j)

        plt.plot(T,S)
        plt.plot(T,I)
        plt.plot(T,R)
        plt.show()


if __name__ == "__main__":
    SIRS = SIRS_Model()
    SIRS.animation_grid()