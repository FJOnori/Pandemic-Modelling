import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from random import random, randint

class SIRS_Model():

    def __init__(self) -> None:
        COVIDdf = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        self.T_final = len(np.array(COVIDdf['reproduction_rate']))
        self.S_initial = 67508935
        self.I_inital = 1
        self.population = self.S_initial + self.I_inital
        
        self.BirthRate = 2.78363e-5
        self.DeathRate = 2.48321e-5
        self.COVIDDeathRate = 8.95027e-3
        self.ImmunityLossRate = 1/365
        self.RecoveryRate = 1/4
        self.InfectionRate = (np.array(COVIDdf['reproduction_rate']))*((1/self.RecoveryRate)/2) - 1


    def NumInt(self):

        S ,I, R, T, D = [self.S_initial], [self.I_inital], [0], [0], [0]
        
        for n in np.arange(0, self.T_final):

               dS = (self.BirthRate*self.population) - self.InfectionRate[n]*(I[n]/self.population)*S[n] \
                + self.ImmunityLossRate*R[n] - self.DeathRate*S[n]

               dI = self.InfectionRate[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n] - self.COVIDDeathRate*I[n] 

               dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n]

               self.population += dI + dR + dS

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)
               D.append(self.COVIDDeathRate*I[n])

        return S, I, R, D, T

    def LinePlot(self):
        S ,I, R, D, T = self.NumInt()
        plt.plot(S, label="Susceptible")
        plt.plot(I, label="Infected")
        plt.plot(R, label="Recovered")
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


    def VectorFieldSR(self,infected_Slice = 400000,spacing = 25):
        
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

    def VectorFieldSI(self, Recovered_Slice = 400000, spacing = 25 ):
        
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

    def VectorFieldIR(self,Sus_Slice = 400000, spacing = 25):
  
        i, r = np.meshgrid(np.linspace(1,self.population,spacing),
                            np.linspace(1,self.population,spacing))

        dI = self.InfectionRate*(i/self.population)*Sus_Slice - self.RecoveryRate*i - self.COVIDDeathRate*i

        dR = self.RecoveryRate*i - self.ImmunityLossRate*r - self.DeathRate*r

        plt.quiver(r,i,dR,dI)
        plt.xlabel("Suseptable")
        plt.ylabel("Infected")
        plt.title("SIRS Model Vector field")
        plt.show()

    def grid_frame(self, sim, mtr):
        
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


    def setup_grid(gridsize=10, max_radius=1):
        row = ["*"]*max_radius + ["S"]*gridsize + ["*"]*max_radius
        edge = ["*"]*(gridsize+(max_radius*2))
        sim = np.array([edge]*max_radius + [row]*gridsize + [edge]*max_radius)
        sim[int(gridsize/2), int(gridsize/2)] = "I"
        return sim

    def grid_sim(self, simtime = 500, gridsize=10, max_radius=1):
        S = I = R = T = list()
        sim = self.setup_grid(gridsize,max_radius)

        for j in range(0,simtime):
            sim = self.grid_frame(sim,max_radius)
            S.append(np.sum(sim == "S"))
            I.append(np.sum(sim == "I"))
            R.append(np.sum(sim == "R"))

        return S,I,R

    def plot_grid_sim(simtime=500, gridsize=10, max_radius=1):
        S,I,R = self.grid_sim(self, simtime, gridsize, max_radius)
        plt.plot(S, label = "Suseptable")
        plt.plot(I, label = "Infectious")
        plt.plot(R, label = "Recovered")
        plt.legend()
        plt.show()
    
    def print_grid_sim(simtime=500, gridsize=10, max_radius=1):

        sim = self.setup_grid(gridsize,max_radius)
        for j in range(0,simtime):
            sim = self.grid_frame(sim,max_radius)
            print()
            print(j)
            print(sim)

    def animate_grid_sim():
        pass



if __name__ == "__main__":
    SIRS = SIRS_Model()
    SIRS.animation_grid()