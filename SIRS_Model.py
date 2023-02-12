import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from random import random, randint

class SIRS_Model():

    def __init__(self) -> None:
        COVIDdf = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        self.T_final = len(np.array(COVIDdf['reproduction_rate']))
        self.S_initial = 1 - 0.00001
        self.I_inital = 0.00001
        self.population = self.S_initial + self.I_inital
       
        self.BirthRate = 2.78363e-5
        self.DeathRate = 2.48321e-5
        self.COVIDDeathRate = 8.95027e-3
        self.ImmunityLossRate = 1/210
        self.RecoveryRate = 1/4
        self.InfectionRate = np.array(COVIDdf['reproduction_rate'])*1/4


    def NumInt(self):

        S ,I, R, D = [self.S_initial], [self.I_inital], [0], [0]
        
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

        return S, I, R, D

    def LinePlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(S, label="Susceptible")
        plt.plot(I, label="Infected")
        plt.plot(R, label="Recovered")
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.title("SIRS Pandemic Model Lineplot")
        plt.legend()
        plt.show()

    def StackPlot(self):
        S ,I, R, D = self.NumInt()
        plt.stackplot(T, S, I, R, labels=["Susceptible","Infected","Recovered"])
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.xlim(0,self.T_final)
        plt.ylim(0,self.population)
        plt.title("SIRS Pandemic model StackPlot")
        plt.legend()
        plt.show()

    def InfectionsPlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(I, label="Infections")
        plt.xlabel("Time (Days)")
        plt.ylabel("Infections from COVID-19")
        plt.title("SIRS Pandemic model Infections")
        plt.show()

    def DeathsPlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(D, label="Deaths", c='k')
        plt.xlabel("Time (Days)")
        plt.ylabel("Deaths from COVID-19")
        plt.title("SIRS Pandemic model Deaths")
        plt.show()

    def PieChart(self):
        S ,I, R, D = self.NumInt()
        data = np.array([S[-1], I[-1], R[-1]])
        l = ["Susceptible", "Infected", "Recovered"]
        plt.pie(data, labels=l)
        plt.title("SIRS Pandemic Final State")
        plt.show()

    def VectorField(self, Slice = 0.1, axis="IS", time= 800, spacing = 25):
        
        x, y = np.meshgrid(np.linspace(0,self.population,spacing),
                           np.linspace(0,self.population,spacing))

        if axis == "SR":
            dS = (self.BirthRate*self.population) - self.InfectionRate[time]*(Slice/self.population)*x + self.ImmunityLossRate*y- self.DeathRate*x
            dR = self.RecoveryRate*Slice - self.ImmunityLossRate*y - self.DeathRate*y
            plt.quiver(x,y,dS,dR)
            plt.xlabel("Susceptible")
            plt.ylabel("Recovered")
        elif axis == "IR":
            dI = self.InfectionRate[time]*(x/self.population)*Slice - self.RecoveryRate*x - self.COVIDDeathRate*x 
            dR = self.RecoveryRate*x - self.ImmunityLossRate*y - self.DeathRate*y
            plt.quiver(x,y,dI,dR)
            plt.xlabel("Infected")
            plt.ylabel("Recovered")
        elif axis == "IS":
            dS = (self.BirthRate*self.population) - self.InfectionRate[time]*(x/self.population)*y + self.ImmunityLossRate*Slice - self.DeathRate*y
            dI = self.InfectionRate[time]*(x/self.population)*y - self.RecoveryRate*x - self.COVIDDeathRate*x
            plt.quiver(x,y,dI,dS)
            plt.xlabel("Infected")
            plt.ylabel("Susceptible")
            
        plt.title("SIRS Model Vector field")
        plt.show()

    def setup_grid(self, gridsize=10, max_radius=1):
        row = [7]*max_radius + [1]*gridsize + [7]*max_radius
        edge = [7]*(gridsize+(max_radius*2))
        sim = np.array([edge]*max_radius + [row]*gridsize + [edge]*max_radius)
        sim[int(gridsize/2), int(gridsize/2)] = 3
        return sim

    def grid_frame(self, sim, max_radius):
        
        InfectionRate = 0.6
        RecoveryRate = 0.1
        ImmunitylossRate = 0.05
        tr = randint(1,max_radius)
        dtr = int(tr * (1/np.sqrt(2)))

        newsim = sim.copy()
        for x in range(len(sim[0])-max_radius):
            for y in range(len(sim[0])-max_radius):                
                if sim[x,y] == 1:
                    if (sim[y,x+tr] == 3 or sim[y,x-tr] == 3 \
                    or sim[y+tr,x] == 3 or sim[y-tr,x] == 3 \
                    or sim[y-dtr,x+dtr] == 3 or sim[y+dtr,x-dtr] == 3 \
                    or sim[y-dtr,x-dtr] == 3 or sim[y+dtr,x+dtr] == 3) \
                    and random() < InfectionRate:
                        newsim[x,y] = 3
                elif sim[x,y] == 3:
                    if random() < RecoveryRate:
                        newsim[x,y] = 5
                elif sim[x,y] == 5:
                    if random() < ImmunitylossRate:
                        newsim[x,y] = 1

        return newsim



    def grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        S = I = R = []
        sim = self.setup_grid(gridsize,max_radius)
        for j in range(0,simtime):
            sim = self.grid_frame(sim,max_radius)
            S.append(np.sum(sim == 1))
            I.append(np.sum(sim == 3))
            R.append(np.sum(sim == 5))

        return S,I,R

    def plot_grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        S,I,R = self.grid_sim(simtime, gridsize, max_radius)
        plt.plot(S, label = "Susceptible")
        plt.plot(I, label = "Infectious")
        plt.plot(R, label = "Recovered")
        plt.title("SIRS Pandemic Grid Model")
        plt.xlabel("Time(days)")
        plt.ylabel("Population")
        plt.legend()
        plt.show()
    
    def stackplot_grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        S,I,R = self.grid_sim(simtime, gridsize, max_radius)
        plt.stackplot(range(0,len(S)), S, I, R,  labels=["Susceptible","Infected","Recovered"])
        plt.title("SIRS Pandemic Grid Model")
        plt.xlabel("Time(days)")
        plt.ylabel("Population")
        plt.legend()
        plt.show()
    
    def print_grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        sim = self.setup_grid(gridsize, max_radius)
        for j in range(0,simtime):
            print()
            print(j)
            print(sim)
            sim = self.grid_frame(sim, max_radius)

    def animate_grid_sim(self,simtime=100,gridsize=100,max_radius=1):
        cmap = colors.ListedColormap(['steelblue','crimson','silver','black'])
        bounds = [0,2,4,6,8]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        sim = self.setup_grid(gridsize, max_radius)
        
        for n in range(0, simtime):
            plt.imshow(sim, cmap=cmap, norm=norm)
            plt.tick_params(bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
            plt.savefig("SIRS_gridframes/frame"+str(n)+".png")
            sim = self.grid_frame(sim, max_radius)

if __name__ == "__main__":
    SIRS = SIRS_Model()
    SIRS.InfectionsPlot()