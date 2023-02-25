import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
from random import random, randint

class SIRS_Model():

    def __init__(self) -> None:
        COVIDdf = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)
        
        self.T_final = 365*5
        self.S_initial = 1 - 0.0000056
        self.I_inital = 0.0000056
        self.population = self.S_initial + self.I_inital
       
        self.BirthRate = 2.78363e-5
        self.DeathRate = 2.48321e-5
        self.COVIDDeathRate = 8.95027e-4
        self.ImmunityLossRate = 1/210
        self.RecoveryRate = 1/10
        self.InfectionRate = 0.19
        self.TDInfectionRate = (self.RecoveryRate)*np.array(COVIDdf['reproduction_rate'][30:1030])

        self.InfectionRateGrid = 0.6
        self.RecoveryRateGrid = 0.1
        self.ImmunitylossRateGrid = 0.05

    def NumInt(self, TD=True):

        S ,I, R, D = [self.S_initial], [self.I_inital], [0], [0]

        if TD: beta = self.TDInfectionRate 
        else: beta = np.array([self.InfectionRate]*len(self.TDInfectionRate))

        for n in np.arange(0, len(self.TDInfectionRate)):

               dS = (self.BirthRate*self.population) - beta[n]*(I[n]/self.population)*S[n] \
                + self.ImmunityLossRate*R[n] - self.DeathRate*S[n]

               dI = beta[n]*(I[n]/self.population)*S[n] - self.RecoveryRate*I[n] - self.COVIDDeathRate*I[n] 

               dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n]

               self.population += dI + dR + dS

               S.append(S[n] + dS)
               I.append(I[n] + dI)
               R.append(R[n] + dR)
               D.append(self.COVIDDeathRate*I[n])

        return np.array(S), np.array(I), np.array(R), np.array(D)

    def LinePlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(np.array(S)*100, label="Susceptible", c="#1d6bc4")
        plt.plot(np.array(I)*100, label="Infected", c="#c41d1d")
        plt.plot(np.array(R)*100, label="Recovered", c="#757575")
        plt.xlabel("Time (Days)")
        plt.ylabel("Percentage of Population")
        plt.title("SIRS Pandemic Model")
        plt.ylim(0,100)
        plt.xlim(0,len(self.TDInfectionRate))
        plt.legend()
        plt.show()

    def StackPlot(self):
        S ,I, R, D = self.NumInt()
        plt.stackplot(range(0,len(S)), S, I, R, labels=["Susceptible","Infected","Recovered"], colors=["#1d6bc4","#c41d1d","#757575"])
        plt.xlabel("Time (Days)")
        plt.ylabel("Population")
        plt.xlim(0,self.T_final)
        plt.ylim(0,self.population)
        plt.title("SIRS Pandemic model StackPlot")
        plt.legend()
        plt.show()

    def InfectionsPlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(I*67000000, label="Infections", c='#a80808')
        plt.xlabel("Time (Days)")
        plt.ylabel("Infections from COVID-19")
        plt.xlim(0,(len(self.TDInfectionRate)))
        plt.ylim(0,120000)
        plt.title("SIRS Pandemic model Infections")
        plt.grid(ls=':',c='grey', axis='y')
        plt.savefig("ReportImgs/FinalResults/SIRS_infectioncurve.png", dpi=227)
        plt.show()

    def DeathsPlot(self):
        S ,I, R, D = self.NumInt()
        plt.plot(D, label="Deaths", c='k')
        plt.xlabel("Time (Days)")
        plt.ylabel("Deaths from COVID-19")
        plt.title("SIRS Deaths")
        plt.show()

    def PieChart(self):
        S ,I, R, D = self.NumInt()
        data = np.array([S[-1], I[-1], R[-1]])
        l = ["Susceptible", "Infected", "Recovered"]
        plt.pie(data, labels=l, colors=["#1d6bc4","#c41d1d","#757575"])
        plt.title("SIRS Pandemic Final State")
        plt.savefig("ReportImgs/SIRS_PieChart",dpi=227)

    def VectorField(self, axis="SIR", time= 400, spacing = 8):
        
        x, y = np.meshgrid(np.linspace(0,self.population,spacing),
                           np.linspace(0,self.population,spacing))

        S ,I, R, D = self.NumInt()

        if axis == "SR":
            dS = (self.BirthRate*self.population) - self.InfectionRate*(I[time]/self.population)*x + self.ImmunityLossRate*y- self.DeathRate*x
            dR = self.RecoveryRate*I[time] - self.ImmunityLossRate*y - self.DeathRate*y
            color = 2 * np.log(np.hypot(dS, dR))
            plt.streamplot(x,y,dS,dR,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.xlabel("Susceptible")
            plt.ylabel("Recovered")
            plt.title("SIRS Vectorfield")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.savefig("ReportImgs/VectorfieldSR_SIRS.png",dpi=227)
        elif axis == "IR":
            dI = self.InfectionRate*(x/self.population)*S[time] - self.RecoveryRate*x - self.COVIDDeathRate*x 
            dR = self.RecoveryRate*x - self.ImmunityLossRate*y - self.DeathRate*y
            color = 2 * np.log(np.hypot(dI, dR))
            plt.streamplot(x,y,dI,dR,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.title("SIRS Vectorfield")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xlabel("Infected")
            plt.ylabel("Recovered")
            plt.savefig("ReportImgs/VectorfieldIR_SIRS.png",dpi=227)
        elif axis == "IS":
            dS = (self.BirthRate*self.population) - self.InfectionRate*(x/self.population)*y + self.ImmunityLossRate*R[time] - self.DeathRate*y
            dI = self.InfectionRate*(x/self.population)*y - self.RecoveryRate*x - self.COVIDDeathRate*x
            color = 2 * np.log(np.hypot(dS, dI))
            plt.streamplot(x,y,dS,dI,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.title("SIRS Vectorfield")
            plt.xlabel("Infected")
            plt.ylabel("Susceptible")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.savefig("ReportImgs/VectorfieldIS_SIRS.png",dpi=227)
        
        elif axis == "SIR":

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            x, y, z = np.meshgrid(np.linspace(0,self.population,spacing),
                                  np.linspace(0,self.population,spacing),
                                  np.linspace(0,self.population,spacing))

            dS = (self.BirthRate*self.population) - self.InfectionRate*(y/self.population)*x + self.ImmunityLossRate*z - self.DeathRate*x
            dI = self.InfectionRate*(y/self.population)*x - self.RecoveryRate*y - self.COVIDDeathRate*y
            dR = self.RecoveryRate*y - self.ImmunityLossRate*z - self.DeathRate*z

            ax.quiver(x,y,z,dS,dI,dR)
            ax.set_title("Proportion of initial population in each compartment")
            ax.set_xlabel('Susceptible')
            ax.set_ylabel('Infected')
            ax.set_zlabel('Recovered')
            plt.savefig("ReportImgs/VectorfieldSIR_SIRS.png",dpi=227)
            plt.show()

    def constant_population_vectorfield(self, spacing=25, axis="IR"):

        x, y = np.meshgrid(np.linspace(0,1,spacing),
                           np.linspace(0,1,spacing))
        
        if axis == "IS":
            dS = - self.InfectionRate*(y/self.population)*x + self.ImmunityLossRate*(1 - x - y)
            dI = self.InfectionRate*(y/self.population)*x - self.RecoveryRate*y
            color = 2 * np.log(np.hypot(dS, dI))
            plt.streamplot(x,y,dS,dI,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.title("SIRS Vectorfield")
            plt.xlabel("Infected")
            plt.ylabel("Susceptible")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.show()
        elif axis == "SR":
            dS = -self.InfectionRate*((1-x-y)/self.population)*x + self.ImmunityLossRate*y
            dR = self.RecoveryRate*(1-x-y) - self.ImmunityLossRate*y
            color = 2 * np.log(np.hypot(dS, dR))
            plt.streamplot(x,y,dS,dR,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.title("SIRS Vectorfield")
            plt.xlabel("Susceptible")
            plt.ylabel("Recovered")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.show()
        elif axis == "IR":
            dI = self.InfectionRate*(x/self.population)*(1-x-y)- self.RecoveryRate*x - self.COVIDDeathRate*x 
            dR = self.RecoveryRate*x - self.ImmunityLossRate*y - self.DeathRate*y
            color = 2 * np.log(np.hypot(dR, dI))
            plt.streamplot(x,y,dR,dI,color=color, linewidth=1,cmap="Blues",density=1.5, arrowstyle='->', arrowsize=1.5)
            plt.title("SIRS Vectorfield")
            plt.xlabel("Infected")
            plt.ylabel("Susceptible")
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.show()
            pass

    def setup_grid(self, gridsize=10, max_radius=1):
        #defines a standard row of the grid
        row = [7]*max_radius + [1]*gridsize + [7]*max_radius
        #defines a top and bottom edge of the grid
        edge = [7]*(gridsize+(max_radius*2))
        #concatinates rows and edges based on gridsize and max_radius
        sim = np.array([edge]*max_radius + [row]*gridsize + [edge]*max_radius)
        #places a single infected individual in the middle of the grid
        centre = int((gridsize+(max_radius*2))/2)
        sim[centre, centre] = 3

        return sim

    def grid_frame_update(self, sim):
        
        max_radius = list(sim).count(list(sim)[0])/2
        print(max_radius)
        tr = randint(1,max_radius)
        dtr = int(tr * (1/np.sqrt(2)))

        newsim = sim.copy()
        for x in range(len(sim[0])-max_radius):
            for y in range(len(sim[0])-max_radius):
                if sim[x,y] == 7:
                    pass                
                elif sim[x,y] == 1:
                    if (sim[y,x+tr] == 3 or sim[y,x-tr] == 3 \
                    or sim[y+tr,x] == 3 or sim[y-tr,x] == 3 \
                    or sim[y-dtr,x+dtr] == 3 or sim[y+dtr,x-dtr] == 3 \
                    or sim[y-dtr,x-dtr] == 3 or sim[y+dtr,x+dtr] == 3) \
                    and random() < self.InfectionRateGrid:
                        newsim[x,y] = 3
                elif sim[x,y] == 3:
                    if random() < self.RecoveryRateGrid:
                        newsim[x,y] = 5
                elif sim[x,y] == 5:
                    if random() < self.ImmunitylossRateGrid:
                        newsim[x,y] = 1
        
        return newsim

    def grid_sim(self, simtime=10, gridsize=10, max_radius=1):
        S = I = R = []
        sim = self.setup_grid(gridsize,max_radius)
        for j in range(0,simtime):
            sim = self.grid_frame_update(sim)
            print(j)
            print(sim)
            S.append(np.sum(sim == 1))
            I.append(np.sum(sim == 3))
            R.append(np.sum(sim == 5))

        return S,I,R,sim

    def calculate_infection_speed_grid(self, simtime=100, gridsize=10, max_radius=1):

        S,I,R,sim = self.grid_sim(simtime, gridsize, max_radius)
        startpos = [int(gridsize/2), int(gridsize/2)]
        furthest_infection = startpos

        for x in range(len(sim[0])-max_radius):
            for y in range(len(sim[0])-max_radius):
                if sim[x,y] == 3 and ((x-startpos[0])**2 + (y-startpos[1])**2) > np.linalg.norm(furthest_infection):
                    furthest_infection = np.array([x,y])

        infection_speed = np.linalg.norm(furthest_infection-startpos)/simtime
        return infection_speed

    def plot_grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        S,I,R,sim = self.grid_sim(simtime, gridsize, max_radius)
        plt.plot(S, label = "Susceptible")
        plt.plot(I, label = "Infectious")
        plt.plot(R, label = "Recovered")
        plt.title("SIRS Pandemic Grid Model")
        plt.xlabel("Time(days)")
        plt.ylabel("Population")
        plt.legend()
        plt.show()
    
    def stackplot_grid_sim(self, simtime=100, gridsize=10, max_radius=1):
        S,I,R,sim = self.grid_sim(simtime, gridsize, max_radius)
        plt.stackplot(range(0,len(S)), S, I, R,  labels=["Susceptible","Infected","Recovered"])
        plt.title("SIRS Pandemic Grid Model")
        plt.xlabel("Time(days)")
        plt.ylabel("Population")
        plt.legend()
        plt.show()
    
    def animate_grid_sim(self,simtime=100,gridsize=100,max_radius=1):
        cmap = colors.ListedColormap(['steelblue','crimson','silver','black'])
        bounds = [0,2,4,6,8]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        sim = self.setup_grid(gridsize, max_radius)
        updateInterval = 1


        fig, ax = plt.subplots()
        img = ax.imshow(sim, interpolation='nearest', cmap=cmap)
        ani = animation.FuncAnimation(fig, self.grid_frame_update, fargs=(sim, gridsize, ),
								frames = 100,
								interval=updateInterval,
								save_count=50)

        
        plt.tick_params(bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
        plt.show()
        
if __name__ == "__main__":
    SIRS = SIRS_Model()
    SIRS.grid_sim()

    