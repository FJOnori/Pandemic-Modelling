import numpy as np
import matplotlib.pyplot as plt


class SIRS_vectorfields():

    def __init__(self) -> None:
        self.I_inital = 0.0001
        self.S_inital = 1 - self.I_inital
        self.beta = 0.05
        self.gamma = 0.01
        self.psi = 0.02
        self.sec = 0.5
        pass

    def SR(self, I=0.5):
    
        population = self.S_inital + self.I_inital
        x,y = np.meshgrid(np.linspace(0.001,1,40),np.linspace(0.001,1,40))
        dS = (-self.beta*(I/population)*x) + self.psi*y
        dR = self.gamma*I - self.psi*y

        color = 2 * np.log(np.hypot(dS, dR))
        plt.streamplot(x,y,dS,dR, cmap="Greys", color=color, density=1.5)
        plt.xlabel("Susceptible")
        plt.ylabel("Recovered")
        plt.title("SIRS Vector Field")
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.savefig("SIRS/SIRSvectorfieldSR.png", dpi=227)
        plt.close()


    def SI(self, R=0.5):
        
        population = self.S_inital + self.I_inital
        x,y = np.meshgrid(np.linspace(0.001,1,40),np.linspace(0.001,1,40))
        dS = (-self.beta*(y/population)*x) + self.psi*R
        dI = (self.beta*(y/population)*x) - self.gamma*y

        color = 2 * np.log(np.hypot(dS, dI))
        plt.streamplot(x,y,dS,dI, cmap="Greys", color=color, density=1.5)
        plt.xlabel("Susceptible")
        plt.ylabel("Infected")
        plt.title("SIRS Vector Field")
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.savefig("SIRS/SIRSvectorfieldSI.png", dpi=227)
        plt.close()

    def IR(self, S=0.5):
        
        population = self.S_inital + self.I_inital
        x,y = np.meshgrid(np.linspace(0.001,1,40),np.linspace(0.001,1,40))
        dI = (self.beta*(x/population)*S) - self.gamma*x
        dR = self.gamma*x - self.psi*y

        color = 2 * np.log(np.hypot(dI, dR))
        plt.streamplot(x,y,dI,dR, cmap="Greys", color=color, density=1.5)
        plt.xlabel("Infected")
        plt.ylabel("Recovered")
        plt.title("SIRS Vector Field")
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.savefig("SIRS/SIRSvectorfieldIR.png", dpi=227)
        plt.close()

    def SIR(self):

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        x, y, z = np.meshgrid(np.linspace(0,1,10),
                              np.linspace(0,1,10),
                              np.linspace(0,1,3))
        
        population = self.S_inital + self.I_inital

        dS = -self.beta*(y/population)*x + self.psi*z
        dI = self.beta*(y/population)*x - self.gamma*y 
        dR = self.gamma*y - self.psi*z

        ax.quiver(x,y,z,dS,dI,dR,length=0.1, normalize=True)
        ax.set_title("SIRS 3D Vectorfield")
        ax.set_xlabel('Susceptible')
        ax.set_ylabel('Infected')
        ax.set_zlabel('Recovered')
        plt.savefig("SIRS/SIRSvectorfieldSIR.png",dpi=227)
        plt.show()

    def runall(self):
        self.SI()
        self.SR()
        self.IR()
        self.SIR()


SIRS = SIRS_vectorfields()
SIRS.SIR()