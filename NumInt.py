import numpy as np
import matplotlib.pyplot as plt


def SIR_NumInt(S_initial, I_inital, t_inital, t_final, dt, beta, gamma):
    
    S ,I ,R, T = [S_initial], [I_inital], [0], [t_inital]
    counter = 0
    population = S_initial + I_inital

    for t in np.arange(t_inital, t_final, dt):
        
        dS = (-beta*(I[counter]/population)*S[counter]) * dt
        dI = (beta*(I[counter]/population)*S[counter] - gamma*I[counter]) * dt
        dR = (gamma*I[counter]) * dt
        
        S.append(S[counter] + dS)
        I.append(I[counter] + dI)
        R.append(R[counter] + dR)
        T.append(t+dt)
        counter += 1

    plt.plot(T, S, label="Susceptible")
    plt.plot(T, I, label="Infected")
    plt.plot(T, R, label="Recovered")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title("SIR Pandemic model")
    plt.xlim(t_inital,t_final)
    plt.ylim(0, 1200)
    plt.legend()
    plt.show()

def SIS_NumInt(S_initial, I_inital, t_inital, t_final, dt, beta, gamma):
    
    S, I, T = [S_initial], [I_inital], [t_inital]
    counter = 0
    population = S_initial + I_inital

    for t in np.arange(t_inital, t_final, dt):
        
        dI = (beta*(I[counter]/population)*S[counter] - gamma*I[counter]) * dt
        dS = -dI

        S.append(S[counter] + dS)
        I.append(I[counter] + dI)
        T.append(t+dt)
        counter += 1

    plt.plot(T, S, label="Susceptible")
    plt.plot(T, I, label="Infected")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title("SIS Pandemic model")
    plt.xlim(t_inital,t_final)
    plt.ylim(0, 1200)
    plt.legend()
    plt.show()

def SIS_VectorField(S_initial, I_inital, beta, gamma):
    
    population = S_initial + I_inital
    x,y = np.meshgrid(np.linspace(1,population,25),np.linspace(1,population,25))
    dI = (beta*(x/population)*y - gamma*x)
    dS = -dI

    plt.quiver(x,y,dI,dS)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIS vector field")
    plt.xlim(0,population)
    plt.ylim(0,population)
    plt.show()

def SIR_VectorField(S_initial, I_inital, beta, gamma):
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    population = S_initial + I_inital
    x,y,z = np.meshgrid(np.linspace(1,population,5),np.linspace(1,population,5),np.linspace(1,population,5) )
    dS = (-beta*(y/population)*x)
    dI = (beta*(y/population)*x - gamma*y)
    dR = (gamma*y)

    ax.quiver(x,y,z,dS,dI,dR)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIR vector field")
    plt.xlim(0,population)
    plt.ylim(0,population)
    plt.show()



SIR_VectorField(S_initial = 1000, I_inital = 200, beta=0.3, gamma=0.1)

