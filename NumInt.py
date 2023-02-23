import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def SIR_NumInt(S_initial, I_inital, t_inital, t_final, dt, beta, gamma, nu):
    
    S ,I ,R, T = [S_initial], [I_inital], [0], [t_inital]
    counter = 0
    population = S_initial + I_inital

    for t in np.arange(t_inital, t_final, dt):
        
        dS = (-beta*(I[counter]/population)*S[counter] + nu*R[counter]) * dt
        dI = (beta*(I[counter]/population)*S[counter] - gamma*I[counter]) * dt
        dR = (gamma*I[counter] - nu*R[counter]) * dt
        
        S.append(S[counter] + dS)
        I.append(I[counter] + dI)
        R.append(R[counter] + dR)
        T.append(t+dt)
        counter += 1

    plt.plot(T, np.array(S)*100, label="Susceptible", c="#1d6bc4")
    plt.plot(T, np.array(I)*100, label="Infected", c="#c41d1d")
    plt.plot(T, np.array(R)*100, label="Recovered", c="#757575")
    plt.xlabel("Time (days)")
    plt.ylabel("Percentage of Population")
    plt.title("SIRS Pandemic model")
    plt.xlim(t_inital,t_final)
    plt.ylim(0, 100)
    plt.legend()
    plt.savefig("SIRS2.png", dpi=227)

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

def SEIRS_NumInt(S_initial, E_inital, t_inital, t_final, dt, beta, gamma, mu, nu, epsilon, sigma):
    
    S ,I, E ,R, T = [S_initial], [0], [E_inital], [0], [t_inital]
    counter = 0
    population = S_initial + E_inital

    for t in np.arange(t_inital, t_final, dt):
        
        dS = (mu*population - beta*((I[counter]+ E[counter])/population)*S[counter] + epsilon*R[counter] - nu*S[counter]) * dt
        dE = (beta*((I[counter] + E[counter])/population)*S[counter] - sigma*E[counter] - nu*E[counter]) * dt
        dI = (sigma*E[counter] - gamma*I[counter] - nu*I[counter]) * dt
        dR = (gamma*I[counter] - epsilon*R[counter] - nu*R[counter]) * dt
        
        S.append(S[counter] + dS)
        E.append(E[counter] + dE)
        I.append(I[counter] + dI)
        R.append(R[counter] + dR)
        T.append(t+dt)
        counter += 1

    plt.plot(T, S, label="Susceptible")
    plt.plot(T, E, label="Exposed")
    plt.plot(T, I, label="Infected")
    plt.plot(T, R, label="Recovered")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title("SEIRS Pandemic model")
    plt.legend()
    plt.show()

def SEVIRS_NumInt(S_initial, E_inital, t_inital, t_final, dt, beta, gamma, mu, nu, epsilon, sigma, upsilon1, upsilon2):
    
    S ,I, E ,R, V, T = [S_initial], [0], [E_inital], [0], [0], [t_inital]
    counter = 0
    population = S_initial + E_inital

    for t in np.arange(t_inital, t_final, dt):
        
        dS = (mu*population - beta*((I[counter]+ E[counter])/population)*S[counter] + epsilon*R[counter] - nu*S[counter] - upsilon1*S[counter]) * dt
        dE = (beta*((I[counter] + E[counter])/population)*S[counter] - sigma*E[counter] - nu*E[counter]) * dt
        dI = (sigma*E[counter] - gamma*I[counter] - nu*I[counter]) * dt
        dR = (gamma*I[counter] - epsilon*R[counter] - nu*R[counter] - upsilon2*R[counter]) * dt
        dV = (upsilon1*S[counter] + upsilon2*R[counter] - nu*V[counter]) * dt
        
        S.append(S[counter] + dS)
        E.append(E[counter] + dE)
        I.append(I[counter] + dI)
        R.append(R[counter] + dR)
        V.append(V[counter] + dV)
        T.append(t+dt)
        counter += 1

    plt.plot(T, S, label="Susceptible")
    plt.plot(T, E, label="Exposed")
    plt.plot(T, I, label="Infected")
    plt.plot(T, R, label="Recovered")
    plt.plot(T, V, label="Vaccinated")
    plt.xlabel("Time (Days)")
    plt.ylabel("Population (Percentage)")
    plt.title("SEVIRS Pandemic model")
    plt.legend()
    plt.show()


def SIS_VectorField(S_initial, I_inital, beta, gamma):
    
    population = S_initial + I_inital
    x,y = np.meshgrid(np.linspace(1,population,25),np.linspace(1,population,25))
    dI = (beta*(x/population)*y - gamma*x)
    dS = -dI

    plt.quiver(y,x,dS,dI)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIS vector field")
    plt.show()

def SIR_VectorField(S_initial, I_inital, beta, gamma):
    population = S_initial + I_inital
    x,y = np.meshgrid(np.linspace(1,population,5),np.linspace(1,population,5))
    dS = (-beta*(y/population)*x)
    dI = (beta*(y/population)*x - gamma*y)

    print(y)

    plt.quiver(x,y,dS,dI)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIR vector field")
    plt.xlim(0,population)
    plt.ylim(0,population)
    plt.show()

def SIRS_VectorField(S_initial, I_inital, beta, gamma):
    population = S_initial + I_inital
    x,y = np.meshgrid(np.linspace(1,population,5),np.linspace(1,population,5))
    dS = (-beta*(y/population)*x)
    dI = (beta*(y/population)*x - gamma*y)

    print(y)

    plt.quiver(x,y,dS,dI)
    plt.xlabel("Susceptible")
    plt.ylabel("Infected")
    plt.title("SIR vector field")
    plt.xlim(0,population)
    plt.ylim(0,population)
    plt.show()

def beta_graph_R():
    df = pd.read_csv('beta.csv')
    beta = list(df['beta'])
    time = np.arange(0, len(beta))
    plt.plot(time, beta)
    plt.xlim(0,len(beta))
    plt.ylim(0,0.25)
    plt.xlabel("Time (days)")
    plt.ylabel("BETA")
    plt.title("UK Contact rate")
    plt.show()

def beta_calculate(S_initial, E_inital, t_inital, t_final, dt, beta, mu, nu, epsilon, upsilon1):
    df = pd.read_csv('owid-covid-data-uk.csv')
    S ,I, E ,R, V, T = [S_initial], [0], [E_inital], [0], [0], [t_inital]
    counter = 0
    population = S_initial + E_inital

def vaccination2_rate_csv():
    df = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
    vrate = []
    fullyVaccinated = list(df['people_fully_vaccinated'])

    for n in range(0, len(df['people_fully_vaccinated'])-1):
        vrate.append((fullyVaccinated[n+1] - fullyVaccinated[n])/67508935)
    
    v = {"vaccination_rate" : vrate}
    tocsv = pd.DataFrame(v)
    tocsv.to_csv('Vaccine2rate.csv')

def vaccination_rate_plot():
    df = (pd.read_csv('owid-covid-data-uk.csv')).replace(np.nan, 0)
    plt.plot(list(df['people_vaccinated']))
    plt.plot(list(df['people_fully_vaccinated']))
    plt.show()

def R_number():
    df = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)
    plt.plot(np.array(df['reproduction_rate'])*1/4)
    
def point():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(0.1, 0.9, 0.2, c='k')

    ax.set_title("Initial State of System")
    ax.set_xlabel('Susceptible')
    ax.set_ylabel('Infected')
    ax.set_zlabel('Recovered')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    plt.savefig("ReportImgs/PointSIR_SIRS.png",dpi=227)
    plt.show()

#SIR_NumInt(S_initial=1-0.0001, I_inital=0.0001, t_inital=0, t_final=1000, dt=1, beta=0.1, gamma=0.05, nu=0.01)
point()
