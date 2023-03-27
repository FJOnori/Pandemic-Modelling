import random
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt


def setup_grid_centre(gridsize=50):
    #concatinates rows and edges based on gridsize and max_radius
    
    sim = []
    for n in range(0, gridsize):
            sim.append(random.choices([1,2], weights=(0.9,0.1), k=gridsize))
    sim = np.array(sim)
    #places a single infected individual in the middle of the grid
    centre = int(gridsize/2)
    sim[centre, centre] = 9
    return sim
    
def progress_bar(current, total, bar_length=100):
        fraction = current / total
        arrow = int(fraction * bar_length - 1) * '-' + '>'
        padding = int(bar_length - len(arrow)) * ' '
        ending = '\n' if current == total else '\r'
        print(f'Progress: [{arrow}{padding}] {round((fraction*100),5)}%', end=ending)

def grid_frame_update(sim):
        
    size = len(sim[0])
    newsim = sim.copy()

    ContactRate = 1/5
    LatencyRate = 1/6
    RecoveryRate = 1/14
    ImmunityLossRate = 1/210
    VaccinationRate1 = 1/365
    VaccinationRate2 = 1/50

    
    for x in range(size):
        for y in range(size):

            neighbours = [sim[y,(x+1)%size],sim[y,(x-1)%size], sim[(y+1)%size,x], sim[(y-1)%size,x], \
                          sim[(y-1)%size,(x+1)%size], sim[(y+1)%size,(x-1)%size], sim[(y-1)%size,(x-1)%size], sim[(y+1)%size,(x+1)%size]]
            contagious = [3,9]
            badsquare = any(x in neighbours for x in contagious)


            #Susceptible
            if sim[x,y] == 1:
                if badsquare: newsim[x,y] = random.choices([1,3], weights=(1 - ContactRate, ContactRate))[0]   
                else: newsim[x,y] = random.choices([1,21], weights=(1 - VaccinationRate1,VaccinationRate1))[0]   

            #Exposed
            elif sim[x,y] == 3: 
                newsim[x,y] = random.choices([3,9], weights=(1 - LatencyRate, LatencyRate))[0]

            #Infected
            elif sim[x,y] == 9:
                newsim[x,y] = random.choices([9,15], weights=(1 - RecoveryRate, RecoveryRate))[0]

            #Recovered
            elif sim[x,y] == 15:
                newsim[x,y] = random.choices([15,1,16], weights=(1 - ImmunityLossRate - VaccinationRate1,ImmunityLossRate,VaccinationRate1))[0]

            #Vaccinated
            elif sim[x,y] == 21:
                if badsquare: newsim[x,y] = random.choices([21,4], weights=(1 - (ContactRate*0.19), (ContactRate*0.19)))[0]
                else: newsim[x,y] = random.choices([21,23], weights=(1 - VaccinationRate2 , VaccinationRate2))[0]
            
            #Vaccinated Infected
            elif sim[x,y] == 10:
                newsim[x,y] = random.choices([10,16], weights=(1 - RecoveryRate, RecoveryRate))[0]

            #Vaccinated Exposed
            elif sim[x,y] == 4:
                newsim[x,y] = random.choices([4,10], weights=(1 - LatencyRate, LatencyRate))[0]

            #Vaccinated Recovered
            elif sim[x,y] == 16:
                newsim[x,y] = random.choices([16,21,17], weights=(1 - ImmunityLossRate - VaccinationRate2,ImmunityLossRate,VaccinationRate2))[0]
            
            #Fully Vaccinated
            elif sim[x,y] == 23:
                if badsquare: newsim[x,y] = random.choices([23,5], weights=(1 - (ContactRate*0.09), (ContactRate*0.09)))[0]
            
            #Fully Vaccinated Exposed
            elif sim[x,y] == 5:
                newsim[x,y] = random.choices([5,11], weights=(1 - LatencyRate, LatencyRate))[0]

            #Fully Vaccinated Infected
            elif sim[x,y] == 11:
                newsim[x,y] = random.choices([11,17], weights=(1 - RecoveryRate, RecoveryRate))[0]

            #Fully Vaccinated Recovered
            elif sim[x,y] == 17:
                newsim[x,y] = random.choices([17,23], weights=(1-ImmunityLossRate, ImmunityLossRate))[0]

            #Susceptible High Risk
            elif sim[x,y] == 2:
                if badsquare: newsim[x,y] = random.choices([2,6], weights=(1 - (ContactRate*2), (ContactRate*2)))[0]
                else: newsim[x,y] = random.choices([2,22], weights=(1-VaccinationRate1,VaccinationRate1))[0]

            #Exposed Hish Risk
            elif sim[x,y] == 6:
                newsim[x,y] = random.choices([6,12], weights=(1 - LatencyRate, LatencyRate))[0]
            
            #Infected Hish Risk
            elif sim[x,y] == 12:
                newsim[x,y] = random.choices([12,18], weights=(1 - RecoveryRate, RecoveryRate))[0]
            
            #Recovered High Risk
            elif sim[x,y] == 18:
                newsim[x,y] = random.choices([18,2,19], weights=(1 - ImmunityLossRate - VaccinationRate1,ImmunityLossRate,VaccinationRate1))[0]

            #Vaccinated High Risk
            elif sim[x,y] == 22:
                if badsquare: newsim[x,y] = random.choices([22,7], weights=(1 - (ContactRate*0.19*2), (ContactRate*0.19*2)))[0]
                else: newsim[x,y] = random.choices([21,24], weights=(1 - VaccinationRate2 , VaccinationRate2))[0]
            
            #Exposed Vaccinated High Risk
            elif sim[x,y] == 7:
                newsim[x,y] = random.choices([7, 13], weights=(1 - LatencyRate, LatencyRate))[0]

            #Infected Vaccinated High Risk
            elif sim[x,y] == 13:
                newsim[x,y] = random.choices([13, 19], weights=(1 - RecoveryRate, RecoveryRate))[0]

            #Recovered Vaccinated High Risk
            elif sim[x,y] == 19:
                newsim[x,y] = random.choices([19,22,20], weights=(1 - ImmunityLossRate - VaccinationRate2,ImmunityLossRate,VaccinationRate2))[0]
            
            #Fully Vaccinated High Risk
            elif sim[x,y] == 24:
                newsim[x,y] = random.choices([24, 8], weights=(1 - (ContactRate*2*0.09), (ContactRate*2*0.09)))[0]
            
            #Fully Vaccinated Exposed High Risk
            elif sim[x,y] == 8:
                newsim[x,y] = random.choices([8, 14], weights=(1 - LatencyRate, LatencyRate))[0]

            #Fully Vaccinated Infected High Risk
            elif sim[x,y] == 14:
                newsim[x,y] = random.choices([14, 20], weights=(1 - RecoveryRate, RecoveryRate))[0]
            
            #Fully Vaccinated Recovered High Risk
            elif sim[x,y] == 20:
                newsim[x,y] = random.choices([20, 24], weights=(1-ImmunityLossRate, ImmunityLossRate))[0]

    return newsim

sim = setup_grid_centre(gridsize=100)
S,I,R,V1,V2 = [],[],[],[],[]

for n in range(1000):
        progress_bar(n,1000)
        sim = grid_frame_update(sim)

        S.append(np.count_nonzero(sim == 1) + np.count_nonzero(sim == 2))
        I.append(np.count_nonzero(sim == 3) + np.count_nonzero(sim == 4) + np.count_nonzero(sim == 5) + np.count_nonzero(sim == 6) + np.count_nonzero(sim == 7) + np.count_nonzero(sim == 8) + \
                 + np.count_nonzero(sim == 9) + np.count_nonzero(sim == 10) + np.count_nonzero(sim == 11) + np.count_nonzero(sim == 12) + np.count_nonzero(sim == 13) + np.count_nonzero(sim == 14) )
        
        R.append(np.count_nonzero(sim == 15) + np.count_nonzero(sim == 16) + np.count_nonzero(sim == 17) + np.count_nonzero(sim == 18) + np.count_nonzero(sim == 19) + np.count_nonzero(sim == 20))
        V1.append(np.count_nonzero(sim == 21) + np.count_nonzero(sim == 22))
        V2.append(np.count_nonzero(sim == 23) + np.count_nonzero(sim == 24))

cmap = colors.ListedColormap(colors=['#1e68b3','#b81111','#aaacad','#5e017d','#439603'])
#plt.plot(S, label ="Susceptible", c='#1e68b3')
#plt.plot(I, label="Infectious", c='#b81111')
#plt.plot(R, label="Recovered", c='#aaacad')
#plt.plot(V1, label="Vaccinated", c='#5e017d')
#plt.plot(V2, label="Fully Vaccinated", c='#439603')
#plt.title("COVIDGE grid simulation")
#plt.xlabel("Time (days)")#
#plt.ylabel("Population")
#plt.legend()
#plt.savefig("COVIDGE_gridsim.png", dpi=227)

#plt.close()

plt.title("COVIDGE Model Stackplot")
plt.xlabel("Time (Days)")
plt.ylabel("Population")
plt.stackplot(1000 , S,I,R,V1,V2, labels=["Susceptible","Infectious","Recovered","Vaccinated","Fully Vaccinated"], colors=['#1e68b3','#b81111','#aaacad','#5e017d','#439603'])
plt.legend()
plt.savefig("COVIDGEstackplot.png", dpi=227)