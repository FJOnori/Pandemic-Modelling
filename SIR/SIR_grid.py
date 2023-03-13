from random import random, randint
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt

def setup_grid_centre(gridsize=10):
        #concatinates rows and edges based on gridsize and max_radius
        sim = np.array([[1]*gridsize]*gridsize)
        #places a single infected individual in the middle of the grid
        centre = int(gridsize/2)
        sim[centre, centre] = 3
        return sim

def setup_grid_random(gridsize=10):
    sim = np.array([[1]*gridsize]*gridsize)
    for x in range(gridsize):
        for y in range(gridsize):
            if random() < 0.2:
                sim[x,y] = 3
    return sim

def grid_frame_update(sim):
        
        size = len(sim[0])
        newsim = sim.copy()
        
        for x in range(size):
            for y in range(size):

                n = random()
                
                if sim[x,y] == 1:

                    if (sim[y,(x+1)%size] == 3 or sim[y,(x-1)%size] == 3 \
                    or sim[(y+1)%size,x] == 3 or sim[(y-1)%size,x] == 3 \
                    or sim[(y-1)%size,(x+1)%size] == 3 or sim[(y+1)%size,(x-1)%size] == 3 \
                    or sim[(y-1)%size,(x-1)%size] == 3 or sim[(y+1)%size,(x+1)%size] == 3):

                        nearsq = [sim[y,(x+1)%size], sim[y,(x-1)%size], 
                                  sim[(y+1)%size,x], sim[(y-1)%size,x],
                                  sim[(y-1)%size,(x+1)%size], sim[(y+1)%size,(x-1)%size],
                                  sim[(y-1)%size,(x-1)%size], sim[(y+1)%size,(x+1)%size]]

                        infecProb = (nearsq.count(3))/8

                        if n < infecProb:
                            newsim[x,y] = 3

                elif sim[x,y] == 3:
                    if n < 0.1:
                        newsim[x,y] = 5
        
        return newsim

def grid_sim(simtime=200, gridsize=100):
        sim = setup_grid_centre(gridsize) #blue, #red, #grey
        cmap = colors.ListedColormap(['#1e68b3','#b81111','#aaacad'])
        bounds = [0,2,4,6]
        plt.tick_params(bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)


        S = []
        I = []
        R = []

        for j in range(0,simtime):
            sim = grid_frame_update(sim)
            
            S.append(np.count_nonzero(sim == 1))
            I.append(np.count_nonzero(sim == 3))
            R.append(np.count_nonzero(sim == 5))

            plt.cla()
            plt.imshow(sim, cmap=cmap)
            plt.draw()
            plt.pause(0.00001)
        
        plt.close()

        plt.plot(S, c='#1e68b3', label="Susceptible")
        plt.plot(I, c='#b81111', label="Infected")
        plt.plot(R, c='#aaacad', label="Recovered")
        plt.title("SIR Grid Simulation")
        plt.ylabel("Population Proportion")
        plt.xlabel("Time (days)")
        plt.xlim(0,simtime)
        plt.ylim(0,gridsize**2)
        plt.legend()
        plt.savefig("SIR/SIRgridsim.png", dpi=227)
           
grid_sim()