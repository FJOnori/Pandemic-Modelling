from random import random, randint
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt

def setup_grid(gridsize=10):
        #concatinates rows and edges based on gridsize and max_radius
        sim = np.array([[1]*gridsize]*gridsize)
        #places a single infected individual in the middle of the grid
        centre = int(gridsize/2)
        sim[centre, centre] = 3
        return sim

def grid_frame_update(sim):
        
        size = len(sim[0])
        newsim = sim.copy()
        
        for x in range(size):
            for y in range(size):
                
                if sim[x,y] == 1:
                    if (sim[y,(x+1)%size] == 3 or sim[y,(x-1)%size] == 3 \
                    or sim[(y+1)%size,x] == 3 or sim[(y-1)%size,x] == 3 \
                    or sim[(y-1)%size,(x+1)%size] == 3 or sim[(y+1)%size,(x-1)%size] == 3 \
                    or sim[(y-1)%size,(x-1)%size] == 3 or sim[(y+1)%size,(x+1)%size] == 3) \
                    and random() < 0.5:
                        newsim[x,y] = 3
                elif sim[x,y] == 3:
                    if random() < 0.25:
                        newsim[x,y] = 5
                elif sim[x,y] == 5:
                    if random() < 0.015:
                        newsim[x,y] = 1
        
        return newsim

def grid_sim(simtime=1000, gridsize=50):
        sim = setup_grid(gridsize)
        cmap = colors.ListedColormap(['steelblue','crimson','silver'])
        bounds = [0,2,4,6]
        plt.tick_params(bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        for j in range(0,simtime):
            sim = grid_frame_update(sim)
            plt.cla()
            plt.imshow(sim, cmap=cmap)
            plt.draw()
            plt.pause(0.01)

            
           
        

grid_sim()