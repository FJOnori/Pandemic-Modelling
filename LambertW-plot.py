import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw


N = 1
gamma = 0.25
beta = 0.5
s0 = N - 0.00001



Wsub = -(beta/(N*gamma))*s0*np.exp(-beta/gamma)
Wmul = -((N*gamma)/beta)

print(Wmul*lambertw(Wsub))

w = []
for n in np.arange(0,300,0.01):
    w.append(lambertw(n))

plt.plot(w)
plt.title("Lambert-W function")
plt.xlabel("x-input")
plt.ylabel("y-output")
plt.ylim(0,1.2)
plt.xlim(0,300)
plt.savefig("W.png",dpi=227)
plt.grid(ls=":",c='grey')
plt.show()

